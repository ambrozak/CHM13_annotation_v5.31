#!/usr/bin/env python3
"""
Compare protein sequences from three sources for each discrepant tmap entry:

  1) protein encoded on GRCh38 by the original MANE transcript
     (MANEv1.5.gff + hg38.p14.fa)
  2) protein encoded on chm13 by the lifted-over MANE transcript
     (fixed_MANE_on_chm13.gff + chm13v2.0.fa)
  3) protein encoded on chm13 by the existing annotation's closest match
     (chm13v2.0_RefSeq_Liftoff_v5.3.gff + chm13v2.0.fa)

For each tmap row with class_code != "=", computes %identity(1,2) and
%identity(1,3) and tallies whether the existing annotation, the lifted-over
annotation, or neither is closer to the true MANE protein.

The tmap is built by gffcompare against a gffread --gene2exon version of
the existing annotation, in which every transcript-less gene (typically a
pseudogene) becomes its own single-exon transcript, so gffcompare itself
detects pseudogene overlap. That file is stripped of attributes, though,
so this script reads the full, original existing GFF (EXISTING_GFF) for
CDS, Dbxref GeneID, pseudo=true and transl_except. A row whose matched
ref_id is flagged pseudo=true there is kept but reported in its own
separate category (a two-way MANE-vs-lifted comparison, same treatment as
a genuinely novel locus) rather than folded into the normal three-way
comparison.

transl_except attributes on CDS lines (selenocysteine-recoded UGA, non-AUG
starts read as Met, etc.) are applied during translation, so a Sec codon
no longer truncates the protein at an apparent premature stop.

Requires: pandas, pyfaidx, biopython, tqdm
    pip install pyfaidx biopython tqdm
"""

import csv
import re
import sys
from collections import defaultdict
from urllib.parse import unquote

import pandas as pd
from pyfaidx import Fasta
from Bio.Seq import Seq
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.Data.IUPACData import protein_letters_3to1_extended
from tqdm import tqdm

# ----------------------------------------------------------------------
# Paths (edit if needed)
# ----------------------------------------------------------------------
TMAP_PATH = "gffcmp.MANE_on_chm13.gff.tmap"
MANE_GFF = "MANEv1.5.gff"
MANE_FASTA = "hg38.p14.fa"
LIFTED_GFF = "MANE_on_chm13.gff"
EXISTING_GFF = "chm13v5.3.gff"
CHM13_FASTA = "chm13v2.0.fa"

OUT_CSV = "protein_comparison_results.csv"
ERROR_LOG = "error.log"

# Optional: one NCBI GeneID per line. If the file exists, every lifted copy
# and transcript of those genes is traced through the run into GENE_TRACE.
GENE_LIST_PATH = "gene_ids_new_in_MANE1.5.txt"
GENE_TRACE = "gene_trace.tsv"

# Add/Review triage for three-way rows where the lifted-over transcript is
# closer to MANE: "Add" requires the lifted and existing transcripts to share
# an NCBI GeneID AND the existing protein to have at least this identity to
# MANE (inclusive). Anything else is "Review".
ADD_MIN_EXISTING_IDENTITY = 0.80

# Add/Review for novel-locus rows (ref_id == '-'): Add requires MANE-vs-lifted
# protein identity >= this (inclusive), and no LiftOn mutation containing any
# of these terms (case-insensitive substring, e.g. 'stop_codon_gain').
# Pseudogene-match rows are always Review.
NOVEL_ADD_MIN_IDENTITY = 0.95
NOVEL_BLOCKING_MUTATION_TERMS = ("codon", "frameshift")

PSEUDO_ATTR_RE = re.compile(r"(?:^|;)pseudo=true(?:;|$)")

# ----------------------------------------------------------------------
# GFF parsing: build transcript_id -> {chrom, strand, cds:[(start,end)]}
# for CDS-bearing transcripts, plus a set of every transcript-like ID in
# the file (coding or not), so a failed CDS lookup can be classified as
# "non-coding transcript" rather than lumped in with genuinely absent IDs.
# NOTE: assumes each CDS is 5'-complete (phase 0 at the true start of the
# CDS). This holds for essentially all complete RefSeq/MANE CDS records;
# partial CDS records (5' truncated) would need phase-based trimming,
# which is not handled here.
# ----------------------------------------------------------------------
TRANSCRIPT_FEATURE_TYPES = {
    "mRNA", "transcript", "primary_transcript", "ncRNA", "lnc_RNA",
    "rRNA", "tRNA", "snRNA", "snoRNA", "miRNA", "pseudogenic_transcript",
    "V_gene_segment", "C_gene_segment", "guide_RNA", "RNase_P_RNA",
    "RNase_MRP_RNA", "telomerase_RNA", "Y_RNA", "antisense_RNA",
    "vault_RNA", "scRNA",
}


def parse_gff(gff_path, errors):
    """errors: list that transl_except rejections are appended to (as
    error-log records, see ERROR_LOG_COLUMNS)."""
    cds_tx = {}
    all_tx_ids = set()
    tx_line_geneid = {}
    tx_parent = {}           # transcript ID -> its Parent (gene) ID
    feature_geneid = {}      # any feature ID (e.g. gene lines) -> GeneID
    with open(gff_path) as f:
        for lineno, line in enumerate(f, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _, ftype, start, end, _, strand, _, attrs = fields
            attr_dict = dict(
                item.split("=", 1) for item in attrs.strip(";").split(";") if "=" in item
            )

            fid = attr_dict.get("ID")
            if fid and ftype not in TRANSCRIPT_FEATURE_TYPES and ftype != "CDS":
                gid = geneid_from_attrs(attr_dict)
                if gid:
                    feature_geneid.setdefault(fid, gid)

            if ftype in TRANSCRIPT_FEATURE_TYPES:
                tid = fid
                if tid:
                    all_tx_ids.add(tid)
                    gid = geneid_from_attrs(attr_dict)
                    if gid:
                        tx_line_geneid[tid] = gid
                    if attr_dict.get("Parent"):
                        tx_parent[tid] = attr_dict["Parent"].split(",")[0]

            if ftype != "CDS":
                continue
            parent = attr_dict.get("Parent")
            if not parent:
                continue
            transl_except = attr_dict.get("transl_except")
            cds_geneid = geneid_from_attrs(attr_dict)
            for p in parent.split(","):
                all_tx_ids.add(p)
                entry = cds_tx.setdefault(
                    p, {"chrom": chrom, "strand": strand, "cds": [], "transl_except_raw": {},
                        "geneid": None, "first_line": lineno}
                )
                entry["cds"].append((int(start), int(end)))
                if transl_except:
                    # raw tag -> line number of the first CDS line carrying it
                    entry["transl_except_raw"].setdefault(transl_except, lineno)
                if cds_geneid and entry["geneid"] is None:
                    entry["geneid"] = cds_geneid

    # NCBI GeneID per coding transcript, first found of: its CDS lines, its
    # transcript line, its parent gene line (some annotations, e.g. the
    # Liftoff CHM13 one, put db_xref only on the gene line). Also covers a
    # CDS parented directly by a gene. None if none of them has one.
    for tid, entry in cds_tx.items():
        if entry["geneid"] is None:
            entry["geneid"] = (tx_line_geneid.get(tid)
                               or feature_geneid.get(tx_parent.get(tid, ""))
                               or feature_geneid.get(tid))

    # Resolve every transl_except to a codon index now, once the transcript's
    # full CDS is known, so get_protein() just substitutes by index.
    n_tx_with_except = 0
    n_applied = 0
    rejected = []
    for tid, entry in cds_tx.items():
        resolved, bad = resolve_transl_except(entry)
        entry["transl_except"] = resolved
        entry["transl_except_rejected"] = [reason for reason, _ in bad]
        if resolved:
            n_tx_with_except += 1
            n_applied += len(resolved)
        for reason, tag_line in bad:
            rejected.append((tid, reason))
            errors.append(error_record(
                scope="gff", role="", source_file=gff_path, transcript_id=tid,
                location=cds_location(entry), gff_line=tag_line,
                error="transl_except_rejected", detail=reason,
            ))
    print(f"  {gff_path}: transl_except applied at {n_applied} codons in "
          f"{n_tx_with_except} transcripts; {len(rejected)} rejected")
    for tid, reason in rejected[:5]:
        print(f"    rejected transl_except: {tid}: {reason}")
    if len(rejected) > 5:
        print(f"    ... and {len(rejected) - 5} more rejected")
    n_sec_fallback = sum(1 for e in cds_tx.values() if e["sec_fallback"])
    print(f"    UGA->Sec fallback (broken aa:Sec tag) enabled for {n_sec_fallback} transcripts")

    return cds_tx, all_tx_ids


# ----------------------------------------------------------------------
# error.log: one tab-separated record per problem, with enough location
# info to find the entry. scope=gff records come from parsing a whole GFF
# (independent of the tmap); scope=row records come from a specific tmap
# row, and row_outcome says which outcome bucket that row ended up in, so
# you can tell whether the error is what caused the row to be skipped.
# ----------------------------------------------------------------------
ERROR_LOG_COLUMNS = [
    "scope", "role", "source_file", "transcript_id", "location", "gff_line",
    "tmap_row", "qry_id", "ref_id", "error", "detail", "row_outcome",
]


def error_record(**kw):
    rec = {c: "" for c in ERROR_LOG_COLUMNS}
    rec.update({k: ("" if v is None else v) for k, v in kw.items()})
    return rec


def cds_location(entry):
    """'chrom:start-end(strand)' spanning a transcript's CDS."""
    s = min(a for a, _ in entry["cds"])
    e = max(b for _, b in entry["cds"])
    return f"{entry['chrom']}:{s}-{e}({entry['strand']})"


# ----------------------------------------------------------------------
# transl_except handling (e.g. UGA recoded as selenocysteine, or a non-AUG
# start codon translated as Met). Values are URL-encoded in GFF3
# ('%2C' for the inner comma) and look like:
#   (pos:complement(30831787..30831789),aa:Sec)
#   (pos:join(100..101,200..200),aa:Sec)       <- codon split by an intron
# Several may be comma-separated on one line. aa:TERM (a partial stop
# completed by polyA) is ignored, since it's a stop either way.
# ----------------------------------------------------------------------
TRANSL_EXCEPT_RE = re.compile(r"\(pos:(.+?),aa:([A-Za-z]+)\)")
TRANSL_EXCEPT_POS_RE = re.compile(r"(\d+)(?:\.\.(\d+))?")


def resolve_transl_except(entry):
    """Map an entry's raw transl_except strings to {codon_index: aa1}.
    Returns (resolved_dict, rejected_list), where each rejected item is
    (reason_str, gff_line_of_the_tag). An exception is rejected
    (and NOT applied) unless its positions are exactly one in-frame codon
    lying within this transcript's CDS -- e.g. a position left in the source
    assembly's coordinates by a liftover tool that didn't update it.
    If a rejected tag is aa:Sec, sets entry["sec_fallback"] = True so that
    get_protein() reads every internal in-frame UGA as Sec instead."""
    raw = entry.pop("transl_except_raw")
    entry["sec_fallback"] = False
    resolved = {}
    rejected = []
    if not raw:
        return resolved, rejected

    strand = entry["strand"]
    # Segments in transcript (5'->3') order.
    segs = sorted(entry["cds"], key=lambda x: x[0], reverse=(strand == "-"))

    def cds_offset(pos):
        cum = 0
        for s, e in segs:
            if s <= pos <= e:
                return cum + (e - pos if strand == "-" else pos - s)
            cum += e - s + 1
        return None

    for value, tag_line in raw.items():
        for m in TRANSL_EXCEPT_RE.finditer(unquote(value)):
            pos_str, aa3 = m.group(1), m.group(2)
            if aa3.upper() == "TERM":
                continue
            if aa3.upper() == "OTHER":
                aa1 = "X"
            else:
                aa1 = protein_letters_3to1_extended.get(aa3.capitalize())
                if aa1 is None:
                    rejected.append((f"unknown amino acid '{aa3}' in {m.group(0)}", tag_line))
                    continue

            positions = set()
            for a, b in TRANSL_EXCEPT_POS_RE.findall(pos_str):
                positions.update(range(int(a), int(b or a) + 1))
            offsets = [cds_offset(p) for p in positions]
            if not offsets or None in offsets:
                rejected.append((f"position outside CDS in {m.group(0)}", tag_line))
                if aa1 == "U":
                    entry["sec_fallback"] = True
                continue
            offsets.sort()
            o = offsets[0]
            if offsets != [o, o + 1, o + 2] or o % 3 != 0:
                rejected.append((f"not a single in-frame codon in {m.group(0)}", tag_line))
                if aa1 == "U":
                    entry["sec_fallback"] = True
                continue
            resolved[o // 3] = aa1

    return resolved, rejected


def parse_multicopy_transcripts(gff_path):
    """Return transcript_id -> {"group", "copy", "geneid"} for every
    transcript whose gene was lifted over in more than one copy -- i.e. the
    primary copy AND every extra copy of any gene for which at least one
    gene line has extra_copy_number > 0 (Liftoff/LiftOn put that field only
    on gene lines: 0 = primary, >= 1 = extra).

    "group" identifies the set of copies of one gene: grouped by NCBI GeneID
    (every copy carries the same one), falling back to the gene ID with any
    trailing _N copy suffix removed if a gene line has no GeneID. "copy" is
    the transcript's own parent gene ID (one per copy). "geneid" is the NCBI
    GeneID, or None. Keys are the exact transcript IDs in the file (e.g.
    'rna-NM_x.1_1'); look up with multicopy_info(), not id_candidates(),
    which strips the copy suffix."""
    gene_group = {}          # gene ID -> group key
    gene_geneid = {}         # gene ID -> NCBI GeneID
    groups_with_extra = set()
    tx_parent = {}
    with open(gff_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            ftype, attrs = fields[2], fields[8]
            attr_dict = dict(
                item.split("=", 1) for item in attrs.strip(";").split(";") if "=" in item
            )
            fid = attr_dict.get("ID")
            if not fid:
                continue
            if ftype in TRANSCRIPT_FEATURE_TYPES:
                if attr_dict.get("Parent"):
                    tx_parent[fid] = attr_dict["Parent"].split(",")[0]
                continue
            if ftype != "gene":
                continue
            gid = geneid_from_attrs(attr_dict)
            key = f"GeneID:{gid}" if gid else "ID:" + re.sub(r"_\d+$", "", fid)
            gene_group[fid] = key
            gene_geneid[fid] = gid
            ecn = attr_dict.get("extra_copy_number", "")
            if ecn.isdigit() and int(ecn) > 0:
                groups_with_extra.add(key)
    return {tid: {"group": gene_group[g], "copy": g, "geneid": gene_geneid[g]}
            for tid, g in tx_parent.items()
            if gene_group.get(g) in groups_with_extra}


def multicopy_info(tx_id, multicopy_tx):
    """Lookup trying only the rna-/no-rna- prefix swap (never stripping the
    _N copy suffix). None if the transcript isn't in a multi-copy gene."""
    alt = tx_id[4:] if tx_id.startswith("rna-") else "rna-" + tx_id
    return multicopy_tx.get(tx_id) or multicopy_tx.get(alt)


def parse_transcript_mutations(gff_path):
    """Return transcript_id -> 'a;b;c' from the transcript line's LiftOn
    'mutation=' attribute (comma-separated in the GFF, e.g.
    'frameshift,stop_codon_gain'), URL-decoded. Transcripts without the
    attribute are absent. Keys are exact IDs; look up with exact_tx_get()."""
    out = {}
    with open(gff_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] not in TRANSCRIPT_FEATURE_TYPES:
                continue
            attr_dict = dict(
                item.split("=", 1) for item in fields[8].strip(";").split(";") if "=" in item
            )
            tid, mut = attr_dict.get("ID"), attr_dict.get("mutation")
            if tid and mut:
                out[tid] = ";".join(m.strip() for m in unquote(mut).split(",") if m.strip())
    return out


def exact_tx_get(tx_id, d, default=None):
    """d[tx_id], trying only the rna-/no-rna- prefix swap (never stripping
    the _N copy suffix, so each liftover copy keeps its own value)."""
    alt = tx_id[4:] if tx_id.startswith("rna-") else "rna-" + tx_id
    return d.get(tx_id, d.get(alt, default))


def parse_feature_geneids(gff_path):
    """Return feature_id -> NCBI GeneID for every gene-like and transcript
    feature in the GFF, coding or not. A transcript without its own
    Dbxref/db_xref GeneID inherits its parent gene's. Used to read the GeneID
    of whatever existing feature a tmap row's ref_id names."""
    own = {}
    parent = {}
    with open(gff_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] in ("exon", "CDS"):
                continue
            attr_dict = dict(
                item.split("=", 1) for item in fields[8].strip(";").split(";") if "=" in item
            )
            fid = attr_dict.get("ID")
            if not fid:
                continue
            gid = geneid_from_attrs(attr_dict)
            if gid:
                own[fid] = gid
            if attr_dict.get("Parent"):
                parent[fid] = attr_dict["Parent"].split(",")[0]
    out = dict(own)
    for fid, p in parent.items():
        if fid not in out and p in own:
            out[fid] = own[p]
    return out


def parse_pseudogene_ids(gff_path):
    """Return the set of every feature ID in the GFF whose own attributes
    contain 'pseudo=true' -- whether that's on a real pseudogene's gene
    line, its transcript line, or a purpose-built reference feature (e.g.
    one transcript whose exons cover every pseudogene). A tmap ref_id is
    checked against this set (via id_candidates) to flag pseudogene
    matches, since gffcompare's own reference-matching already handles
    pseudogene overlap detection as long as the GFF marks them this way."""
    ids = set()
    with open(gff_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            attrs = fields[8]
            if not PSEUDO_ATTR_RE.search(attrs):
                continue
            attr_dict = dict(
                item.split("=", 1) for item in attrs.strip(";").split(";") if "=" in item
            )
            tid = attr_dict.get("ID")
            if tid:
                ids.add(tid)
    return ids


GENEID_RE = re.compile(r"GeneID:(\d+)")
# Attribute names the cross-reference list appears under: GFF3's "Dbxref"
# (NCBI RefSeq, MANE) and "db_xref" (e.g. the Liftoff CHM13 annotation).
DBXREF_KEYS = ("Dbxref", "db_xref")


def geneid_from_attrs(attr_dict):
    """NCBI GeneID (string of digits) from a feature's Dbxref/db_xref, or None."""
    for key in DBXREF_KEYS:
        m = GENEID_RE.search(attr_dict.get(key, ""))
        if m:
            return m.group(1)
    return None


def parse_gene_order(gff_path):
    """Parse 'gene' features and return internal_gene_id ("gene-XXX") ->
    (left_geneid_or_END, right_geneid_or_END), based on coordinate order
    (start position) within each chromosome. Strand is deliberately
    ignored -- "left"/"right" here just means immediately-adjacent-by-
    position, since the same physical neighbor relationship can appear on
    either strand depending on the assembly.

    Genes lacking a resolvable Dbxref GeneID are skipped entirely: they
    never appear as a queryable gene, and never appear as a neighbor
    (their immediate neighbors become adjacent to each other instead).
    This keeps every value in the returned mapping a plain NCBI GeneID
    (or the "END" sentinel for a chromosome terminus), so two orderings
    built from different assemblies stay directly comparable."""
    raw = []  # (chrom, start, internal_id, geneid)
    with open(gff_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            chrom, start, attrs = fields[0], int(fields[3]), fields[8]
            attr_dict = dict(
                item.split("=", 1) for item in attrs.strip(";").split(";") if "=" in item
            )
            internal_id = attr_dict.get("ID")
            geneid = geneid_from_attrs(attr_dict)
            if not internal_id or not geneid:
                continue
            raw.append((chrom, start, internal_id, geneid))

    by_chrom = defaultdict(list)
    for chrom, start, internal_id, geneid in raw:
        by_chrom[chrom].append((start, internal_id, geneid))

    neighbors_by_internal_id = {}
    for chrom, entries in by_chrom.items():
        entries.sort(key=lambda x: x[0])
        n = len(entries)
        for idx, (start, internal_id, geneid) in enumerate(entries):
            left = entries[idx - 1][2] if idx > 0 else "END"
            right = entries[idx + 1][2] if idx < n - 1 else "END"
            neighbors_by_internal_id[internal_id] = (left, right)
    return neighbors_by_internal_id


def id_candidates(tx_id):
    """Generate ID variants to try, in order: as-given; rna-/no-rna-prefix
    swap; and the same two with a trailing liftoff copy suffix (e.g. '_1',
    '_2') stripped, since duplicated-copy IDs won't appear verbatim in the
    source annotation."""
    candidates = [tx_id]
    if tx_id.startswith("rna-"):
        candidates.append(tx_id[4:])
    else:
        candidates.append("rna-" + tx_id)

    stripped = re.sub(r"_\d+$", "", tx_id)
    if stripped != tx_id:
        candidates.append(stripped)
        if stripped.startswith("rna-"):
            candidates.append(stripped[4:])
        else:
            candidates.append("rna-" + stripped)

    return candidates


def lookup_tx(tx_id, tx_dict):
    for cand in id_candidates(tx_id):
        if cand in tx_dict:
            return tx_dict[cand]
    return None


def lookup_geneid(tx_id, cds_tx):
    """NCBI GeneID (string) of a coding transcript, or None if the
    transcript isn't found or carries no Dbxref GeneID."""
    entry = lookup_tx(tx_id, cds_tx)
    return entry["geneid"] if entry is not None else None


def lookup_neighbors(gene_id, neighbors_by_internal_id):
    for cand in id_candidates(gene_id):
        if cand in neighbors_by_internal_id:
            return neighbors_by_internal_id[cand]
    return None


def pseudo_id_candidates(ref_id):
    """id_candidates() plus a gene-/no-gene- prefix swap on each. Needed
    because gffread --gene2exon turns a transcript-less gene (typically a
    pseudogene) into its own transcript, so the tmap ref_id for those rows
    is the gene's ID, while pseudo=true sits on the gene line of the full
    GFF, possibly with a different gene- prefix."""
    out = []
    for c in id_candidates(ref_id):
        out.append(c)
        out.append(c[5:] if c.startswith("gene-") else "gene-" + c)
    return out


def cds_sequence(entry, fasta):
    """Spliced CDS DNA in transcript (5'->3') orientation. Caller must have
    checked entry['chrom'] is in the FASTA."""
    chrom_seq = fasta[entry["chrom"]]
    cds_sorted = sorted(entry["cds"], key=lambda x: x[0])
    # gff coords are 1-based inclusive; pyfaidx slicing is 0-based half-open
    full_seq = "".join(str(chrom_seq[s - 1:e]) for (s, e) in cds_sorted)
    if entry["strand"] == "-":
        full_seq = str(Seq(full_seq).reverse_complement())
    return full_seq


def get_protein(tx_id, cds_tx, all_tx_ids, fasta):
    """Returns (protein_str_or_None, status_str, location_str_or_None) where
    status is one of: 'ok', 'non_coding', 'not_found', 'chrom_not_found',
    'not_multiple_of_3'. location_str is 'chrom:start-end(strand)' spanning
    the transcript's CDS, given whenever a CDS entry was found at all (even
    if translation later fails for one of the reasons above) -- it's None
    only when the transcript itself couldn't be found in the GFF at all
    ('not_found'/'non_coding').
    Pure function -- no side effects -- so it's safe to call speculatively
    without corrupting any tallies; all counting happens once, in main(),
    after a row's outcome is fully decided."""
    entry = lookup_tx(tx_id, cds_tx)
    if entry is None:
        if any(cand in all_tx_ids for cand in id_candidates(tx_id)):
            return None, "non_coding", None
        return None, "not_found", None

    chrom = entry["chrom"]
    strand = entry["strand"]
    cds_sorted = sorted(entry["cds"], key=lambda x: x[0])
    location = f"{chrom}:{cds_sorted[0][0]}-{cds_sorted[-1][1]}({strand})"

    if chrom not in fasta:
        return None, "chrom_not_found", location

    full_seq = cds_sequence(entry, fasta)

    if len(full_seq) % 3 != 0:
        return None, "not_multiple_of_3", location

    # Translate every codon, apply any transl_except substitutions (so e.g. a
    # Sec-recoded UGA becomes 'U' instead of a stop, and a non-AUG start
    # becomes 'M'), then truncate at the first *remaining* stop codon. A
    # genuine premature stop still yields the truncated peptide as before.
    aa = translate_cds(entry, full_seq)
    protein = "".join(aa).split("*", 1)[0]
    return protein, "ok", location


def translate_cds(entry, full_seq):
    """Per-codon amino acids for a CDS (length a multiple of 3), with the
    UGA->Sec fallback and any valid transl_except applied. Stops are '*'."""
    aa = list(str(Seq(full_seq).translate()))
    for idx in sec_fallback_codons(entry, full_seq):
        aa[idx] = "U"
    for idx, aa1 in entry["transl_except"].items():
        if idx < len(aa):
            aa[idx] = aa1
    return aa


def cds_offset_to_genome(entry, offset):
    """Genomic (1-based) position of a 0-based offset into the spliced CDS,
    in transcript orientation."""
    strand = entry["strand"]
    for s, e in sorted(entry["cds"], key=lambda x: x[0], reverse=(strand == "-")):
        n = e - s + 1
        if offset < n:
            return e - offset if strand == "-" else s + offset
        offset -= n
    return None


def premature_stops(entry, full_seq):
    """[(codon_index, codon)] of every in-frame stop before the CDS's final
    codon, after Sec/transl_except substitutions -- i.e. the stops that
    truncate the translated protein."""
    aa = translate_cds(entry, full_seq)
    return [(i, full_seq[3 * i:3 * i + 3].upper())
            for i in range(len(aa) - 1) if aa[i] == "*"]


def sec_fallback_codons(entry, full_seq):
    """Codon indices read as Sec by the fallback for a broken aa:Sec
    transl_except tag: every in-frame UGA (TGA in the DNA; UAA/UAG are left
    as stops) except the CDS's final codon, which is the real terminal stop
    even when it's UGA. Empty unless entry["sec_fallback"] is set. Any valid
    tags on the same transcript are still applied on top by get_protein()."""
    if not entry["sec_fallback"]:
        return []
    n_codons = len(full_seq) // 3
    return [i for i in range(n_codons - 1) if full_seq[3 * i:3 * i + 3].upper() == "TGA"]


# ----------------------------------------------------------------------
# Percent identity: matches / max(len1, len2), from a global alignment.
# This guarantees a sequence < 50% the length of another can never score
# > 0.5 identity against it.
# ----------------------------------------------------------------------
_aligner = Align.PairwiseAligner()
_aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
_aligner.mode = "global"
_aligner.open_gap_score = -10
_aligner.extend_gap_score = -0.5

# BLOSUM62 has no Sec (U) or Pyl (O), and the aligner raises on unknown
# letters. Score them as their closest standard analogs (Cys, Lys) for the
# alignment only; identity is still counted on the original letters, so
# U only ever matches U.
_ALIGN_SUBST = str.maketrans({"U": "C", "O": "K"})


def percent_identity(seq1, seq2):
    if not seq1 or not seq2:
        return 0.0
    alignment = _aligner.align(seq1.translate(_ALIGN_SUBST), seq2.translate(_ALIGN_SUBST))[0]
    matches = 0
    for (s1, e1), (s2, e2) in zip(*alignment.aligned):
        a, b = seq1[s1:e1], seq2[s2:e2]
        matches += sum(x == y for x, y in zip(a, b))
    return matches / max(len(seq1), len(seq2))


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
class _OutcomeCounter(defaultdict):
    """defaultdict(int) that also remembers the last key written, so each
    row's error-log records can be tagged with the outcome bucket the row
    ended up in (every row increments exactly one bucket, exactly once)."""

    def __init__(self):
        super().__init__(int)
        self.last_key = None

    def __setitem__(self, key, value):
        super().__setitem__(key, value)
        self.last_key = key


ROW_ERROR_STATUSES = {"not_found", "chrom_not_found", "not_multiple_of_3"}


def row_error_records(i, qry_id, ref_id, role, tx_id, status, loc,
                      cds_dict, gff_path, fasta_path, fasta):
    """Error-log records for one transcript fetched for tmap row i: its
    get_protein() status if that's an error, plus any transl_except of its
    that was rejected (with what was done instead)."""
    recs = []
    entry = lookup_tx(tx_id, cds_dict)
    gff_line = entry["first_line"] if entry is not None else None
    base = dict(scope="row", role=role, source_file=gff_path, transcript_id=tx_id,
                location=loc, gff_line=gff_line, tmap_row=i, qry_id=qry_id, ref_id=ref_id)

    if status in ROW_ERROR_STATUSES:
        if status == "not_found":
            detail = f"not in GFF (tried IDs: {', '.join(id_candidates(tx_id))})"
        elif status == "chrom_not_found":
            detail = f"chromosome {entry['chrom']} not in {fasta_path}"
        else:
            cds_len = sum(e - s + 1 for s, e in entry["cds"])
            detail = (f"CDS length {cds_len} nt, {cds_len % 3} nt over a multiple of 3, "
                      f"{len(entry['cds'])} CDS segment(s)")
        recs.append(error_record(**base, error=status, detail=detail))

    if status == "ok" and entry is not None:
        full_seq = cds_sequence(entry, fasta)
        stops = premature_stops(entry, full_seq)
        if stops:
            n_codons = len(full_seq) // 3
            first_idx, first_codon = stops[0]
            p1 = cds_offset_to_genome(entry, 3 * first_idx)
            p3 = cds_offset_to_genome(entry, 3 * first_idx + 2)
            lo, hi = min(p1, p3), max(p1, p3)
            listed = ", ".join(f"{i + 1}:{c}" for i, c in stops[:10])
            more = f", +{len(stops) - 10} more" if len(stops) > 10 else ""
            recs.append(error_record(
                **base, error="premature_stop",
                detail=(f"{len(stops)} in-frame stop(s) before the last codon; first at codon "
                        f"{first_idx + 1} of {n_codons} ({first_codon}, "
                        f"{entry['chrom']}:{lo}-{hi}({entry['strand']})), protein truncated to "
                        f"{first_idx} aa; all (codon:triplet, 1-based): {listed}{more}"),
            ))

    if entry is not None and entry["transl_except_rejected"]:
        if not entry["sec_fallback"]:
            handling = "translated without it"
        elif status != "ok":
            handling = "UGA->Sec fallback not reached (no protein translated)"
        else:
            idxs = sec_fallback_codons(entry, cds_sequence(entry, fasta))
            shown = ", ".join(str(x + 1) for x in idxs) if idxs else "none"
            handling = (f"UGA->Sec fallback: {len(idxs)} internal in-frame UGA(s) "
                        f"read as Sec (codon {shown}, 1-based)")
        recs.append(error_record(
            **base, error="transl_except_rejected",
            detail="; ".join(entry["transl_except_rejected"]) + " -- " + handling,
        ))
    return recs


GENE_TRACE_COLUMNS = [
    "geneid", "gene_name", "copy_gene_id", "extra_copy_number", "transcript_id",
    "location", "class_code", "ref_id", "tmap_row", "outcome", "action",
]


def parse_genes_and_transcripts(gff_path):
    """Return (genes, transcripts):
      genes:       NCBI GeneID -> list of {"id", "name", "extra_copy_number",
                   "location"}, one per gene line (i.e. per liftover copy)
      transcripts: transcript ID -> {"parent", "geneid", "location"}; a
                   transcript without its own GeneID inherits its parent's."""
    genes = defaultdict(list)
    gene_geneid = {}
    transcripts = {}
    with open(gff_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _, ftype, start, end, _, strand, _, attrs = fields
            if ftype != "gene" and ftype not in TRANSCRIPT_FEATURE_TYPES:
                continue
            attr_dict = dict(
                item.split("=", 1) for item in attrs.strip(";").split(";") if "=" in item
            )
            fid = attr_dict.get("ID")
            if not fid:
                continue
            loc = f"{chrom}:{start}-{end}({strand})"
            gid = geneid_from_attrs(attr_dict)
            if ftype == "gene":
                gene_geneid[fid] = gid
                if gid:
                    genes[gid].append({
                        "id": fid,
                        "name": attr_dict.get("Name") or attr_dict.get("gene") or "",
                        "extra_copy_number": attr_dict.get("extra_copy_number", ""),
                        "location": loc,
                    })
            else:
                transcripts[fid] = {"parent": (attr_dict.get("Parent") or "").split(",")[0],
                                    "geneid": gid, "location": loc}
    for t in transcripts.values():
        if t["geneid"] is None:
            t["geneid"] = gene_geneid.get(t["parent"])
    return genes, transcripts


def _norm_tx(tx_id):
    return tx_id[4:] if tx_id.startswith("rna-") else tx_id


def write_gene_trace(tmap_all, tmap, outcome_by_row, action_by_row):
    """For each GeneID in GENE_LIST_PATH, write one GENE_TRACE row per lifted
    transcript (every copy) saying where it ended up, plus a row for genes
    that were never lifted over at all. Prints a short per-gene summary."""
    try:
        with open(GENE_LIST_PATH) as f:
            wanted = [l.strip() for l in f if l.strip() and not l.startswith("#")]
    except FileNotFoundError:
        print(f"\nGene trace skipped: {GENE_LIST_PATH} not found.")
        return

    mane_genes, _ = parse_genes_and_transcripts(MANE_GFF)
    lifted_genes, lifted_tx = parse_genes_and_transcripts(LIFTED_GFF)

    # tmap rows by transcript ID (rna- prefix normalized). Filtered-tmap row
    # indices follow tmap_all order with the '=' rows removed.
    tmap_rows = defaultdict(list)
    filtered_idx = 0
    for qry_id, ref_id, cc in zip(tmap_all["qry_id"].astype(str),
                                   tmap_all["ref_id"].astype(str),
                                   tmap_all["class_code"].astype(str)):
        idx = None
        if cc != "=":
            idx = filtered_idx
            filtered_idx += 1
        tmap_rows[_norm_tx(qry_id)].append((cc, ref_id, idx))
    assert filtered_idx == len(tmap)

    tx_by_copy = defaultdict(list)
    for tid, t in lifted_tx.items():
        tx_by_copy[t["parent"]].append(tid)

    records = []
    print(f"\nGene trace for {len(wanted)} GeneIDs from {GENE_LIST_PATH}:")
    for gid in wanted:
        name = next((g["name"] for g in mane_genes.get(gid, []) + lifted_genes.get(gid, [])
                     if g["name"]), "")
        base = {"geneid": gid, "gene_name": name}
        copies = lifted_genes.get(gid, [])
        gene_records = []
        if not copies:
            outcome = ("not_mapped (in MANE, no lifted gene)" if gid in mane_genes
                       else "not_in_MANE (GeneID not in MANE GFF)")
            gene_records.append({**base, "outcome": outcome})
        for copy in copies:
            cbase = {**base, "copy_gene_id": copy["id"],
                     "extra_copy_number": copy["extra_copy_number"]}
            tids = sorted(tx_by_copy.get(copy["id"], []))
            if not tids:
                gene_records.append({**cbase, "location": copy["location"],
                                     "outcome": "no_transcripts (gene line only)"})
            for tid in tids:
                tbase = {**cbase, "transcript_id": tid, "location": lifted_tx[tid]["location"]}
                rows = tmap_rows.get(_norm_tx(tid), [])
                if not rows:
                    gene_records.append({**tbase, "outcome": "not_in_tmap"})
                for cc, ref_id, idx in rows:
                    if idx is None:
                        outcome, action = "exact_match (class_code '=', not compared)", ""
                    else:
                        outcome = outcome_by_row.get(idx, "")
                        action = action_by_row.get(idx, "")
                    gene_records.append({**tbase, "class_code": cc, "ref_id": ref_id,
                                         "tmap_row": "" if idx is None else idx,
                                         "outcome": outcome, "action": action})
        records.extend(gene_records)

        parts = []
        for r in gene_records:
            label = r["outcome"] + (f" -> {r['action']}" if r.get("action") else "")
            parts.append(f"{r['transcript_id']}: {label}" if r.get("transcript_id") else label)
        print(f"  {gid} ({name or '?'}), {len(copies)} lifted copies: " + "; ".join(parts))

    with open(GENE_TRACE, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=GENE_TRACE_COLUMNS, delimiter="\t",
                           lineterminator="\n", restval="")
        w.writeheader()
        w.writerows(records)
    print(f"Gene trace written to: {GENE_TRACE} ({len(records)} records)")


# ----------------------------------------------------------------------
# Same-ID check for the three-way branch: compare against the existing
# transcript with our own accession on the gene gffcompare paired us with,
# when there is one, rather than gffcompare's (possibly different) pick.
# ----------------------------------------------------------------------
REF_MATCH_LABELS = [
    ("same_id_on_gene", "Same transcript ID found on the paired gene (compared against it)"),
    ("id_on_other_gene", "Transcript ID exists in existing annotation, but not on the paired gene"),
    ("id_not_in_annotation", "No transcript with this ID in existing annotation"),
    ("paired_gene_unknown", "Paired gene could not be determined"),
]

def tx_accession(tx_id):
    """Transcript ID as compared against the existing annotation: only an
    'rna-' prefix is dropped. Copy suffixes are kept, so 'rna-NM_1.2_1'
    only matches an existing 'NM_1.2_1', not 'NM_1.2'."""
    return tx_id[4:] if tx_id.startswith("rna-") else tx_id


def parse_transcript_index(gff_path):
    """Index the transcripts of a GFF for the same-ID check:
      "parent":    transcript ID -> parent gene ID
      "location":  transcript ID -> 'chrom:start-end(strand)'
      "by_acc":    accession (tx_accession) -> [transcript IDs]
      "gene_ids":  every gene-like feature ID (transcript parents included)"""
    parent, location, by_acc, gene_ids = {}, {}, defaultdict(list), set()
    with open(gff_path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom, _, ftype, start, end, _, strand, _, attrs = fields
            if ftype in ("exon", "CDS"):
                continue
            attr_dict = dict(
                item.split("=", 1) for item in attrs.strip(";").split(";") if "=" in item
            )
            fid = attr_dict.get("ID")
            if not fid:
                continue
            if ftype in TRANSCRIPT_FEATURE_TYPES:
                p = (attr_dict.get("Parent") or "").split(",")[0]
                if p:
                    parent[fid] = p
                    gene_ids.add(p)
                location[fid] = f"{chrom}:{start}-{end}({strand})"
                by_acc[tx_accession(fid)].append(fid)
            else:
                gene_ids.add(fid)
    return {"parent": parent, "location": location, "by_acc": by_acc, "gene_ids": gene_ids}


def choose_existing_transcript(qry_id, ref_id, ref_gene_id, index):
    """Return (transcript_to_compare, status, detail). status is one of the
    REF_MATCH_LABELS keys. The paired gene is the parent of gffcompare's
    ref_id in the existing GFF, falling back to the tmap's ref_gene_id if
    that names a gene in the file. Only 'same_id_on_gene' changes the
    transcript compared; every other status keeps gffcompare's ref_id."""
    parent = index["parent"]
    gene = next((parent[c] for c in id_candidates(ref_id) if c in parent), None)
    if gene is None:
        gene = next((c for c in pseudo_id_candidates(ref_gene_id) if c in index["gene_ids"]), None)

    acc = tx_accession(qry_id)
    same_acc = index["by_acc"].get(acc, [])

    if gene is None:
        return ref_id, "paired_gene_unknown", (
            f"gffcompare pick {ref_id} (gene {ref_gene_id}) not placed on a gene in existing annotation")

    on_gene = [t for t in same_acc if parent.get(t) == gene]
    if on_gene:
        # Prefer gffcompare's own pick if it is one of them (only possible
        # with both 'X' and 'rna-X' present), then the first by ID.
        ref_cands = set(id_candidates(ref_id))
        chosen = (next((t for t in on_gene if t in ref_cands), None)
                  or sorted(on_gene)[0])
        if chosen in ref_cands:
            detail = f"already gffcompare's pick, on gene {gene}"
        else:
            detail = f"replaced gffcompare's pick {ref_id}, on gene {gene}"
        return chosen, "same_id_on_gene", detail

    if same_acc:
        shown = "; ".join(f"{t} on gene {parent.get(t, '?')} {index['location'].get(t, '')}"
                          for t in sorted(same_acc)[:3])
        more = f"; +{len(same_acc) - 3} more" if len(same_acc) > 3 else ""
        return ref_id, "id_on_other_gene", f"paired gene {gene}; {acc} is at: {shown}{more}"

    return ref_id, "id_not_in_annotation", f"paired gene {gene}; no {acc} in existing annotation"


def accession_in_annotation(qry_id, index):
    """For novel-locus rows: ('id_in_annotation', where) if a transcript with
    our accession exists anywhere in the existing annotation, else
    ('id_not_in_annotation', '')."""
    acc = tx_accession(qry_id)
    hits = index["by_acc"].get(acc, [])
    if not hits:
        return "id_not_in_annotation", f"no {acc} in existing annotation"
    shown = "; ".join(f"{t} on gene {index['parent'].get(t, '?')} {index['location'].get(t, '')}"
                      for t in sorted(hits)[:3])
    more = f"; +{len(hits) - 3} more" if len(hits) > 3 else ""
    return "id_in_annotation", f"novel locus; {acc} is at: {shown}{more}"


def find_copies_matching_own_gene(tmap_all, lifted_multicopy, existing_geneids):
    """From the FULL tmap (including class_code '=' rows), return the set of
    (group, copy) for multi-copy lifted transcripts whose matched existing
    feature (ref_id) has the same NCBI GeneID as the lifted gene -- i.e.
    copies that found their real gene. Copies without a GeneID never match."""
    found = set()
    for qry_id, ref_id in zip(tmap_all["qry_id"].astype(str), tmap_all["ref_id"].astype(str)):
        if ref_id == "-":
            continue
        info = multicopy_info(qry_id, lifted_multicopy)
        if info is None or info["geneid"] is None:
            continue
        ref_geneid = next((existing_geneids[c] for c in pseudo_id_candidates(ref_id)
                           if c in existing_geneids), None)
        if ref_geneid == info["geneid"]:
            found.add((info["group"], info["copy"]))
    return found


def main():
    tmap_all = pd.read_csv(TMAP_PATH, sep="\t")
    tmap = tmap_all[tmap_all["class_code"] != "="].reset_index(drop=True)

    print(f"Loaded tmap: {len(tmap)} rows with class_code != '='")

    # Everything that ends up in ERROR_LOG, in order: whole-GFF transl_except
    # rejections first (from parsing), then per-row records.
    errors = []

    print("Parsing GFFs...")
    mane_cds, mane_all_ids = parse_gff(MANE_GFF, errors)
    lifted_cds, lifted_all_ids = parse_gff(LIFTED_GFF, errors)
    existing_cds, existing_all_ids = parse_gff(EXISTING_GFF, errors)
    existing_pseudo_ids = parse_pseudogene_ids(EXISTING_GFF)
    existing_tx_index = parse_transcript_index(EXISTING_GFF)
    lifted_multicopy = parse_multicopy_transcripts(LIFTED_GFF)
    lifted_mutations = parse_transcript_mutations(LIFTED_GFF)
    copies_matching_own_gene = find_copies_matching_own_gene(
        tmap_all, lifted_multicopy, parse_feature_geneids(EXISTING_GFF))

    print("Parsing gene order for neighbor-swap checks...")
    mane_gene_neighbors = parse_gene_order(MANE_GFF)
    lifted_gene_neighbors = parse_gene_order(LIFTED_GFF)

    print("Indexing genome FASTAs (pyfaidx)...")
    mane_fasta = Fasta(MANE_FASTA)
    chm13_fasta = Fasta(CHM13_FASTA)

    # Every row is assigned to exactly one outcome bucket, exactly once, so
    # row_outcomes.values() always sums to len(tmap) -- there is no other
    # counting anywhere else in the loop, and get_protein() has no side
    # effects, so nothing can be tallied twice or silently dropped.
    row_outcomes = _OutcomeCounter()
    # Purely informational (doesn't affect row_outcomes / the total): for
    # rows dropped because MANE itself has no CDS, what the mapped and
    # existing transcripts' CDS status looked like anyway -- since those
    # rows short-circuit before ever reaching the mapped-vs-existing
    # crosstab below. Only applies to the normal (three-way) branch, since
    # that's the only branch that ever looks at the existing transcript.
    mane_non_coding_diagnostic = defaultdict(int)
    # Also informational: how many three-way rows had the neighbor check
    # skipped outright because the gene couldn't be resolved (missing
    # GeneID) in one or both orderings -- these are NOT discarded, just
    # left unchecked, since missing data isn't evidence of a real swap.
    neighbor_check_unavailable = 0
    # Informational: Add/Review split of the lifted_better rows (a sub-split
    # of compared_lifted_better, so NOT part of row_outcomes / the total),
    # plus how often each Review reason fired (a row can have both).
    row_actions = defaultdict(int)
    review_reason_counts = defaultdict(int)
    results = []

    # Per-row error records are held here until the row's outcome is known
    # (the next iteration's start, or after the loop), then tagged with it.
    pending_row_errors = []
    # For the gene trace: tmap row index -> outcome bucket, and -> Add/Review.
    outcome_by_row = {}
    action_by_row = {}
    # Informational: same-ID check result for every three-way row that got
    # past the neighbor check (independent of the row's final outcome).
    ref_match_counts = defaultdict(int)
    # Informational: for novel_compared rows, whether the transcript ID is
    # already somewhere in the existing annotation.
    novel_id_counts = defaultdict(int)
    # Add/Review for the two-way rows, keyed (prefix, action); and per-reason
    # counts for novel-locus Reviews (a row can have both reasons).
    twoway_actions = defaultdict(int)
    novel_review_reason_counts = defaultdict(int)
    current_row = [None]

    def flush_row_errors():
        for rec in pending_row_errors:
            rec["row_outcome"] = row_outcomes.last_key or ""
            errors.append(rec)
        pending_row_errors.clear()
        if current_row[0] is not None:
            outcome_by_row[current_row[0]] = row_outcomes.last_key or ""

    for i, row in tqdm(tmap.iterrows(), total=len(tmap), desc="Comparing proteins"):
        flush_row_errors()
        row_outcomes.last_key = None
        current_row[0] = i
        qry_id = str(row["qry_id"])
        qry_gene_id = str(row["qry_gene_id"])
        ref_id = str(row["ref_id"])

        is_pseudogene_ref = ref_id != "-" and any(
            cand in existing_pseudo_ids for cand in pseudo_id_candidates(ref_id)
        )

        if ref_id == "-" or is_pseudogene_ref:
            # Either gffcompare found no reference match at all, or the
            # match it found is flagged pseudo=true. Either way there's no
            # meaningful existing-annotation protein to compare against --
            # just compare protein1 (MANE) vs protein2 (lifted).

            # If the gene was lifted over in multiple copies and a DIFFERENT
            # copy of it (primary or extra) is matched in the tmap to an
            # existing gene with the same NCBI GeneID -- i.e. another copy
            # found the real gene -- then this copy landing on a pseudogene
            # is the spurious one: discard it outright, before any protein
            # is fetched (so regardless of identity or CDS status).
            if is_pseudogene_ref:
                info = multicopy_info(qry_id, lifted_multicopy)
                if info is not None and any(
                    g == info["group"] and c != info["copy"]
                    for g, c in copies_matching_own_gene
                ):
                    row_outcomes["pseudogene_discarded_extra_copy"] += 1
                    continue

            protein1, status1, loc1 = get_protein(qry_id, mane_cds, mane_all_ids, mane_fasta)
            protein2, status2, loc2 = get_protein(qry_id, lifted_cds, lifted_all_ids, chm13_fasta)
            pending_row_errors.extend(
                row_error_records(i, qry_id, ref_id, "mane", qry_id, status1, loc1,
                                  mane_cds, MANE_GFF, MANE_FASTA, mane_fasta)
                + row_error_records(i, qry_id, ref_id, "lifted", qry_id, status2, loc2,
                                    lifted_cds, LIFTED_GFF, CHM13_FASTA, chm13_fasta)
            )

            prefix = "pseudogene" if is_pseudogene_ref else "novel"

            if status1 == "non_coding":
                row_outcomes[f"{prefix}_mane_non_coding"] += 1
                continue
            if status1 == "not_found":
                tqdm.write(f"WARNING: qry_id '{qry_id}' not found in MANE GFF ({MANE_GFF}); skipping row {i}")
                row_outcomes[f"{prefix}_mane_not_found"] += 1
                continue
            if status2 == "non_coding":
                row_outcomes[f"{prefix}_lifted_non_coding"] += 1
                continue
            if status2 == "not_found":
                tqdm.write(f"WARNING: qry_id '{qry_id}' not found in lifted-over GFF ({LIFTED_GFF}); skipping row {i}")
                row_outcomes[f"{prefix}_lifted_not_found"] += 1
                continue
            if protein1 is None:
                row_outcomes[f"{prefix}_mane_{status1}"] += 1  # not_multiple_of_3 / chrom_not_found
                continue
            if protein2 is None:
                row_outcomes[f"{prefix}_lifted_{status2}"] += 1
                continue

            id_12 = percent_identity(protein1, protein2)

            # Remaining pseudogene matches are only kept if the lifted
            # protein is identical to MANE's.
            if is_pseudogene_ref and id_12 < 1.0:
                row_outcomes["pseudogene_discarded_not_identical"] += 1
                continue

            # Novel loci: does our transcript ID already exist anywhere in
            # the existing annotation? (Reported only; no comparison.)
            novel_status = novel_detail = None
            if not is_pseudogene_ref:
                novel_status, novel_detail = accession_in_annotation(qry_id, existing_tx_index)

            # Two-way Add/Review: pseudogene matches always go to Review;
            # novel loci are Add if MANE-vs-lifted identity >= the threshold
            # and no LiftOn mutation mentions 'codon' or 'frameshift'.
            mutations = exact_tx_get(qry_id, lifted_mutations, "")
            if is_pseudogene_ref:
                tw_reasons = ["pseudogene_match"]
            else:
                tw_reasons = []
                if id_12 < NOVEL_ADD_MIN_IDENTITY:
                    tw_reasons.append(f"identity_below_{NOVEL_ADD_MIN_IDENTITY:.2f}")
                bad = [m for m in mutations.split(";")
                       if m and any(term in m.lower() for term in NOVEL_BLOCKING_MUTATION_TERMS)]
                if bad:
                    tw_reasons.append("mutation:" + ",".join(bad))
            tw_action = "Review" if tw_reasons else "Add"
            tw_review_reason = ";".join(tw_reasons) if tw_reasons else None
            action_by_row[i] = tw_action + (f" ({tw_review_reason})" if tw_review_reason else "")
            twoway_actions[(prefix, tw_action)] += 1
            if not is_pseudogene_ref:
                for r in tw_reasons:
                    novel_review_reason_counts[r.split(":", 1)[0]] += 1

            results.append({
                "ref_gene_id": row["ref_gene_id"],
                "ref_id": ref_id,
                "class_code": row["class_code"],
                "qry_gene_id": qry_gene_id,
                "qry_id": qry_id,
                "identity_mane_vs_lifted": id_12,
                "identity_mane_vs_existing": None,
                "verdict": "pseudogene_match" if is_pseudogene_ref else "novel_locus",
                "len_mane": len(protein1),
                "len_lifted": len(protein2),
                "len_existing": None,
                "location_mane": loc1,
                "location_lifted": loc2,
                "location_existing": None,
                "geneid_lifted": None,
                "geneid_existing": None,
                "action": tw_action,
                "review_reason": tw_review_reason,
                "mutations_lifted": mutations,
                "ref_id_compared": None,
                "ref_match_status": novel_status,
                "ref_match_detail": novel_detail,
            })
            if novel_status is not None:
                novel_id_counts[novel_status] += 1
            row_outcomes[f"{prefix}_compared"] += 1
            continue

        # --- Normal case: a real, non-pseudogene reference transcript ---

        # Gene-order neighbor check: if this gene's two immediate neighbors
        # (by GeneID, position-based, strand-independent) in the original
        # MANE ordering share NEITHER member with its two immediate
        # neighbors in the lifted-over ordering, it's likely been swapped
        # with a closely-related paralog during mapping -- discard it.
        orig_neighbors = lookup_neighbors(qry_gene_id, mane_gene_neighbors)
        new_neighbors = lookup_neighbors(qry_gene_id, lifted_gene_neighbors)
        if orig_neighbors is None or new_neighbors is None:
            neighbor_check_unavailable += 1
        elif set(orig_neighbors).isdisjoint(new_neighbors):
            row_outcomes["discarded_neighbor_mismatch"] += 1
            continue

        # Same-ID check: if the gene gffcompare paired us with carries a
        # transcript with our own accession, compare against that one
        # instead of gffcompare's pick. Otherwise keep gffcompare's pick.
        ref_used, ref_match_status, ref_match_detail = choose_existing_transcript(
            qry_id, ref_id, str(row["ref_gene_id"]), existing_tx_index)
        ref_match_counts[ref_match_status] += 1

        protein1, status1, loc1 = get_protein(qry_id, mane_cds, mane_all_ids, mane_fasta)
        protein2, status2, loc2 = get_protein(qry_id, lifted_cds, lifted_all_ids, chm13_fasta)
        protein3, status3, loc3 = get_protein(ref_used, existing_cds, existing_all_ids, chm13_fasta)
        pending_row_errors.extend(
            row_error_records(i, qry_id, ref_id, "mane", qry_id, status1, loc1,
                              mane_cds, MANE_GFF, MANE_FASTA, mane_fasta)
            + row_error_records(i, qry_id, ref_id, "lifted", qry_id, status2, loc2,
                                lifted_cds, LIFTED_GFF, CHM13_FASTA, chm13_fasta)
            + row_error_records(i, qry_id, ref_id, "existing", ref_used, status3, loc3,
                                existing_cds, EXISTING_GFF, CHM13_FASTA, chm13_fasta)
        )

        if status1 == "non_coding":
            row_outcomes["mane_non_coding"] += 1
            # Informational only -- record what the mapped/existing CDS
            # status was anyway, even though it didn't determine the outcome.
            mane_non_coding_diagnostic[f"lifted={status2}, existing={status3}"] += 1
            continue
        if status1 == "not_found":
            tqdm.write(f"WARNING: qry_id '{qry_id}' not found in MANE GFF ({MANE_GFF}); skipping row {i}")
            row_outcomes["mane_not_found"] += 1
            continue

        # Mapped-vs-existing CDS presence -- checked before either "not
        # found" case so it takes priority whenever both apply.
        if status2 == "non_coding" and status3 == "non_coding":
            row_outcomes["both_lack_cds"] += 1
            continue
        if status2 == "non_coding":
            row_outcomes["existing_has_cds_mapped_does_not"] += 1
            continue
        if status3 == "non_coding":
            print(f"Mapped has CDS existing does not: {qry_id}")
            # Written to the CSV as a novel locus for Review (the paired
            # existing transcript has no CDS, so there's no protein to
            # compare against). Still counted in its own bucket below.
            id_12_nc = (percent_identity(protein1, protein2)
                        if protein1 is not None and protein2 is not None else None)
            nc_reason = "existing_transcript_non_coding"
            action_by_row[i] = f"Review ({nc_reason})"
            results.append({
                "ref_gene_id": row["ref_gene_id"],
                "ref_id": ref_id,
                "class_code": row["class_code"],
                "qry_gene_id": qry_gene_id,
                "qry_id": qry_id,
                "identity_mane_vs_lifted": id_12_nc,
                "identity_mane_vs_existing": None,
                "verdict": "novel_locus",
                "len_mane": len(protein1) if protein1 is not None else None,
                "len_lifted": len(protein2) if protein2 is not None else None,
                "len_existing": None,
                "location_mane": loc1,
                "location_lifted": loc2,
                "location_existing": None,
                "geneid_lifted": lookup_geneid(qry_id, lifted_cds),
                "geneid_existing": None,
                "action": "Review",
                "review_reason": nc_reason,
                "mutations_lifted": exact_tx_get(qry_id, lifted_mutations, ""),
                "ref_id_compared": ref_used,
                "ref_match_status": ref_match_status,
                "ref_match_detail": ref_match_detail,
            })
            row_outcomes["mapped_has_cds_existing_does_not"] += 1
            continue

        if status2 == "not_found":
            tqdm.write(f"WARNING: qry_id '{qry_id}' not found in lifted-over GFF ({LIFTED_GFF}); skipping row {i}")
            row_outcomes["lifted_not_found"] += 1
            continue
        if status3 == "not_found":
            tqdm.write(f"WARNING: ref_id '{ref_used}' not found in existing annotation ({EXISTING_GFF}); skipping row {i}")
            row_outcomes["existing_not_found"] += 1
            continue

        if protein1 is None:
            row_outcomes[f"mane_{status1}"] += 1  # not_multiple_of_3 / chrom_not_found
            continue
        if protein2 is None:
            row_outcomes[f"lifted_{status2}"] += 1
            continue
        if protein3 is None:
            row_outcomes[f"existing_{status3}"] += 1
            continue

        id_12 = percent_identity(protein1, protein2)  # MANE vs lifted-over-on-chm13
        id_13 = percent_identity(protein1, protein3)  # MANE vs existing-on-chm13

        if id_12 > id_13:
            verdict = "lifted_better"
        elif id_13 > id_12:
            verdict = "existing_better"
        else:
            verdict = "equal"

        # Add/Review triage, lifted_better rows only.
        geneid_lifted = lookup_geneid(qry_id, lifted_cds)
        geneid_existing = lookup_geneid(ref_used, existing_cds)
        action = None
        review_reason = None
        if verdict == "lifted_better":
            reasons = []
            if geneid_lifted is None or geneid_existing is None:
                reasons.append("geneid_missing")
            elif geneid_lifted != geneid_existing:
                reasons.append("geneid_mismatch")
            if id_13 < ADD_MIN_EXISTING_IDENTITY:
                reasons.append(f"existing_identity_below_{ADD_MIN_EXISTING_IDENTITY:.2f}")
            action = "Review" if reasons else "Add"
            # An Add that was compared against an existing transcript with
            # our own ID on the paired gene replaces that transcript.
            if action == "Add" and ref_match_status == "same_id_on_gene":
                action = "Replace"
            review_reason = ";".join(reasons) if reasons else None
            action_by_row[i] = action + (f" ({review_reason})" if review_reason else "")
            row_actions[action] += 1
            for r in reasons:
                review_reason_counts[r] += 1

        results.append({
            "ref_gene_id": row["ref_gene_id"],
            "ref_id": ref_id,
            "class_code": row["class_code"],
            "qry_gene_id": qry_gene_id,
            "qry_id": qry_id,
            "identity_mane_vs_lifted": id_12,
            "identity_mane_vs_existing": id_13,
            "verdict": verdict,
            "len_mane": len(protein1),
            "len_lifted": len(protein2),
            "len_existing": len(protein3),
            "location_mane": loc1,
            "location_lifted": loc2,
            "location_existing": loc3,
            "geneid_lifted": geneid_lifted,
            "geneid_existing": geneid_existing,
            "action": action,
            "review_reason": review_reason,
            "mutations_lifted": exact_tx_get(qry_id, lifted_mutations, ""),
            "ref_id_compared": ref_used,
            "ref_match_status": ref_match_status,
            "ref_match_detail": ref_match_detail,
        })
        row_outcomes[f"compared_{verdict}"] += 1

    flush_row_errors()

    results_df = pd.DataFrame(results)
    results_df.to_csv(OUT_CSV, index=False)

    with open(ERROR_LOG, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=ERROR_LOG_COLUMNS, delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(errors)

    n_existing_better = row_outcomes["compared_existing_better"]
    n_lifted_better = row_outcomes["compared_lifted_better"]
    n_equal = row_outcomes["compared_equal"]
    n_novel = row_outcomes["novel_compared"]
    n_pseudogene = row_outcomes["pseudogene_compared"]
    n_three_way = n_existing_better + n_lifted_better + n_equal

    print("\n" + "=" * 60)
    print(f"Discarded outright: gene-order neighbor mismatch (three-way branch only -- "
          f"neither original MANE neighbor is still adjacent after mapping): "
          f"{row_outcomes['discarded_neighbor_mismatch']}")
    if neighbor_check_unavailable:
        print(f"  (Neighbor check skipped -- gene not resolvable via GeneID in one or "
              f"both orderings -- for {neighbor_check_unavailable} further rows; not discarded.)")

    print(f"\nThree-way comparisons (MANE vs. lifted vs. existing): {n_three_way}")
    print(f"  Existing annotation closer to MANE:    {n_existing_better}")
    print(f"  Lifted-over annotation closer to MANE: {n_lifted_better}")
    print(f"  Equal identity:                        {n_equal}")
    print(f"  Lifted-better split -- Add (same GeneID as paired existing transcript, "
          f"existing identity >= {ADD_MIN_EXISTING_IDENTITY:.2f}): {row_actions['Add']}")
    print(f"  Lifted-better split -- Replace (Add criteria met, and compared against the "
          f"existing transcript with the same ID on the paired gene): {row_actions['Replace']}")
    print(f"  Lifted-better split -- Review: {row_actions['Review']}")
    for reason, count in sorted(review_reason_counts.items()):
        print(f"    Review reason {reason}: {count}")
    print("  Existing transcript used, three-way rows past the neighbor check "
          "(gffcompare's pick is used unless the first applies):")
    for status, label in REF_MATCH_LABELS:
        print(f"    {label}: {ref_match_counts[status]}")

    print(f"\nTwo-way comparisons for genuinely novel loci (ref_id == '-'; "
          f"MANE vs. lifted only): {n_novel}")
    print(f"  Transcript ID already present in existing annotation (elsewhere): "
          f"{novel_id_counts['id_in_annotation']}")
    print(f"  Transcript ID not in existing annotation: "
          f"{novel_id_counts['id_not_in_annotation']}")
    print(f"  Novel-locus split -- Add (MANE identity >= {NOVEL_ADD_MIN_IDENTITY:.2f}, no mutation "
          f"containing {'/'.join(NOVEL_BLOCKING_MUTATION_TERMS)}): {twoway_actions[('novel', 'Add')]}")
    print(f"  Novel-locus split -- Review: {twoway_actions[('novel', 'Review')]}")
    for reason, count in sorted(novel_review_reason_counts.items()):
        print(f"    Review reason {reason}: {count}")

    print(f"\nTwo-way comparisons where the matched reference is a pseudogene "
          f"(MANE vs. lifted only, kept separate from the three-way comparisons "
          f"above): {n_pseudogene}")
    print(f"  Pseudogene-match split -- Review (all pseudogene matches): "
          f"{twoway_actions[('pseudogene', 'Review')]}")
    print(f"  Pseudogene matches discarded -- gene lifted in multiple copies "
          f"and another copy matched an existing gene with the same GeneID: "
          f"{row_outcomes['pseudogene_discarded_extra_copy']}")
    print(f"  Pseudogene matches discarded -- lifted protein not identical to MANE "
          f"(identity < 1.0): {row_outcomes['pseudogene_discarded_not_identical']}")

    print("\nRows skipped because mapped and/or existing transcript lacks a CDS "
          "(three-way branch only):")
    print(f"  Mapped has CDS, existing does not:  {row_outcomes['mapped_has_cds_existing_does_not']}")
    print(f"  Existing has CDS, mapped does not:  {row_outcomes['existing_has_cds_mapped_does_not']}")
    print(f"  Neither has a CDS:                  {row_outcomes['both_lack_cds']}")
    print("  (These three counts only cover rows where MANE itself has a CDS -- see below for what")
    print("   the mapped/existing CDS status looked like in the MANE-non-coding rows dropped first.)")

    n_mane_non_coding_total = (
        row_outcomes["mane_non_coding"]
        + row_outcomes["novel_mane_non_coding"]
        + row_outcomes["pseudogene_mane_non_coding"]
    )
    print(f"\nRows skipped, MANE transcript is non-coding (no CDS), across all branches: "
          f"{n_mane_non_coding_total}")
    if mane_non_coding_diagnostic:
        print("  Informational breakdown of mapped/existing CDS status for these rows")
        print("  (three-way branch only; did NOT determine the outcome -- MANE non-coding alone did):")
        for label, count in sorted(mane_non_coding_diagnostic.items(), key=lambda kv: -kv[1]):
            print(f"    {label}: {count}")
    print(f"Rows skipped, ref_id not found in existing annotation: {row_outcomes['existing_not_found']}")

    bad_cds_reasons = {k: v for k, v in row_outcomes.items()
                        if "not_multiple_of_3" in k or "chrom_not_found" in k}
    if bad_cds_reasons:
        print("Rows skipped, CDS malformed (not a multiple of 3, or its chromosome wasn't found), by source:")
        for label, count in sorted(bad_cds_reasons.items()):
            print(f"  {label}: {count}")

    print(f"\nFull per-row results written to: {OUT_CSV}")
    print(f"Errors written to: {ERROR_LOG} ({len(errors)} records)")

    write_gene_trace(tmap_all, tmap, outcome_by_row, action_by_row)

    # Sanity check: every row was assigned to exactly one outcome, so this
    # must equal the tmap row count printed at the start of the run.
    accounted = sum(row_outcomes.values())
    print(f"\nSanity check: {accounted} rows accounted for across all categories above "
          f"(should equal the {len(tmap)} loaded at the start).")
    if accounted != len(tmap):
        print("  MISMATCH -- please report this, it indicates a bug in the row accounting.")


if __name__ == "__main__":
    sys.exit(main())