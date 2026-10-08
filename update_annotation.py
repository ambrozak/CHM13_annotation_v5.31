#!/usr/bin/env python3
"""
Update the existing CHM13 annotation with MANE-on-CHM13 (LiftOn) models,
using the per-transcript actions from compare_proteins.py plus manual
review decisions.

Inputs
  RESULTS_CSV   protein_comparison_results.csv (compare_proteins.py output)
  MANUAL_CSV    MANUAL_REVIEW_DECISIONS.csv: qry_id (lifted gene ID),
                decision, matchedwithpseudo, novellocus, notes
  LIFTED_GFF    MANE lifted onto CHM13 by LiftOn
  EXISTING_GFF  the existing CHM13 annotation (Liftoff RefSeq)

Decisions
  Add            add the lifted transcript to the existing gene it was
                 compared against (must have the same NCBI GeneID). On a
                 novel locus / pseudogene match there is no such gene, so
                 Add there means AddGene.
  Replace        remove the existing transcript with the same ID (wherever
                 it is) and insert the lifted one in its place, under the
                 same gene, keeping the ID. That gene must have the lifted
                 transcript's GeneID.
  Review         use the manual decision for the lifted gene; with no
                 manual decision the row is left alone (logged).
  AddGene        new gene feature with the lifted transcript(s) as children.
  AddPseudo      new gene feature, pseudo=true, no children.
  ReplacePseudo  remove the existing pseudogene the lifted gene overlaps,
                 then AddGene.
  Leave          nothing.
  A manual decision for a gene takes precedence over the automatic action
  of its transcripts; overriding an automatic Add/Replace is logged.

Format (hard-coded from the CHM13 Liftoff annotation, see FORMAT below)
  gene        ID;gene_name;db_xref;description;gbkey;gene;gene_biotype;
              gene_synonym;[pseudo;]sequence_ID;extra_copy_number;copy_num_ID
  transcript  ID;Parent;db_xref;gbkey;gene;product;tag;transcript_biotype;
              matches_ref_protein;extra_copy_number
  exon        Parent;db_xref;gene;product;tag;transcript_biotype;
              exon_number;extra_copy_number
  CDS         Parent;db_xref;gbkey;gene;note;product;protein_id;tag;
              exon_number;extra_copy_number[;transl_except]
  tag is always NEW_TAG ("MANE Select") on added transcripts/exons/CDS;
  CDS note= comes from LiftOn's CDS Note=.
  Children: all exons, then all CDS, each in ascending coordinate order.
  Values come from LiftOn; source column = NEW_SOURCE. Liftoff QC values:
  sequence_ID <- LiftOn dna_identity, matches_ref_protein <- (LiftOn
  protein_identity == 1); coverage and valid_ORF omitted. Empty values
  are left out. transl_except is appended to CDS lines that have one
  (dropping it would change the encoded protein).

IDs
  New transcripts keep their accession without 'rna-'. Any new ID (gene or
  transcript) that already exists in the file is FLAGged and suffixed
  _1 (then _2, ...) until unique. A Replace reuses the removed ID.

Outputs
  OUT_GFF        updated annotation (every untouched line byte-identical)
  CHANGES_TSV    one row per planned action: what was done or why not
"""

import csv
import re
import sys
from collections import defaultdict, Counter

RESULTS_CSV = "protein_comparison_results.csv"
MANUAL_CSV = "manual_review_decisions.csv"
LIFTED_GFF = "MANE_on_chm13.gff"
EXISTING_GFF = "./chm13v5.3.gff"

OUT_GFF = "chm13v5.31_MANE_updated.gff"
CHANGES_TSV = "annotation_changes.tsv"

NEW_SOURCE = "LiftOn"
# Feature type used for mRNA transcripts in the CHM13 annotation.
TRANSCRIPT_TYPE = "transcript"
# tag= written on every added transcript, exon and CDS line.
NEW_TAG = "MANE Select"

KNOWN_DECISIONS = {"Add", "Replace", "Review", "Leave", "AddGene", "AddPseudo", "ReplacePseudo"}

# Parts of a transcript. A transcript is any feature with an ID whose Parent
# is a gene and whose type is not one of these -- no list of transcript
# type names needed ('transcript' in CHM13, 'mRNA' in LiftOn, ...).
SUBFEATURE_TYPES = {"exon", "CDS", "five_prime_UTR", "three_prime_UTR",
                    "start_codon", "stop_codon"}
GENEID_RE = re.compile(r"GeneID:(\d+)")


def parse_attrs(s):
    """Ordered [(key, value)] of a GFF3 column 9 (values left encoded)."""
    out = []
    for item in s.strip().strip(";").split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out.append((k.strip(), v))
    return out


def dbxref_of(a):
    return a.get("Dbxref") or a.get("db_xref") or ""


def geneid_of(a):
    m = GENEID_RE.search(dbxref_of(a))
    return m.group(1) if m else None


def overlaps(a_start, a_end, b_start, b_end):
    return a_start <= b_end and b_start <= a_end


# ----------------------------------------------------------------------
# Existing annotation: raw lines (written back untouched) and indexes.
# ----------------------------------------------------------------------
class Existing:
    def __init__(self, path):
        self.lines = []
        self.feat = {}            # line idx -> {type, chrom, start, end, strand, id, parent}
        self.attrs = {}           # line idx -> dict, gene & transcript lines only
        self.id_idx = {}
        self.children = defaultdict(list)   # parent id -> [line idx]
        self.genes_by_chrom = defaultdict(list)

        id_groups = defaultdict(set)

        with open(path) as f:
            for idx, line in enumerate(f):
                self.lines.append(line)
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue
                chrom, _, ftype, start, end, _, strand, _, attr_s = fields
                a = dict(parse_attrs(attr_s))
                fid, parent = a.get("ID"), (a.get("Parent") or "").split(",")[0] or None
                self.feat[idx] = {"type": ftype, "chrom": chrom, "start": int(start),
                                  "end": int(end), "strand": strand, "id": fid, "parent": parent}
                if fid:
                    # One ID on several lines of the same type and parent is
                    # one multi-line feature (e.g. CDS segments), not a dup.
                    id_groups[fid].add((ftype, parent))
                    self.id_idx.setdefault(fid, idx)
                if parent:
                    self.children[parent].append(idx)

                if ftype == "gene" and not parent:
                    self.attrs[idx] = a
                    self.genes_by_chrom[chrom].append(idx)
                elif fid and parent and ftype not in SUBFEATURE_TYPES:
                    self.attrs[idx] = a   # candidate transcript; checked below

        self.dup_ids = {k for k, g in id_groups.items() if len(g) > 1}
        self._index_blocks()
        self._index_gene_geneids()

    def _index_blocks(self):
        """Root gene of every feature, and the last line of each root's
        block (the line after which children may be appended)."""
        root_cache = {}

        def root_of(fid):
            seen = []
            cur = fid
            while cur is not None and cur not in root_cache:
                seen.append(cur)
                idx = self.id_idx.get(cur)
                p = self.feat[idx]["parent"] if idx is not None else None
                if p is None:
                    root_cache[cur] = cur
                    break
                cur = p
            r = root_cache.get(cur, cur)
            for s in seen:
                root_cache[s] = r
            return r

        self.root = {}
        self.block_end = {}
        for idx, f in self.feat.items():
            r = root_of(f["id"]) if f["id"] else root_of(f["parent"]) if f["parent"] else None
            if r is None:
                continue
            self.root[idx] = r
            if idx > self.block_end.get(r, -1):
                self.block_end[r] = idx
        for c in self.genes_by_chrom:
            self.genes_by_chrom[c].sort(key=lambda i: (self.feat[i]["start"], i))

    def _index_gene_geneids(self):
        """NCBI GeneID of each gene: own db_xref, else its transcripts'."""
        self.gene_geneid = {}
        for idx, a in self.attrs.items():
            f = self.feat[idx]
            if f["type"] == "gene":
                g = geneid_of(a)
                if g:
                    self.gene_geneid[f["id"]] = g
        gene_ids = {self.feat[i]["id"] for i, a in self.attrs.items() if self.feat[i]["type"] == "gene"}
        for idx, a in self.attrs.items():
            f = self.feat[idx]
            if f["parent"] in gene_ids and f["parent"] not in self.gene_geneid:
                g = geneid_of(a)
                if g:
                    self.gene_geneid[f["parent"]] = g

    def descendants(self, fid):
        out, stack = [], [fid]
        while stack:
            for i in self.children.get(stack.pop(), []):
                out.append(i)
                if self.feat[i]["id"]:
                    stack.append(self.feat[i]["id"])
        return out

    def resolve(self, some_id):
        """Line idx of an ID as written in the results CSV (tolerates rna-/
        gene- prefix differences)."""
        if not some_id:
            return None
        cands = [some_id]
        for p in ("rna-", "gene-"):
            cands.append(some_id[len(p):] if some_id.startswith(p) else p + some_id)
        for c in cands:
            if c in self.id_idx:
                return self.id_idx[c]
        return None


# ----------------------------------------------------------------------
# Lifted (LiftOn) annotation
# ----------------------------------------------------------------------
class Lifted:
    def __init__(self, path):
        self.feats = {}                     # id -> feature dict
        self.children = defaultdict(list)   # parent id -> [feature]
        with open(path) as f:
            for line in f:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue
                chrom, _, ftype, start, end, _, strand, phase, attr_s = fields
                a = dict(parse_attrs(attr_s))
                feat = {"type": ftype, "chrom": chrom, "start": int(start), "end": int(end),
                        "strand": strand, "phase": phase, "a": a}
                if a.get("ID"):
                    self.feats.setdefault(a["ID"], feat)
                for p in (a.get("Parent") or "").split(","):
                    if p:
                        self.children[p].append(feat)

    def get(self, some_id):
        for c in (some_id, some_id[4:] if some_id.startswith("rna-") else "rna-" + some_id):
            if c in self.feats:
                return c, self.feats[c]
        return None, None

    def transcripts_of(self, gene_id):
        return [(c["a"]["ID"], c) for c in self.children.get(gene_id, [])
                if c["type"] not in SUBFEATURE_TYPES and c["a"].get("ID")]


def tx_accession(tx_id):
    return tx_id[4:] if tx_id.startswith("rna-") else tx_id


def gene_symbol(gene):
    a = gene["a"]
    return a.get("Name") or a.get("gene") or re.sub(r"^gene-", "", a["ID"])


# ----------------------------------------------------------------------
# Formatting new features in the CHM13 format
# ----------------------------------------------------------------------
def pick_xref(dbxref, prefixes):
    """First entry of a comma-separated Dbxref with one of the prefixes,
    tried in order (e.g. 'HGNC' -> 'HGNC:HGNC:3025')."""
    entries = dbxref.split(",") if dbxref else []
    for pref in prefixes:
        for e in entries:
            if e.split(":")[0] == pref:
                return e
    return None


def attr_line(f, ftype, start, end, phase, pairs):
    """One GFF line; (key, value) pairs with an empty value are left out."""
    attrs = ";".join(f"{k}={v}" for k, v in pairs if v not in (None, ""))
    return "\t".join([f["chrom"], NEW_SOURCE, ftype, str(start), str(end), ".",
                       f["strand"], phase, attrs]) + "\n"


def format_gene(gene, new_id, txs, pseudo=False):
    a = gene["a"]
    symbol = gene_symbol(gene)
    ecn = a.get("extra_copy_number", "0")
    seq_id = next((t["a"]["dna_identity"] for _, t in txs if t["a"].get("dna_identity")), None)
    return attr_line(gene, "gene", gene["start"], gene["end"], ".", [
        ("ID", new_id),
        ("gene_name", symbol),
        # CHM13 gene lines carry the HGNC xref, or GeneID when there is none
        ("db_xref", pick_xref(dbxref_of(a), ("HGNC", "GeneID"))),
        ("description", a.get("description")),
        ("gbkey", "Gene"),
        ("gene", a.get("gene") or symbol),
        ("gene_biotype", "pseudogene" if pseudo else a.get("gene_biotype")),
        ("gene_synonym", a.get("gene_synonym")),
        ("pseudo", "true" if pseudo else None),
        ("sequence_ID", seq_id),
        ("extra_copy_number", ecn),
        ("copy_num_ID", f"{symbol}_{ecn}"),
    ])


def format_transcript(lifted, tx_id, tx, new_id, parent_id, symbol, ecn):
    """Transcript line followed by its exon lines then its CDS lines."""
    a = tx["a"]
    geneid_x = pick_xref(dbxref_of(a), ("GeneID",))
    biotype = a.get("transcript_biotype") or "mRNA"
    tag = NEW_TAG
    product = a.get("product")
    pid = a.get("protein_identity")
    try:
        matches = None if pid is None else ("True" if float(pid) == 1.0 else "False")
    except ValueError:
        matches = None
    lines = [attr_line(tx, TRANSCRIPT_TYPE, tx["start"], tx["end"], ".", [
        ("ID", new_id), ("Parent", parent_id), ("db_xref", geneid_x),
        ("gbkey", a.get("gbkey") or "mRNA"), ("gene", symbol), ("product", product),
        ("tag", tag), ("transcript_biotype", biotype),
        ("matches_ref_protein", matches), ("extra_copy_number", ecn),
    ])]

    kids = lifted.children.get(tx_id, [])
    exons = sorted([k for k in kids if k["type"] == "exon"], key=lambda k: k["start"])
    cds = sorted([k for k in kids if k["type"] == "CDS"], key=lambda k: k["start"])
    # exon_number counts in transcript (5'->3') order
    n = len(exons)
    exon_no = {id(e): (k + 1 if tx["strand"] != "-" else n - k) for k, e in enumerate(exons)}

    def exon_number_of(seg):
        for e in exons:
            if e["start"] <= seg["start"] and seg["end"] <= e["end"]:
                return exon_no[id(e)]
        return None

    for e in exons:
        lines.append(attr_line(e, "exon", e["start"], e["end"], ".", [
            ("Parent", new_id), ("db_xref", geneid_x), ("gene", symbol),
            ("product", e["a"].get("product") or product), ("tag", tag),
            ("transcript_biotype", biotype), ("exon_number", exon_no[id(e)]),
            ("extra_copy_number", ecn),
        ]))
    for c in cds:
        lines.append(attr_line(c, "CDS", c["start"], c["end"], c["phase"], [
            ("Parent", new_id), ("db_xref", geneid_x), ("gbkey", "CDS"), ("gene", symbol),
            # CHM13 writes RefSeq's CDS Note (e.g. "isoform a is encoded by
            # transcript variant 2") as lowercase note=
            ("note", c["a"].get("Note") or c["a"].get("note")),
            ("product", c["a"].get("product")), ("protein_id", c["a"].get("protein_id")),
            ("tag", tag), ("exon_number", exon_number_of(c)),
            ("extra_copy_number", ecn), ("transl_except", c["a"].get("transl_except")),
        ]))
    return lines


# ----------------------------------------------------------------------
# Planning
# ----------------------------------------------------------------------
def load_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def build_plans(results, manual, lifted, log):
    """One plan per (decision, lifted gene, [transcripts]). Gene-level
    decisions (AddGene, AddPseudo, ReplacePseudo) are merged per gene."""
    per_tx = []          # (decision, source, gene_id, tx_id, row, note)
    genes_with_rows = set()

    for row in results:
        qry = row.get("qry_id", "")
        tx_id, tx = lifted.get(qry)
        if tx is None:
            log(dict(lifted_transcript=qry, decision=row.get("action", ""), status="FLAG",
                     flags=f"qry_id not found in {LIFTED_GFF}"))
            continue
        gene_id = (tx["a"].get("Parent") or "").split(",")[0]
        genes_with_rows.add(gene_id)
        auto = row.get("action", "") or ""
        verdict = row.get("verdict", "")
        man = manual.get(gene_id)
        note = ""
        if man:
            decision, source = man["decision"], "manual"
            if auto in ("Add", "Replace") and decision != auto:
                note = f"manual {decision} overrides automatic {auto}"
        elif auto == "Review":
            log(dict(lifted_gene=gene_id, lifted_transcript=tx_id, decision="Review",
                     source="auto", status="not applied",
                     notes=f"Review with no manual decision ({row.get('review_reason', '')})"))
            continue
        elif auto in ("Add", "Replace"):
            decision, source = auto, "auto"
        else:
            continue
        if decision == "Add" and verdict in ("novel_locus", "pseudogene_match"):
            note = (note + "; " if note else "") + f"Add on {verdict} -> AddGene"
            decision = "AddGene"
        per_tx.append((decision, source, gene_id, tx_id, row, note))

    for gene_id, m in manual.items():
        if gene_id in genes_with_rows:
            continue
        gene = lifted.feats.get(gene_id)
        if gene is None:
            log(dict(lifted_gene=gene_id, decision=m["decision"], source="manual", status="FLAG",
                     flags=f"manual qry_id not a gene ID in {LIFTED_GFF}"))
            continue
        txs = lifted.transcripts_of(gene_id)
        if m["decision"] == "AddPseudo" or not txs:
            per_tx.append((m["decision"], "manual", gene_id, None, None,
                           "gene not in results CSV"))
        for tx_id, _ in txs:
            per_tx.append((m["decision"], "manual", gene_id, tx_id, None,
                           "gene not in results CSV"))

    plans, gene_level = [], {}
    for decision, source, gene_id, tx_id, row, note in per_tx:
        if decision in ("AddGene", "AddPseudo", "ReplacePseudo"):
            key = (decision, gene_id)
            p = gene_level.get(key)
            if p is None:
                p = gene_level[key] = dict(decision=decision, source=source, gene_id=gene_id,
                                           txs=[], rows=[], notes=[])
                plans.append(p)
            if tx_id and tx_id not in p["txs"]:
                p["txs"].append(tx_id)
            if row is not None:
                p["rows"].append(row)
            if note:
                p["notes"].append(note)
        else:
            plans.append(dict(decision=decision, source=source, gene_id=gene_id,
                              txs=[tx_id] if tx_id else [], rows=[row] if row else [],
                              notes=[note] if note else []))
    return plans


# ----------------------------------------------------------------------
# Applying
# ----------------------------------------------------------------------
class Editor:
    def __init__(self, ex, lifted, log):
        self.ex, self.lifted, self.log = ex, lifted, log
        self.removed = set()                 # line idx
        self.inserts = defaultdict(list)     # line idx -> [(priority, sort key, seq, lines)]
        self.overrides = {}                  # line idx -> replacement line
        self.used_ids = set(ex.id_idx)
        self.seq = 0
        self.flags_total = 0

    # -- helpers -------------------------------------------------------
    def alloc_id(self, wanted, kind, flags):
        new = wanted
        n = 0
        while new in self.used_ids:
            n += 1
            new = f"{wanted}_{n}"
        if new != wanted:
            flags.append(f"{kind} ID {wanted} already exists; inserted as {new}")
        self.used_ids.add(new)
        return new

    def insert_before(self, idx, lines, priority, sort_key=0):
        self.seq += 1
        self.inserts[idx].append((priority, sort_key, self.seq, lines))

    def remove_block(self, fid):
        idx = self.ex.id_idx[fid]
        gone = [idx] + self.ex.descendants(fid)
        self.removed.update(gone)
        ids = list(dict.fromkeys(self.ex.feat[i]["id"] for i in gone if self.ex.feat[i]["id"]))
        self.used_ids.difference_update(ids)
        return ids

    def extend_gene(self, gene_id, start, end, notes):
        gidx = self.ex.id_idx[gene_id]
        line = self.overrides.get(gidx, self.ex.lines[gidx])
        fields = line.rstrip("\n").split("\t")
        s, e = int(fields[3]), int(fields[4])
        if start < s or end > e:
            fields[3], fields[4] = str(min(s, start)), str(max(e, end))
            self.overrides[gidx] = "\t".join(fields) + "\n"
            notes.append(f"gene {gene_id} extended {s}-{e} -> {fields[3]}-{fields[4]}")

    def new_gene_position(self, chrom, start):
        """Line idx to insert a new top-level gene before, keeping the file
        position-sorted: before the first remaining gene on the chromosome
        starting after it, else after the chromosome's last gene block."""
        genes = [i for i in self.ex.genes_by_chrom.get(chrom, []) if i not in self.removed]
        for i in genes:
            if self.ex.feat[i]["start"] > start:
                return i
        if genes:
            return self.ex.block_end[self.ex.feat[genes[-1]]["id"]] + 1
        return len(self.ex.lines)

    def target_gene(self, tx_id, tx, gene_id, row, flags):
        """Existing gene to add to / replace in: the parent of the existing
        transcript the row was compared against; its GeneID must equal the
        lifted gene's. Without a results row: the unique existing gene with
        that GeneID overlapping the lifted transcript."""
        lifted_gid = geneid_of(tx["a"]) or geneid_of(self.lifted.feats.get(gene_id, {"a": {}})["a"])
        if lifted_gid is None:
            flags.append("lifted transcript has no GeneID")
            return None
        ref = (row.get("ref_id_compared") or row.get("ref_id")) if row else None
        if ref and ref != "-":
            ridx = self.ex.resolve(ref)
            if ridx is None:
                flags.append(f"compared existing transcript {ref} not found in {EXISTING_GFF}")
                return None
            g = self.ex.feat[ridx]["parent"] or self.ex.feat[ridx]["id"]
            ggid = self.ex.gene_geneid.get(g)
            if ggid != lifted_gid:
                flags.append(f"GeneID mismatch: lifted {lifted_gid} vs existing gene {g} "
                             f"({ggid}) -- likely a bug in compare_proteins.py; skipped")
                return None
            return g
        cands = [self.ex.feat[i]["id"] for i in self.ex.genes_by_chrom.get(tx["chrom"], [])
                 if i not in self.removed
                 and self.ex.gene_geneid.get(self.ex.feat[i]["id"]) == lifted_gid
                 and overlaps(self.ex.feat[i]["start"], self.ex.feat[i]["end"], tx["start"], tx["end"])]
        if len(cands) == 1:
            return cands[0]
        flags.append(f"{len(cands)} existing genes with GeneID {lifted_gid} overlap the lifted "
                     f"transcript (need exactly 1)" + (f": {', '.join(cands)}" if cands else ""))
        return None

    # -- actions -------------------------------------------------------
    def apply(self, p):
        flags, notes = [], list(p["notes"])
        rec = dict(lifted_gene=p["gene_id"], lifted_transcript=",".join(p["txs"]),
                   decision=p["decision"], source=p["source"])
        d = p["decision"]
        try:
            if d in ("Leave", "Review"):
                rec["status"] = "not applied"
                notes.append("Leave" if d == "Leave" else "manual decision is Review")
            elif d in ("Add", "Replace"):
                self._add_or_replace(p, rec, flags, notes)
            elif d in ("AddGene", "ReplacePseudo"):
                self._add_gene(p, rec, flags, notes, replace_pseudo=(d == "ReplacePseudo"))
            elif d == "AddPseudo":
                self._add_pseudo(p, rec, flags, notes)
            else:
                flags.append(f"unhandled decision {d}")
        finally:
            if flags:
                rec["status"] = "FLAG" if not rec.get("new_ids") else "applied (flagged)"
                self.flags_total += 1
            rec.setdefault("status", "applied")
            rec["flags"] = "; ".join(flags)
            rec["notes"] = "; ".join(dict.fromkeys(n for n in notes if n))
            self.log(rec)

    def _add_or_replace(self, p, rec, flags, notes):
        if not p["txs"]:
            flags.append("no lifted transcript")
            return
        tx_id = p["txs"][0]
        _, tx = self.lifted.get(tx_id)
        row = p["rows"][0] if p["rows"] else None
        if p["decision"] == "Replace":
            self._replace(p, tx_id, tx, rec, flags, notes)
            return

        gene = self.target_gene(tx_id, tx, p["gene_id"], row, flags)
        if gene is None:
            return
        rec["target_gene"] = gene
        gidx = self.ex.id_idx[gene]
        ecn = self.ex.attrs[gidx].get("extra_copy_number", "0")
        symbol = self.ex.attrs[gidx].get("gene") or gene
        new_id = self.alloc_id(tx_accession(tx_id), "transcript", flags)
        lines = format_transcript(self.lifted, tx_id, tx, new_id, gene, symbol, ecn)
        self.insert_before(self.ex.block_end[self.ex.root[gidx]] + 1, lines, priority=0)
        rec["new_ids"] = new_id
        self.extend_gene(gene, tx["start"], tx["end"], notes)

    def _replace(self, p, tx_id, tx, rec, flags, notes):
        """Replace the existing transcript whose ID is this transcript's ID
        (without 'rna-'), wherever it is in the file -- independent of which
        gene gffcompare paired the row with, so it works the same for
        automatic and manual Replace, novel-locus and pseudogene rows, and
        genes with no results row. The new transcript keeps that ID, goes
        under the same gene, in the same place in the file. The gene's
        GeneID must equal the lifted transcript's."""
        acc = tx_accession(tx_id)
        oidx = self.ex.id_idx.get(acc)
        if oidx is None:
            oidx = self.ex.id_idx.get("rna-" + acc)
        if oidx is None:
            flags.append(f"no existing transcript with ID {acc} to replace")
            return
        old = self.ex.feat[oidx]["id"]
        if oidx in self.removed:
            flags.append(f"existing transcript {old} was already removed by another action")
            return
        gene = self.ex.feat[oidx]["parent"]
        gidx = self.ex.id_idx.get(gene) if gene else None
        if gidx is None or self.ex.feat[gidx]["type"] != "gene":
            flags.append(f"existing {old} is not a transcript of a gene (Parent={gene})")
            return
        if old in self.ex.dup_ids:
            flags.append(f"ID {old} is used by more than one feature in {EXISTING_GFF}; "
                         f"replacing the first one (line {oidx + 1})")
        lifted_gid = geneid_of(tx["a"]) or geneid_of(self.lifted.feats.get(p["gene_id"], {"a": {}})["a"])
        gene_gid = self.ex.gene_geneid.get(gene)
        if lifted_gid != gene_gid:
            flags.append(f"GeneID mismatch: lifted {lifted_gid} vs gene {gene} ({gene_gid}) "
                         f"of existing {old}; skipped")
            return
        o = self.ex.feat[oidx]
        if (o["chrom"], o["strand"]) != (tx["chrom"], tx["strand"]) or \
                not overlaps(o["start"], o["end"], tx["start"], tx["end"]):
            flags.append(f"existing {old} is at {o['chrom']}:{o['start']}-{o['end']}({o['strand']}) "
                         f"but the lifted transcript is at {tx['chrom']}:{tx['start']}-{tx['end']}"
                         f"({tx['strand']}) -- not overlapping on the same strand; skipped")
            return
        rec["target_gene"] = gene
        ecn = self.ex.attrs[gidx].get("extra_copy_number", "0")
        symbol = self.ex.attrs[gidx].get("gene") or gene
        removed_ids = self.remove_block(old)
        new_id = self.alloc_id(old, "transcript", flags)
        lines = format_transcript(self.lifted, tx_id, tx, new_id, gene, symbol, ecn)
        self.insert_before(oidx, lines, priority=0)
        rec["removed_ids"] = ",".join(removed_ids)
        rec["new_ids"] = new_id
        self.extend_gene(gene, tx["start"], tx["end"], notes)

    def _find_pseudogene(self, p, gene, flags):
        for row in p["rows"]:
            ridx = self.ex.resolve(row.get("ref_id", ""))
            if ridx is None:
                continue
            f = self.ex.feat[ridx]
            gid = f["id"] if f["type"] == "gene" else f["parent"]
            gidx = self.ex.id_idx.get(gid)
            if gidx is not None and self.ex.attrs.get(gidx, {}).get("pseudo") == "true":
                return gid
        cands = [self.ex.feat[i]["id"] for i in self.ex.genes_by_chrom.get(gene["chrom"], [])
                 if i not in self.removed and self.ex.attrs[i].get("pseudo") == "true"
                 and overlaps(self.ex.feat[i]["start"], self.ex.feat[i]["end"], gene["start"], gene["end"])]
        if len(cands) == 1:
            return cands[0]
        flags.append(f"{len(cands)} existing pseudogenes overlap {p['gene_id']} and the results "
                     f"row did not identify one (need exactly 1)"
                     + (f": {', '.join(cands)}" if cands else ""))
        return None

    def _add_gene(self, p, rec, flags, notes, replace_pseudo):
        gene = self.lifted.feats.get(p["gene_id"])
        if gene is None:
            flags.append(f"gene {p['gene_id']} not in {LIFTED_GFF}")
            return
        if replace_pseudo:
            pg = self._find_pseudogene(p, gene, flags)
            if pg is None:
                return
            rec["removed_ids"] = ",".join(self.remove_block(pg))
        txs = [(t, self.lifted.feats[t]) for t in p["txs"]] or self.lifted.transcripts_of(p["gene_id"])
        if not txs:
            flags.append("lifted gene has no transcripts")
            return
        symbol = gene_symbol(gene)
        new_gene = self.alloc_id(symbol, "gene", flags)
        ecn = gene["a"].get("extra_copy_number", "0")
        lines = [format_gene(gene, new_gene, txs)]
        new_ids = [new_gene]
        for tx_id, tx in txs:
            tid = self.alloc_id(tx_accession(tx_id), "transcript", flags)
            new_ids.append(tid)
            lines += format_transcript(self.lifted, tx_id, tx, tid, new_gene, symbol, ecn)
        self.insert_before(self.new_gene_position(gene["chrom"], gene["start"]), lines,
                           priority=1, sort_key=gene["start"])
        rec["new_ids"] = ",".join(new_ids)
        if gene["chrom"] not in self.ex.genes_by_chrom:
            flags.append(f"chromosome {gene['chrom']} not in {EXISTING_GFF}; appended at end")

    def _add_pseudo(self, p, rec, flags, notes):
        gene = self.lifted.feats.get(p["gene_id"])
        if gene is None:
            flags.append(f"gene {p['gene_id']} not in {LIFTED_GFF}")
            return
        txs = self.lifted.transcripts_of(p["gene_id"])
        new_gene = self.alloc_id(gene_symbol(gene), "gene", flags)
        line = format_gene(gene, new_gene, txs, pseudo=True)
        self.insert_before(self.new_gene_position(gene["chrom"], gene["start"]), [line],
                           priority=1, sort_key=gene["start"])
        rec["new_ids"] = new_gene

    # -- output --------------------------------------------------------
    def write(self, path):
        with open(path, "w") as out:
            n = len(self.ex.lines)
            for i in range(n + 1):
                for _, _, _, lines in sorted(self.inserts.get(i, []), key=lambda x: x[:3]):
                    out.writelines(lines)
                if i < n and i not in self.removed:
                    out.write(self.overrides.get(i, self.ex.lines[i]))


def validate(path, ex):
    """Re-read the output: duplicate IDs, dangling Parents, children outside
    their parent's span, CDS lengths not a multiple of 3."""
    ids, feats = defaultdict(set), []
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            fl = line.rstrip("\n").split("\t")
            if len(fl) < 9:
                continue
            a = dict(parse_attrs(fl[8]))
            feats.append((fl[1], fl[2], int(fl[3]), int(fl[4]), a))
            if a.get("ID"):
                ids[a["ID"]].add((fl[2], a.get("Parent")))
    spans = {}
    for src, t, s, e, a in feats:
        if a.get("ID"):
            spans[a["ID"]] = (s, e)
    dup = [k for k, g in ids.items() if len(g) > 1 and k not in ex.dup_ids]
    dangling, outside = [], []
    cds_len = defaultdict(int)
    for src, t, s, e, a in feats:
        for p in (a.get("Parent") or "").split(","):
            if not p:
                continue
            if p not in spans:
                dangling.append(f"{t}@{s}-{e} -> {p}")
            elif src == NEW_SOURCE:
                ps, pe = spans[p]
                if s < ps or e > pe:
                    outside.append(f"{a.get('ID') or t}@{s}-{e} outside {p} {ps}-{pe}")
        if src == NEW_SOURCE and t == "CDS":
            cds_len[a.get("Parent")] += e - s + 1
    bad_cds = [f"{k} ({v} nt)" for k, v in cds_len.items() if v % 3]
    print("\nValidation of output:")
    print(f"  duplicate IDs introduced: {len(dup)}" + (f" ({', '.join(dup[:10])})" if dup else ""))
    print(f"  Parent references not found: {len(dangling)}" + (f" ({'; '.join(dangling[:5])})" if dangling else ""))
    print(f"  new features outside their parent's span: {len(outside)}" + (f" ({'; '.join(outside[:5])})" if outside else ""))
    print(f"  new CDSs not a multiple of 3: {len(bad_cds)}" + (f" ({', '.join(bad_cds[:5])})" if bad_cds else ""))


CHANGE_COLUMNS = ["lifted_gene", "lifted_transcript", "decision", "source", "status",
                  "target_gene", "new_ids", "removed_ids", "flags", "notes"]


def main():
    print(f"Reading {EXISTING_GFF} ...")
    ex = Existing(EXISTING_GFF)
    if ex.dup_ids:
        print(f"  note: {len(ex.dup_ids)} IDs already used by more than one feature in the existing file (left as is)")

    print(f"Reading {LIFTED_GFF} ...")
    lifted = Lifted(LIFTED_GFF)

    results = load_csv(RESULTS_CSV)
    manual = {}
    for r in load_csv(MANUAL_CSV):
        q, d = r.get("qry_id", "").strip(), r.get("decision", "").strip()
        if not q:
            continue
        if d not in KNOWN_DECISIONS:
            print(f"FLAG: unknown manual decision '{d}' for {q}; ignored")
            continue
        manual[q] = {"decision": d, "notes": r.get("notes", "")}
    print(f"Loaded {len(results)} result rows, {len(manual)} manual decisions")

    changes = []
    plans = build_plans(results, manual, lifted, changes.append)
    ed = Editor(ex, lifted, changes.append)
    # Removals first (Replace, ReplacePseudo), so new-gene placement and
    # ID allocation see the post-removal file.
    order = {"ReplacePseudo": 0, "Replace": 1, "Add": 2, "AddGene": 3, "AddPseudo": 4, "Leave": 5}
    for p in sorted(plans, key=lambda p: order.get(p["decision"], 9)):
        ed.apply(p)

    ed.write(OUT_GFF)
    with open(CHANGES_TSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=CHANGE_COLUMNS, delimiter="\t", lineterminator="\n",
                           restval="", extrasaction="ignore")
        w.writeheader()
        w.writerows(changes)

    tally = Counter((c.get("decision", ""), c.get("status", "")) for c in changes)
    print("\nChanges by decision and status:")
    for (d, s), n in sorted(tally.items()):
        print(f"  {d:<14} {s:<20} {n}")
    flagged = [c for c in changes if c.get("flags")]
    for c in flagged:
        print(f"FLAG: {c.get('lifted_gene') or c.get('lifted_transcript')} [{c.get('decision')}]: {c['flags']}")
    print(f"\nUpdated annotation written to: {OUT_GFF}")
    print(f"Change log written to: {CHANGES_TSV} ({len(changes)} records, {len(flagged)} flagged)")
    validate(OUT_GFF, ex)


if __name__ == "__main__":
    sys.exit(main())