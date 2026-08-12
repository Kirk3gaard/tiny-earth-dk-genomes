"""
Subset the GTDB-Tk bac120 user MSA to the genomes that passed assembly QC.

gtdbtk classify_wf already identifies the 120 bacterial marker genes in every
genome and produces a trimmed, concatenated alignment of them. That alignment is
a far better basis for a tree of a taxonomically broad isolate collection than a
de novo singleton search: the markers are curated and universal, so there is no
"no locus is present in enough genomes" failure mode. This script takes that
alignment, keeps only the genomes we are willing to defend, renames the records
so the tip labels match SEQID everywhere else in the pipeline, and records what
was kept and why.

Deliberately a plain command-line program invoked from a `shell:` directive
rather than a snakemake `script:` rule. The `script:` mechanism prepends a
preamble that extends sys.path into snakemake's own site-packages and unpickles
the snakemake object using the *rule's conda environment* Python. When that
Python differs from the one snakemake runs under — as it does for the pinned
gtdbtk environment — the preamble fails before any user code runs, so nothing
reaches the log and the failure is undiagnosable. A CLI has no such coupling,
and the calling shell block handles logging.

This module uses only the standard library, so it needs no conda environment.

Usage:
    python3 scripts/subset_msa.py \\
        --gtdbtk-dir temp/gtdbtk_out \\
        --assembly-qc results/assembly_qc_summary.tsv \\
        --out-msa temp/tree/bac120_subset.faa \\
        --out-membership results/tree_membership.tsv \\
        [--min-markers 0] [--genome-suffix .flye]
"""

import argparse
import glob
import gzip
import os
import sys


def describe_tree(root, limit=60):
    """Render the directory tree, so a 'file not found' says what IS present."""
    lines = []
    if not os.path.isdir(root):
        return [f"  {root} does not exist or is not a directory"]
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames.sort()
        rel = os.path.relpath(dirpath, root)
        prefix = "" if rel == "." else rel + "/"
        for name in sorted(filenames):
            lines.append(f"  {prefix}{name}")
            if len(lines) >= limit:
                lines.append(f"  ... (listing truncated at {limit} entries)")
                return lines
        if not filenames and not dirnames:
            lines.append(f"  {prefix}  (empty directory)")
    return lines or [f"  {root} is empty"]


# ── IO helpers ───────────────────────────────────────────────────────────────
def open_text(path):
    """Open plain or gzipped text transparently."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def read_fasta(path):
    """
    Yield (name, sequence) pairs. Tolerates wrapped sequence lines and header
    descriptions after the first whitespace token.
    """
    name, chunks = None, []
    with open_text(path) as handle:
        for line in handle:
            line = line.rstrip("\n").rstrip("\r")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].strip().split()[0] if line[1:].strip() else ""
                chunks = []
            elif line:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def find_one(gtdbtk_dir, pattern):
    """Locate a single file under gtdbtk_dir, searching recursively."""
    hits = sorted(glob.glob(os.path.join(gtdbtk_dir, "**", pattern), recursive=True))
    return hits[0] if hits else None


def read_tsv(path):
    """Read a TSV with a header row into a list of dicts."""
    rows = []
    with open_text(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            if not line.strip():
                continue
            rows.append(dict(zip(header, line.rstrip("\n").split("\t"))))
    return rows


# ── genome id mapping ────────────────────────────────────────────────────────
def to_seqid(genome_id, known, suffix):
    """
    Map a GTDB-Tk genome id onto the pipeline's SEQID.

    gtdbtk ran with --extension fa.gz over files named {sample}.flye.fa.gz, so
    it strips only ".fa.gz" and reports "{sample}.flye". Prefer an exact match,
    then fall back to stripping the configured suffix, so this keeps working if
    the assembly filenames ever change.
    """
    if genome_id in known:
        return genome_id
    if suffix and genome_id.endswith(suffix):
        stripped = genome_id[: -len(suffix)]
        if stripped in known:
            return stripped
        return stripped
    return genome_id


def main(gtdbtk_dir, asm_qc_path, out_msa, out_membership,
         min_markers=0, genome_suffix=".flye", min_coverage=0.8):

    # ── which genomes are we willing to put on a tree ────────────────────────
    asm_rows = read_tsv(asm_qc_path)
    passing = [r["sample"] for r in asm_rows if r.get("status") == "PASS"]
    passing_set = set(passing)
    asm_reason = {r["sample"]: r.get("reason", "-") for r in asm_rows}
    print(f"assembly QC: {len(passing)} of {len(asm_rows)} genomes passed", flush=True)

    # ── locate the bac120 alignment ──────────────────────────────────────────
    msa_path = (find_one(gtdbtk_dir, "gtdbtk.bac120.user_msa.fasta.gz")
                or find_one(gtdbtk_dir, "gtdbtk.bac120.user_msa.fasta"))
    if msa_path is None:
        print("ERROR: could not find gtdbtk.bac120.user_msa.fasta[.gz] under "
              f"{gtdbtk_dir}.", file=sys.stderr)
        print("Expected it at align/gtdbtk.bac120.user_msa.fasta.gz. Note that "
              "classify_wf removes align/intermediate_results but keeps align/ "
              "itself, so the MSA should survive a default run.", file=sys.stderr)
        print(f"What is actually under {gtdbtk_dir}:", file=sys.stderr)
        for line in describe_tree(gtdbtk_dir):
            print(line, file=sys.stderr)
        sys.exit(1)
    print(f"bac120 MSA: {msa_path}", flush=True)

    # Archaeal genomes get a separate marker set and cannot share a bac120 tree.
    ar_msa = (find_one(gtdbtk_dir, "gtdbtk.ar53.user_msa.fasta.gz")
              or find_one(gtdbtk_dir, "gtdbtk.ar53.user_msa.fasta"))
    if ar_msa:
        n_ar = sum(1 for _ in read_fasta(ar_msa))
        print(
            f"NOTE: {n_ar} genome(s) were classified as archaea and appear in "
            f"{os.path.basename(ar_msa)}. They use the ar53 marker set and are "
            "not included in this bacterial tree.",
            flush=True,
        )

    # ── per-genome bac120 marker counts, for reporting ───────────────────────
    markers = {}
    markers_path = find_one(gtdbtk_dir, "gtdbtk.bac120.markers_summary.tsv")
    if markers_path:
        for row in read_tsv(markers_path):
            name = row.get("name", "")
            seqid = to_seqid(name, passing_set, genome_suffix)
            try:
                unique = int(row.get("number_unique_genes", 0) or 0)
                multi_unique = int(row.get("number_multiple_unique_genes", 0) or 0)
                missing = int(row.get("number_missing_genes", 0) or 0)
            except ValueError:
                continue
            # GTDB-Tk counts "multiple unique" (identical duplicate hits) as
            # unique for downstream purposes, so include them.
            markers[seqid] = (unique + multi_unique, missing)
        print(f"marker counts: {markers_path}", flush=True)
    else:
        print("NOTE: markers_summary.tsv not found; marker counts unavailable.",
              flush=True)

    # ── read and filter the alignment ────────────────────────────────────────
    in_msa = {}
    for genome_id, seq in read_fasta(msa_path):
        in_msa[to_seqid(genome_id, passing_set, genome_suffix)] = seq
    covered = len(passing_set & set(in_msa))
    coverage = covered / len(passing_set) if passing_set else 0.0
    print(f"bac120 MSA contains {len(in_msa)} genomes; {covered} of "
          f"{len(passing_set)} QC-passing genomes are in it "
          f"({coverage:.1%} coverage)", flush=True)

    # A broken-input check, not a quality filter. If most QC-passing genomes are
    # missing from the alignment, the MSA was not built for this collection and
    # any tree from it would silently under-represent the data.
    if min_coverage > 0 and coverage < min_coverage:
        print("", file=sys.stderr)
        print(f"ERROR: only {coverage:.1%} of assembly-QC-passing genomes are "
              f"present in the bac120 MSA (threshold {min_coverage:.0%}).",
              file=sys.stderr)
        print("This usually means the MSA came from `gtdbtk classify_wf`, whose "
              "skani ANI pre-screen removes every genome it can classify to "
              "species level from the identify step, so those genomes never "
              "reach the alignment. Build the MSA with standalone `gtdbtk "
              "identify` + `gtdbtk align` instead (rules gtdbtk_identify and "
              "gtdbtk_align), which have no ANI screen.", file=sys.stderr)
        print(f"MSA searched: {msa_path}", file=sys.stderr)
        print("Set tree_min_msa_coverage: 0 in config.yaml to override.",
              file=sys.stderr)
        sys.exit(1)

    kept, membership = [], []
    for sample in sorted(passing_set):
        n_markers, n_missing = markers.get(sample, ("", ""))
        seq = in_msa.get(sample)
        reason = ""
        if seq is None:
            reason = "absent_from_bac120_msa"
        elif min_markers and isinstance(n_markers, int) and n_markers < min_markers:
            reason = f"bac120_markers<{min_markers}"
        if not reason:
            kept.append((sample, seq))
        membership.append({
            "sample": sample,
            "assembly_qc": "PASS",
            "assembly_qc_reason": asm_reason.get(sample, "-"),
            "in_bac120_msa": "TRUE" if seq is not None else "FALSE",
            "n_bac120_markers": str(n_markers),
            "n_bac120_missing": str(n_missing),
            "in_bac120_tree": "FALSE" if reason else "TRUE",
            "tree_exclusion_reason": reason or "-",
        })

    # Genomes that failed assembly QC never had a chance; record them too so the
    # table accounts for every genome, not just the ones that got close.
    for row in asm_rows:
        if row.get("status") == "PASS":
            continue
        membership.append({
            "sample": row["sample"],
            "assembly_qc": "FAIL",
            "assembly_qc_reason": row.get("reason", "-"),
            "in_bac120_msa": "NA",
            "n_bac120_markers": "",
            "n_bac120_missing": "",
            "in_bac120_tree": "FALSE",
            "tree_exclusion_reason": "failed_assembly_qc",
        })

    # ── report ───────────────────────────────────────────────────────────────
    dropped = [m for m in membership
               if m["assembly_qc"] == "PASS" and m["in_bac120_tree"] == "FALSE"]
    print(f"kept {len(kept)} genomes for the tree", flush=True)
    if dropped:
        print(f"{len(dropped)} assembly-QC-passing genome(s) still excluded:",
              flush=True)
        for m in dropped[:20]:
            print(f"    {m['sample']}: {m['tree_exclusion_reason']}", flush=True)
        if len(dropped) > 20:
            print(f"    ... and {len(dropped) - 20} more "
                  "(see results/tree_membership.tsv)", flush=True)

    if len(kept) < 4:
        sys.exit(
            f"ERROR: only {len(kept)} genomes available for tree inference; "
            "FastTree needs at least 4. See results/tree_membership.tsv."
        )

    # Every record in an alignment must be the same length.
    lengths = {len(seq) for _, seq in kept}
    if len(lengths) != 1:
        sys.exit(
            "ERROR: subset alignment is ragged; sequence lengths present: "
            f"{sorted(lengths)}. The bac120 user MSA should be uniform width."
        )
    print(f"alignment width: {lengths.pop()} columns", flush=True)

    # ── write ────────────────────────────────────────────────────────────────
    os.makedirs(os.path.dirname(out_msa) or ".", exist_ok=True)
    with open(out_msa, "w") as handle:
        for sample, seq in kept:
            handle.write(f">{sample}\n{seq}\n")

    os.makedirs(os.path.dirname(out_membership) or ".", exist_ok=True)
    columns = ["sample", "assembly_qc", "assembly_qc_reason", "in_bac120_msa",
               "n_bac120_markers", "n_bac120_missing", "in_bac120_tree",
               "tree_exclusion_reason"]
    with open(out_membership, "w") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in sorted(membership, key=lambda r: r["sample"]):
            handle.write("\t".join(row[c] for c in columns) + "\n")


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Subset the GTDB-Tk bac120 user MSA to QC-passing genomes.")
    p.add_argument("--gtdbtk-dir", required=True,
                   help="gtdbtk classify_wf output directory (temp/gtdbtk_out)")
    p.add_argument("--assembly-qc", required=True,
                   help="results/assembly_qc_summary.tsv")
    p.add_argument("--out-msa", required=True,
                   help="path to write the subset alignment (FASTA)")
    p.add_argument("--out-membership", required=True,
                   help="path to write the tree membership table (TSV)")
    p.add_argument("--min-markers", type=int, default=0,
                   help="minimum bac120 unique markers required; 0 disables")
    p.add_argument("--genome-suffix", default=".flye",
                   help="suffix gtdbtk appended to genome ids (default: .flye)")
    p.add_argument("--min-coverage", type=float, default=0.8,
                   help="fail if fewer than this fraction of QC-passing genomes "
                        "are in the MSA; 0 disables (default: 0.8)")
    return p.parse_args(argv)


if __name__ == "__main__":
    args = parse_args()
    print("subset_msa.py running under Python "
          + sys.version.split()[0] + " at " + sys.executable, flush=True)
    main(
        gtdbtk_dir=args.gtdbtk_dir,
        asm_qc_path=args.assembly_qc,
        out_msa=args.out_msa,
        out_membership=args.out_membership,
        min_markers=args.min_markers,
        genome_suffix=args.genome_suffix,
        min_coverage=args.min_coverage,
    )
