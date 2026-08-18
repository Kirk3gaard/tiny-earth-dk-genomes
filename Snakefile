configfile: "config.yaml"

import os
import re

# ── Input discovery ──────────────────────────────────────────────────────────
# Two entry points:
#   data/{sample}.fastq.gz                    raw nanopore reads -> filtlong -> flye
#   data/input_assemblies/{sample}.fa.gz      a finished assembly, reads unavailable
#   data/input_assemblies/{sample}.flye.fa.gz (the original naming, also accepted)
#
# Pre-made assemblies are imported into temp/all_assemblies/ next to flye's
# output and treated identically from checkm2 onwards — same QC gate, same
# annotation, same tree.
#
# NOTE ON THE DIRECTORY NAME: the import source used to be a top-level
# all_assemblies/, one "temp/" away from the pipeline's own output directory
# temp/all_assemblies/. Genomes dropped into the output directory by hand were
# never discovered (nothing globs it) while still being swept up by gtdbtk's
# --genome_dir, so they were classified but never assembled, QC'd, annotated or
# placed on the tree. The source is now data/input_assemblies/, which cannot be
# confused with an output path, and both the old location and stray files in
# temp/all_assemblies/ are warned about below.
(input_np_raw,) = glob_wildcards("data/{filename}.fastq.gz")

IMPORT_DIR = "data/input_assemblies"


def _discover_assemblies():
    """Map sample name -> path of a pre-made assembly under data/input_assemblies/."""
    sources = {}
    for name in glob_wildcards(IMPORT_DIR + "/{f}.flye.fa.gz").f:
        sources[name] = IMPORT_DIR + "/" + name + ".flye.fa.gz"
    for name in glob_wildcards(IMPORT_DIR + "/{f}.fa.gz").f:
        # {f}.fa.gz also matches the .flye.fa.gz files, leaving ".flye" on the
        # captured name; those are already recorded above.
        if name.endswith(".flye"):
            continue
        sources.setdefault(name, IMPORT_DIR + "/" + name + ".fa.gz")
    return sources


ASSEMBLY_SOURCES = _discover_assemblies()

# Reads win over a pre-made assembly for the same sample, so a stale file in
# data/input_assemblies/ cannot quietly shadow a genome the pipeline can assemble.
_read_based = set(input_np_raw)
_both = sorted(set(ASSEMBLY_SOURCES) & _read_based)
for _name in _both:
    del ASSEMBLY_SOURCES[_name]
if _both:
    print("NOTE: {} sample(s) have both reads and a pre-made assembly; "
          "assembling from reads: {}{}".format(
              len(_both), ", ".join(_both[:5]), " ..." if len(_both) > 5 else ""))

IMPORTED_SAMPLES = sorted(ASSEMBLY_SOURCES)
ALL_SAMPLES = sorted(_read_based | set(IMPORTED_SAMPLES))
print("Input: {} sample(s) with reads, {} imported assembly/assemblies, "
      "{} total.".format(len(input_np_raw), len(IMPORTED_SAMPLES), len(ALL_SAMPLES)))


# ── Misplaced-input warnings ─────────────────────────────────────────────────
# Both of these are the failure mode that motivated moving the import directory:
# a genome sitting somewhere the pipeline does not look, producing no error and
# no work. Neither aborts the run — they just make the situation impossible to
# miss in the log.
def _warn_misplaced_inputs():
    if os.path.isdir("all_assemblies"):
        legacy = sorted(f for f in os.listdir("all_assemblies") if f.endswith(".fa.gz"))
        if legacy:
            print("WARNING: {} assembly/assemblies found in the legacy directory "
                  "all_assemblies/ — this location is NO LONGER READ. Move them to "
                  "{}/ to have them processed: {}{}".format(
                      len(legacy), IMPORT_DIR, ", ".join(legacy[:5]),
                      " ..." if len(legacy) > 5 else ""))

    if os.path.isdir("temp/all_assemblies"):
        expected = set()
        for name in ALL_SAMPLES:
            expected.add(name + ".flye.fa.gz")
            expected.add(name + ".assembly_info.txt")
        stray = sorted(f for f in os.listdir("temp/all_assemblies") if f not in expected)
        if stray:
            print("WARNING: {} unrecognised file(s) in temp/all_assemblies/. That is a "
                  "pipeline OUTPUT directory, not an input one: files placed there are "
                  "picked up by gtdbtk's genome scan but by nothing else, so they are "
                  "classified while never being QC'd, annotated or placed on the tree. "
                  "Pre-made assemblies belong in {}/ as <sample>.fa.gz. Stray: {}{}".format(
                      len(stray), IMPORT_DIR, ", ".join(stray[:5]),
                      " ..." if len(stray) > 5 else ""))


_warn_misplaced_inputs()


def _only(names):
    """Regex matching exactly these sample names, or nothing when empty."""
    return "|".join(re.escape(n) for n in sorted(names)) if names else r"$^"


# 1 metadata file

# Output files
# genome assembly
# genome QC
# genome classification
# genome annotation
# antiobiotic profile
# tree with names from metadata

# ── QC-gated sample lists ────────────────────────────────────────────────────
# checkpoint read_qc (defined below, after filtlong_subset) inspects each
# sample's filtered reads and marks it PASS/FAIL. These helper functions
# resolve to the PASS-only sample list once read_qc has run for everyone,
# and are used as inputs wherever the pipeline previously used the static
# input_np_raw list for flye and everything downstream of it.
def passing_samples(wildcards):
    """
    Samples whose filtered reads passed the QC thresholds in read_qc.
    Calling checkpoints.read_qc.get() forces Snakemake to run read_qc for
    every raw sample first, then re-evaluates the DAG so flye and everything
    downstream of it only ever gets requested for samples worth assembling.
    """
    passed = []
    for name in input_np_raw:
        qc_file = checkpoints.read_qc.get(sample=name).output[0]
        with open(qc_file) as f:
            next(f)  # header
            fields = f.readline().strip().split("\t")
            if len(fields) == 5 and fields[4] == "PASS":
                passed.append(name)
    # Imported assemblies have no reads to QC, so they skip this gate and enter
    # at the assembly stage. They still face assembly_qc like everything else.
    passed.extend(IMPORTED_SAMPLES)
    return sorted(passed)


def flye_outputs(wildcards):
    return expand("temp/all_assemblies/{name}.flye.fa.gz", name=passing_samples(wildcards))

def flye_info_outputs(wildcards):
    return expand("temp/all_assemblies/{name}.assembly_info.txt", name=passing_samples(wildcards))

def checkm2_outputs(wildcards):
    return expand("temp/checkm2/{name}", name=passing_samples(wildcards))

def bakta_outputs(wildcards):
    return expand("temp/bakta_out/{name}", name=passing_samples(wildcards))

def antismash_outputs(wildcards):
    return expand("temp/antismash_out/{name}", name=passing_samples(wildcards))

def antismash_csv_outputs(wildcards):
    return expand("temp/antismash_csv/{name}_overview.csv", name=passing_samples(wildcards))


def assembly_qc_outputs(wildcards):
    return expand("temp/assembly_qc/{name}.asmqc.tsv", name=passing_samples(wildcards))


# ── Assembly-level QC gate (hard filter, tree only) ──────────────────────────
# checkpoint assembly_qc (defined below, after checkm2) combines each sample's
# CheckM2 completeness/contamination with its contig count and assembly size
# into a single PASS/FAIL verdict.
#
# This gate is deliberately NOT applied to the pipeline at large. Genomes that
# fail are still annotated and still contribute BGC results; they are only kept
# off the tree. Branch lengths and topology are only as trustworthy as the
# weakest genome contributing to the alignment, so the tree gets a clean input
# set. The BGC catalogue does not need one.
def tree_samples(wildcards):
    """
    Samples good enough to put on the tree: passed read QC (or were imported as
    finished assemblies), then passed the assembly QC thresholds in checkpoint
    assembly_qc. Calling checkpoints.assembly_qc.get() forces flye + checkm2 to
    complete for every candidate before the tree branch of the DAG is resolved.
    """
    passed = []
    for name in passing_samples(wildcards):
        qc_file = checkpoints.assembly_qc.get(sample=name).output[0]
        with open(qc_file) as f:
            next(f)  # header
            fields = f.readline().rstrip("\n").split("\t")
            if len(fields) == 7 and fields[5] == "PASS":
                passed.append(name)
    return passed


def tree_assembly_inputs(wildcards):
    return expand("temp/all_assemblies/{name}.flye.fa.gz", name=tree_samples(wildcards))


# ── Tree ─────────────────────────────────────────────────────────────────────
# Built from the GTDB-Tk bac120 marker alignment: stage_tree_genomes ->
# gtdbtk_identify -> gtdbtk_align -> subset_bac120_msa -> infer_tree_gtdbtk.
#
# This replaced getphylo, which was removed on 2026-08-07. getphylo searched for
# single-copy orthologues de novo and required each locus to be present in a
# high percentage of genomes; across a taxonomically broad soil isolate
# collection nothing survived that threshold, and it never produced a tree here
# (three runs, all InsufficientLociError). bac120 is a curated universal marker
# set, so that failure mode cannot occur, and classify_wf has already done the
# expensive part.
TREE_FILE = "temp/tree/gtdbtk_tree.nwk"


rule all:
    input:
        expand("temp/filtered_reads/{name}.filtlong.fastq.gz", name=input_np_raw),
        expand("temp/qc/{name}.qc.tsv", name=input_np_raw),
        "results/qc_summary.tsv",
        flye_outputs,
        checkm2_outputs,
        "results/assembly_qc_summary.tsv",
        "temp/gtdbtk_out",
        bakta_outputs,
        antismash_outputs,
        antismash_csv_outputs,
        TREE_FILE,
        "results/tree_membership.tsv",
        "results/antismash_failures.tsv",
        "results/aggregated_results.tsv",
        "results/phylogenetic_tree.pdf"

rule filtlong_subset:
    input:
        "data/{sample}.fastq.gz"
    output:
        "temp/filtered_reads/{sample}.filtlong.fastq.gz"
    threads: 4
    resources:
        mem_mb=lambda wc, input: max(3 * input.size_mb, 4096),
        node_type="general",
        time="00-05:00:00",
    conda:
        "envs/filtlong_env.yml"
    params:
        bases=config.get('filtlong_target', 300000000)
    shell:
        """
        filtlong --target_bases {params.bases} {input} | gzip > {output}
        """


checkpoint read_qc:
    # Evaluates the filtered reads for each sample and flags samples whose
    # filtered read set is too small/short to realistically assemble.
    # This stops near-empty barcodes or ones dominated by short adapter/
    # primer-dimer reads from consuming a flye/checkm2/bakta/antismash slot
    # and then aborting the whole run once nothing else can proceed.
    input:
        "temp/filtered_reads/{sample}.filtlong.fastq.gz"
    output:
        "temp/qc/{sample}.qc.tsv"
    threads: 1
    resources:
        mem_mb=2048,
        node_type="general",
        time="00-00:30:00",
    params:
        min_bases=config.get("qc_min_total_bases", 10000000),
        min_mean_len=config.get("qc_min_mean_length", 300),
    shell:
        """
        mkdir -p temp/qc
        zcat {input} | awk -v min_bases={params.min_bases} -v min_mean={params.min_mean_len} -v sample={wildcards.sample} '
            NR % 4 == 2 {{ c++; b += length($0) }}
            END {{
                mean = (c > 0) ? b / c : 0
                status = (b >= min_bases && mean >= min_mean) ? "PASS" : "FAIL"
                printf "sample\\treads\\ttotal_bp\\tmean_length\\tstatus\\n"
                printf "%s\\t%d\\t%d\\t%.2f\\t%s\\n", sample, c, b, mean, status
            }}' > {output}
        """


rule qc_summary:
    # Combines every sample's QC verdict into one table so skipped samples
    # are visible at a glance instead of only surfacing as a cryptic flye
    # error buried in the cluster log.
    input:
        expand("temp/qc/{name}.qc.tsv", name=input_np_raw)
    output:
        "results/qc_summary.tsv"
    threads: 1
    resources:
        mem_mb=1024,
        node_type="general",
        time="00-00:10:00",
    shell:
        """
        mkdir -p results
        head -n1 {input[0]} > {output}
        for f in {input}; do tail -n +2 "$f" >> {output}; done
        """


rule import_assembly:
    # Bring a pre-made assembly from data/input_assemblies/ into
    # temp/all_assemblies/ so it is indistinguishable from flye output
    # downstream. Used for genomes whose reads are no longer available.
    #
    # flye also writes assembly_info.txt, which assembly_qc reads for contig
    # count and total length and compile_results reads for genome size. There is
    # no such file for an imported assembly, so an equivalent one is derived
    # from the FASTA itself. Coverage and circularity are genuinely unknown and
    # are recorded as NA rather than invented.
    input:
        asm = lambda wc: ASSEMBLY_SOURCES[wc.sample]
    output:
        asm     = "temp/all_assemblies/{sample}.flye.fa.gz",
        asminfo = "temp/all_assemblies/{sample}.assembly_info.txt"
    wildcard_constraints:
        sample = _only(IMPORTED_SAMPLES)
    threads: 1
    resources:
        mem_mb=2048,
        node_type="general",
        time="00-00:20:00",
    log:
        "logs/import_assembly/{sample}.log"
    shell:
        """
        mkdir -p temp/all_assemblies logs/import_assembly
        exec > >(tee -a "{log}") 2>&1
        echo "=== import_assembly {wildcards.sample} $(date -Is) ==="
        echo "source: {input.asm}"

        cp "{input.asm}" "{output.asm}"

        zcat "{output.asm}" | awk -v OFS='\\t' '
            /^>/ {{
                if (n != "") print n, L, "NA", "N", "N", 1, "*", "*"
                split(substr($0, 2), a, /[ \\t]/); n = a[1]; L = 0; next
            }}
            {{ L += length($0) }}
            END {{
                if (n != "") print n, L, "NA", "N", "N", 1, "*", "*"
            }}
        ' > "{output.asminfo}.body"

        printf '#seq_name\\tlength\\tcov.\\tcirc.\\trepeat\\tmult.\\talt_group\\tgraph_path\\n' \
            > "{output.asminfo}"
        cat "{output.asminfo}.body" >> "{output.asminfo}"
        rm -f "{output.asminfo}.body"

        NCONTIG=$(($(wc -l < "{output.asminfo}") - 1))
        TOTAL=$(awk -F'\\t' 'NR>1 {{s+=$2}} END {{print s+0}}' "{output.asminfo}")
        echo "imported: $NCONTIG contigs, $TOTAL bp"
        if [ "$NCONTIG" -lt 1 ] || [ "$TOTAL" -lt 1 ]; then
            echo "ERROR: no sequence found in {input.asm}" >&2
            exit 1
        fi
        """


rule flye:
    input:
        NPreads="temp/filtered_reads/{sample}.filtlong.fastq.gz"
    output:
        asm="temp/all_assemblies/{sample}.flye.fa.gz",
        asminfo="temp/all_assemblies/{sample}.assembly_info.txt"
    wildcard_constraints:
        # Keeps flye and import_assembly from both matching the same target.
        sample = _only(_read_based)
    threads: config["assembly_threads"]
    resources:
        mem_mb=config["assembly_mb"],
        node_type="general",
        time="00-05:00:00",
    conda:
        "envs/env_flye.yml"
    shell:
        """
        mkdir -p temp/flye/{wildcards.sample}
        flye --nano-hq {input.NPreads} --threads $(nproc) --meta --out-dir temp/flye/{wildcards.sample}
        cat temp/flye/{wildcards.sample}/assembly.fasta | gzip > {output.asm}
        cp temp/flye/{wildcards.sample}/assembly_info.txt {output.asminfo}
        """

rule checkm2:
    input:
        asm="temp/all_assemblies/{sample}.flye.fa.gz",
    output:
        directory("temp/checkm2/{sample}")
    threads: 10
    resources:
        mem_mb=lambda wc, input: max(5 * input.size_mb, 10240),
        node_type="general",
        time="00-05:00:00",
    params:
        db=config.get("checkm2_db")
    conda:
        "envs/env_checkm2.yml"
    shell:
        """
        TMP_DIR=$(mktemp -d -t checkm2_{wildcards.sample}_XXXXXX)
        
        # Ensure it gets cleaned up on exit, even if the script fails
        trap "rm -rf '$TMP_DIR'" EXIT
        
        # 3. Copy the file into the temporary directory
        cp {input.asm} "$TMP_DIR"/
        checkm2 predict \
            --threads {threads} \
            --database_path {params.db} \
            -x fa.gz \
            --input "$TMP_DIR"/ \
            --output-directory {output}
        """

checkpoint assembly_qc:
    # Combines CheckM2 completeness/contamination with the contig count and
    # total assembly length into one PASS/FAIL verdict per sample.
    # PASS samples are the only ones that reach the tree (see tree_samples).
    #
    # Both inputs are already produced by the pipeline, so this adds no new
    # tools and costs seconds per sample. Reading contig count and size from
    # flye rather than CheckM2 keeps us on a column layout we control.
    input:
        checkm2    = "temp/checkm2/{sample}",
        flye_info  = "temp/all_assemblies/{sample}.assembly_info.txt"
    output:
        "temp/assembly_qc/{sample}.asmqc.tsv"
    threads: 1
    resources:
        mem_mb=1024,
        node_type="general",
        time="00-00:10:00",
    params:
        min_completeness  = config.get("asm_min_completeness", 90),
        max_contamination = config.get("asm_max_contamination", 5),
        max_contigs       = config.get("asm_max_contigs", 200),
        min_genome_bp     = config.get("asm_min_genome_bp", 1500000),
        max_genome_bp     = config.get("asm_max_genome_bp", 15000000),
    shell:
        """
        mkdir -p temp/assembly_qc

        # CheckM2 can exit cleanly but write an empty report for degenerate
        # assemblies. Substituting /dev/null makes that a FAIL with an explicit
        # reason rather than an awk crash under bash strict mode.
        CHECKM2_TSV="{input.checkm2}/quality_report.tsv"
        if [ ! -s "$CHECKM2_TSV" ]; then CHECKM2_TSV=/dev/null; fi

        awk -F'\\t' -v OFS='\\t' \
            -v sample="{wildcards.sample}" \
            -v min_comp={params.min_completeness} \
            -v max_cont={params.max_contamination} \
            -v max_contigs={params.max_contigs} \
            -v min_bp={params.min_genome_bp} \
            -v max_bp={params.max_genome_bp} '
            # flye assembly_info.txt: one row per contig, $2 = length
            FILENAME ~ /assembly_info/ {{
                if (FNR == 1) next
                contigs++; total += $2; next
            }}
            # checkm2 quality_report.tsv: header + one data row
            FILENAME ~ /quality_report/ {{
                if (FNR == 1) {{
                    for (i = 1; i <= NF; i++) {{
                        if ($i == "Completeness")  ci = i
                        if ($i == "Contamination") ni = i
                    }}
                    next
                }}
                if (FNR == 2 && ci && ni) {{
                    comp = $ci + 0; cont = $ni + 0; have_checkm2 = 1
                }}
                next
            }}
            END {{
                reason = ""
                if (!have_checkm2)                       reason = reason "no_checkm2_report;"
                if (have_checkm2 && comp <  min_comp)    reason = reason "completeness<" min_comp ";"
                if (have_checkm2 && cont >= max_cont)    reason = reason "contamination>=" max_cont ";"
                if (contigs > max_contigs)               reason = reason "contigs>" max_contigs ";"
                if (total   < min_bp)                    reason = reason "size<" min_bp ";"
                if (total   > max_bp)                    reason = reason "size>" max_bp ";"
                status = (reason == "") ? "PASS" : "FAIL"
                if (reason == "") reason = "-"
                print "sample","completeness","contamination","contigs","total_bp","status","reason"
                printf "%s\\t%.2f\\t%.2f\\t%d\\t%d\\t%s\\t%s\\n", \
                       sample, comp, cont, contigs, total, status, reason
            }}' "{input.flye_info}" "$CHECKM2_TSV" > {output}
        """


rule assembly_qc_summary:
    # One table of every assembly's QC verdict. Serves two purposes: a
    # human-readable record of what was excluded and why, and an input to both
    # subset_bac120_msa (which genomes may go on the tree) and
    # compile_results.R (so exclusions are documented alongside the results).
    input:
        assembly_qc_outputs
    output:
        "results/assembly_qc_summary.tsv"
    threads: 1
    resources:
        mem_mb=1024,
        node_type="general",
        time="00-00:10:00",
    shell:
        """
        mkdir -p results
        head -n1 {input[0]} > {output}
        for f in {input}; do tail -n +2 "$f" >> {output}; done
        """


rule gtdb:
    # Taxonomic classification of every genome the pipeline is carrying —
    # assembled and imported alike.
    #
    # Two things this rule used to get wrong, both of which mattered once
    # imported genomes entered the picture:
    #
    #  1. It pointed --genome_dir straight at temp/all_assemblies/, so gtdbtk
    #     classified whatever happened to be in that directory rather than the
    #     genomes snakemake had actually declared as inputs. Leftovers from a
    #     killed run, and hand-dropped files, were classified while being absent
    #     from every other stage. The declared inputs are now staged into
    #     temp/gtdbtk_in/ as symlinks — the same pattern stage_tree_genomes
    #     already uses — so {input} and --genome_dir cannot disagree.
    #
    #  2. It never cleared its own output directory. Adding a genome makes this
    #     rule rerun, and classify_wf refuses a non-empty --out_dir, so the
    #     rerun died on a directory the previous successful run had created.
    input:
        assemblies = flye_outputs,
    output:
        directory("temp/gtdbtk_out")
    threads: 32
    resources:
        mem_mb=250000,
        node_type="highmem",
        time="00-05:00:00",
        tmpdir="/tmp",
    conda:
        "envs/env_gtdbtk.yml"
    params:
        db=config.get("gtdb_db")
    log:
        "logs/gtdb.log"
    shell:
        """
        mkdir -p logs temp
        exec > >(tee -a "{log}") 2>&1
        echo "=== gtdb classify_wf started $(date -Is) ==="

        export GTDBTK_DATA_PATH="{params.db}"

        rm -rf temp/gtdbtk_in {output}
        mkdir -p temp/gtdbtk_in
        for f in {input.assemblies}; do
            ln -s "$(readlink -f "$f")" "temp/gtdbtk_in/$(basename "$f")"
        done

        N_GENOMES=$(ls temp/gtdbtk_in | wc -l)
        echo "classifying $N_GENOMES genomes"
        if [ "$N_GENOMES" -lt 1 ]; then
            echo "ERROR: no genomes staged for classification." >&2
            exit 1
        fi

        gtdbtk classify_wf \
            --genome_dir temp/gtdbtk_in \
            --out_dir {output} \
            --cpus {threads} \
            --extension fa.gz

        rm -rf temp/gtdbtk_in
        echo "=== gtdb classify_wf finished $(date -Is) ==="
        """

# rule prokka was removed on 2026-08-07. Its only consumer was getphylo, which
# is also gone; bakta provides the annotation the pipeline actually uses (the
# .gbff feeding antiSMASH, and the rRNA/tRNA counts in compile_results). It had
# continued to run on every genome because rule all requested it directly —
# 549 jobs at 8 threads, ~515 CPU-hours per run, producing output nothing read.
rule bakta:
    input:
        asm="temp/all_assemblies/{sample}.flye.fa.gz",
    output:
        directory("temp/bakta_out/{sample}")
    threads: 16
    resources:
        mem_mb=lambda wc, input: max(5 * input.size_mb, 46080),
        disk_mb=1000,
        node_type="general",
        time="00-05:00:00",
    params:
        db=config.get("bakta_db")
    conda:
        "envs/env_bakta.yaml"
    shell:
        """
        gzip -dc {input.asm} > "temp/{wildcards.sample}_temp.fa"
        
        bakta --db {params.db} \
              --output {output} \
              --prefix {wildcards.sample} \
              --threads {threads} \
              --tmp-dir "{resources.tmpdir}" \
               "temp/{wildcards.sample}_temp.fa"
        
        rm "temp/{wildcards.sample}_temp.fa"
        """

rule antismash:
    input:
        rules.bakta.output
    output:
        directory("temp/antismash_out/{sample}")
    threads: 16
    resources:
        mem_mb=lambda wc, input: max(5 * input.size_mb, 20240),
        node_type="general",
        time="00-05:00:00",
    params:
        db=config.get("antismash_db")
    log:
        "logs/antismash/{sample}.log"
    conda:
        "envs/env_antismash.yaml"
    shell:
        """
        mkdir -p logs/antismash {output}
        exec > >(tee -a "{log}") 2>&1
        echo "=== antismash {wildcards.sample} started $(date -Is) ==="

        # Record where this ran. On 2026-08-05, 526 of 549 jobs died instantly
        # with SIGILL because MOODS (a pip dependency of antiSMASH) hardcodes
        # -march=native in its setup.py and had been compiled on an AVX-512
        # machine; jobs landing on nodes without AVX-512 executed an illegal
        # instruction on import. Diagnosing that took days because nothing
        # recorded which node a job ran on, or whether the CPU could run the
        # binaries. Two lines here would have made it obvious.
        HOST=$(hostname)
        CPU=$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | sed 's/^ *//')
        HAS_AVX512=$(grep -qm1 avx512f /proc/cpuinfo && echo yes || echo no)
        echo "host=$HOST"
        echo "cpu=$CPU"
        echo "avx512f=$HAS_AVX512"

        # Fail fast and unmistakably if the interpreter cannot even load the
        # extension modules on this node, rather than burning a 16-cpu slot.
        if ! python -c "import MOODS.scan" 2>/dev/null; then
            echo "ERROR: 'import MOODS.scan' crashes on $HOST (avx512f=$HAS_AVX512)." >&2
            echo "The antismash conda env contains binaries this CPU cannot run." >&2
            echo "Rebuild MOODS on the oldest available architecture:" >&2
            echo "  pip install --force-reinstall --no-binary :all: --no-cache-dir MOODS-python==<version>" >&2
            printf 'MOODS import crashed on %s (avx512f=%s) at %s\\n' \
                "$HOST" "$HAS_AVX512" "$(date -Is)" > "{output}/ANTISMASH_FAILED"
            exit 0
        fi

        # antiSMASH failures are recorded rather than fatal. A single genome
        # that antiSMASH cannot parse should not block a 550-genome run: two
        # genomes have consistently died with
        #   ValueError: Cannot convert protein positions without protein length
        #               when sequence start coordinate is ambiguous
        # which is antiSMASH choking on CDS features with fuzzy coordinates
        # (e.g. <1..500) that bakta emits at contig edges. Genomes that fail
        # get an ANTISMASH_FAILED marker in their output directory; downstream
        # rules treat them as having no BGC data and compile_results reports
        # them, so they are visible rather than silently missing.
        set +e
        antismash --databases {params.db} \
        --output-dir {output} \
        --asf \
        --cc-mibig \
        --cb-subclusters \
        --cb-knownclusters \
        --rre \
        --cb-general \
        --tfbs \
        --pfam2go \
        -c {threads} \
        --genefinding-tool prodigal \
        "{input}/{wildcards.sample}.gbff"
        STATUS=$?
        set -e

        if [ "$STATUS" -ne 0 ]; then
            echo "WARNING: antismash exited $STATUS for {wildcards.sample}; recording as failed." >&2
            printf 'antismash exited %s at %s\\n' "$STATUS" "$(date -Is)" \
                > "{output}/ANTISMASH_FAILED"
        elif [ ! -f "{output}/regions.js" ] && [ ! -f "{output}/index.html" ]; then
            echo "WARNING: antismash exited 0 for {wildcards.sample} but produced no result files; recording as failed." >&2
            printf 'antismash produced no output at %s\\n' "$(date -Is)" \
                > "{output}/ANTISMASH_FAILED"
        fi

        echo "=== antismash {wildcards.sample} finished $(date -Is) ==="
        """


rule antismash_failures:
    # One table of the genomes antiSMASH could not process, so failures are
    # visible in results/ rather than only in the cluster logs.
    input:
        antismash_outputs
    output:
        "results/antismash_failures.tsv"
    threads: 1
    resources:
        mem_mb=1024,
        node_type="general",
        time="00-00:10:00",
    params:
        max_frac = config.get("antismash_max_failure_fraction", 0.1)
    log:
        "logs/antismash_failures.log"
    shell:
        """
        mkdir -p results logs
        exec > >(tee -a "{log}") 2>&1
        printf 'sample\\treason\\n' > {output}
        TOTAL=0
        for d in {input}; do
            TOTAL=$((TOTAL + 1))
            if [ -f "$d/ANTISMASH_FAILED" ]; then
                printf '%s\\t%s\\n' "$(basename "$d")" \
                    "$(tr -d '\\n' < "$d/ANTISMASH_FAILED")" >> {output}
            fi
        done
        FAILED=$(($(wc -l < {output}) - 1))
        echo "antismash: $FAILED of $TOTAL genomes failed"

        # Individual failures are tolerated so one unparseable genome cannot
        # block the run. A large fraction failing is a different thing entirely
        # — an environment or cluster problem — and must not pass as success.
        # On 2026-08-05, 526 of 549 died with exit 132 (SIGILL) and the run
        # still reported completion, leaving 17 of 376 tree genomes with BGC
        # data. Tune with antismash_max_failure_fraction in config.yaml.
        if [ "$TOTAL" -gt 0 ]; then
            awk -v f="$FAILED" -v t="$TOTAL" -v m={params.max_frac} '
                BEGIN {{
                    frac = f / t
                    printf "failure fraction: %.1f%% (threshold %.0f%%)\\n", frac*100, m*100
                    exit (frac > m) ? 1 : 0
                }}' || {{
                echo "ERROR: too many antismash failures. See {output} for the list" >&2
                echo "and logs/antismash/<sample>.log for individual errors." >&2
                exit 1
            }}
        fi
        """


rule antismash_processing:
    input:
        data=rules.antismash.output
    output:
        overview_path      = "temp/antismash_csv/{sample}_overview.csv",
        cluster_blast_path = "temp/antismash_csv/{sample}_cblast_general.csv",
        known_cblast_path  = "temp/antismash_csv/{sample}_cblast_known.csv",
        areas_path         = "temp/antismash_csv/{sample}_areas.csv"
    threads: 1
    resources:
        mem_mb=4096,
        node_type="general",
        time="00-01:00:00",
    conda:
        "envs/env_antismash_processing.yaml"
    script:
        "scripts/process_antismash.R"

# ── bac120 alignment for the tree ────────────────────────────────────────────
# NOT reused from rule gtdb. `classify_wf` runs a skani ANI pre-screen and then
# REMOVES every genome that got a species-level ANI match from the identify
# step (gtdbtk/main.py: identify(..., process_classified_species=False) drops
# them). Those genomes therefore never enter align/ and never appear in
# gtdbtk.bac120.user_msa.fasta.gz. On the 2026-08-04 run that left 52 of 549
# genomes in the MSA and a 42-tip tree.
#
# identify + align run standalone have no ANI screen, so every genome gets its
# markers found and aligned. This also skips pplacer entirely — the memory-hungry
# part — and leaves the existing classification in temp/gtdbtk_out untouched.
rule stage_tree_genomes:
    # gtdbtk needs a directory of genomes. Symlink just the assembly-QC-passing
    # ones, keeping the original {sample}.flye.fa.gz names so the genome ids
    # gtdbtk reports match what subset_msa.py expects.
    input:
        tree_assembly_inputs
    output:
        directory("temp/tree/genomes")
    threads: 1
    resources:
        mem_mb=1024,
        node_type="general",
        time="00-00:20:00",
    log:
        "logs/stage_tree_genomes.log"
    shell:
        """
        mkdir -p logs
        exec > >(tee -a "{log}") 2>&1
        rm -rf {output}
        mkdir -p {output}
        for f in {input}; do
            ln -s "$(readlink -f "$f")" "{output}/$(basename "$f")"
        done
        echo "staged $(ls {output} | wc -l) genomes for the bac120 alignment"
        """


rule gtdbtk_identify:
    input:
        genomes = "temp/tree/genomes"
    output:
        directory("temp/tree/gtdbtk_msa/identify")
    threads: 32
    resources:
        mem_mb=64000,
        node_type="general",
        time="01-00:00:00",
        tmpdir="/tmp",
    params:
        db = config.get("gtdb_db")
    log:
        "logs/gtdbtk_identify.log"
    conda:
        "envs/env_gtdbtk.yml"
    shell:
        """
        mkdir -p logs
        exec > >(tee -a "{log}") 2>&1
        echo "=== gtdbtk_identify started $(date -Is) ==="
        export GTDBTK_DATA_PATH="{params.db}"

        echo "genomes to process: $(ls {input.genomes} | wc -l)"
        rm -rf temp/tree/gtdbtk_msa/identify
        gtdbtk identify \
            --genome_dir {input.genomes} \
            --out_dir temp/tree/gtdbtk_msa \
            --extension fa.gz \
            --cpus {threads}
        echo "=== gtdbtk_identify finished $(date -Is) ==="
        """


rule gtdbtk_align:
    input:
        identify = "temp/tree/gtdbtk_msa/identify"
    output:
        directory("temp/tree/gtdbtk_msa/align")
    threads: 16
    resources:
        mem_mb=64000,
        node_type="general",
        time="00-12:00:00",
        tmpdir="/tmp",
    params:
        db = config.get("gtdb_db")
    log:
        "logs/gtdbtk_align.log"
    conda:
        "envs/env_gtdbtk.yml"
    shell:
        """
        mkdir -p logs
        exec > >(tee -a "{log}") 2>&1
        echo "=== gtdbtk_align started $(date -Is) ==="
        export GTDBTK_DATA_PATH="{params.db}"

        # --skip_gtdb_refs: we want an alignment of our isolates only, not the
        # thousands of GTDB reference genomes. user_msa is unaffected either
        # way, but skipping them makes this far lighter.
        rm -rf temp/tree/gtdbtk_msa/align
        gtdbtk align \
            --identify_dir temp/tree/gtdbtk_msa \
            --out_dir temp/tree/gtdbtk_msa \
            --skip_gtdb_refs \
            --cpus {threads}
        echo "=== gtdbtk_align finished $(date -Is) ==="
        """


rule subset_bac120_msa:
    # Subset the bac120 alignment to the assembly-QC-passing genomes and rename
    # records to bare SEQID so tip labels line up with aggregated_results.tsv.
    input:
        align_dir    = "temp/tree/gtdbtk_msa/align",
        identify_dir = "temp/tree/gtdbtk_msa/identify",
        asm_qc       = "results/assembly_qc_summary.tsv"
    output:
        msa        = "temp/tree/bac120_subset.faa",
        membership = "results/tree_membership.tsv"
    threads: 1
    resources:
        mem_mb=8192,
        node_type="general",
        time="00-00:30:00",
    params:
        gtdbtk_dir    = "temp/tree/gtdbtk_msa",
        min_markers   = config.get("tree_min_bac120_markers", 0),
        # gtdbtk ran with --extension fa.gz over {sample}.flye.fa.gz, so it
        # reports genome ids as "{sample}.flye". Strip that to recover SEQID.
        genome_suffix = ".flye",
        # Fail if the alignment covers far less of the QC-passing set than
        # expected. This is a broken-input check, not a quality filter: the
        # 2026-08-04 run silently produced a 42-tip tree from a 376-genome
        # collection because classify_wf's ANI screen had excluded almost
        # everything from the MSA, and nothing objected.
        min_coverage  = config.get("tree_min_msa_coverage", 0.8),
    log:
        "logs/subset_bac120_msa.log"
    # Invoked as a plain command rather than via `script:`, and with no conda
    # environment. The `script:` mechanism extends sys.path into snakemake's
    # site-packages and unpickles the snakemake object using the rule's conda
    # Python; when that Python differs from snakemake's own — as it does for the
    # pinned gtdbtk env — it fails in the preamble before any user code runs, so
    # the log stays empty and the failure cannot be diagnosed. subset_msa.py
    # needs only the standard library, so the interpreter running snakemake
    # (guaranteed present, since the jobscript invokes snakemake) is sufficient.
    shell:
        """
        mkdir -p logs temp/tree results
        exec > >(tee -a "{log}") 2>&1
        echo "=== subset_bac120_msa started $(date -Is) ==="
        echo "interpreter: $(command -v python3 || echo 'python3 NOT FOUND')"

        python3 scripts/subset_msa.py \
            --gtdbtk-dir "{params.gtdbtk_dir}" \
            --assembly-qc "{input.asm_qc}" \
            --out-msa "{output.msa}" \
            --out-membership "{output.membership}" \
            --min-markers {params.min_markers} \
            --genome-suffix "{params.genome_suffix}" \
            --min-coverage {params.min_coverage}

        echo "=== subset_bac120_msa finished $(date -Is) ==="
        """


rule infer_tree_gtdbtk:
    # FastTree (WAG + gamma, SH-like support) over the bac120 subset via
    # `gtdbtk infer`. Minutes rather than hours, because the expensive part —
    # marker identification and alignment — was already done by classify_wf.
    input:
        msa = "temp/tree/bac120_subset.faa"
    output:
        tree = "temp/tree/gtdbtk_tree.nwk"
    threads: 16
    resources:
        mem_mb=32768,
        node_type="general",
        time="00-08:00:00",
        tmpdir="/tmp",
    params:
        db        = config.get("gtdb_db"),
        prot_model = config.get("tree_prot_model", "WAG"),
    log:
        "logs/infer_tree_gtdbtk.log"
    conda:
        "envs/env_gtdbtk.yml"
    shell:
        """
        mkdir -p logs temp/tree
        exec > >(tee -a "{log}") 2>&1
        echo "=== infer_tree_gtdbtk started $(date -Is) ==="

        export GTDBTK_DATA_PATH="{params.db}"

        N_SEQS=$(grep -c '^>' {input.msa})
        echo "inferring tree from $N_SEQS genomes"

        rm -rf temp/tree/infer
        gtdbtk infer \
            --msa_file {input.msa} \
            --out_dir temp/tree/infer \
            --prot_model {params.prot_model} \
            --gamma \
            --cpus {threads} \
            --prefix gtdbtk

        # `gtdbtk infer` writes [prefix].unrooted.tree, but its exact depth in
        # the output directory has moved between versions. Locate it rather than
        # hardcoding a path, and fail loudly if there is not exactly one.
        mapfile -t TREES < <(find temp/tree/infer -name '*.unrooted.tree' -type f)
        echo "found ${{#TREES[@]}} candidate tree file(s): ${{TREES[*]:-none}}"
        if [ "${{#TREES[@]}}" -ne 1 ]; then
            echo "ERROR: expected exactly one *.unrooted.tree under temp/tree/infer." >&2
            exit 1
        fi
        cp "${{TREES[0]}}" {output.tree}

        echo "tip count in final tree: $(tr ',' '\\n' < {output.tree} | grep -c ':' || true)"
        echo "=== infer_tree_gtdbtk finished $(date -Is) ==="
        """


# rule getphylo was removed on 2026-08-07, replaced by the bac120 tree above.
# It searched for single-copy orthologues de novo and required each locus to be
# present in a high percentage of the input genomes. Across a taxonomically
# broad soil-isolate collection nothing cleared that bar: every run ended in
# InsufficientLociError, even after hard-filtering the input set and lowering
# the presence threshold to 95%. It never produced a tree in this project.
# Recover it from git history if a de novo marker search is ever wanted.

rule compile_results:
    input:
        checkm2_reports   = checkm2_outputs,
        flye_info         = flye_info_outputs,
        gtdbtk_summary    = "temp/gtdbtk_out",
        bakta_dirs        = bakta_outputs,
        antismash_overviews = antismash_csv_outputs,
        qc_summary        = "results/qc_summary.tsv",
        asm_qc_summary    = "results/assembly_qc_summary.tsv",
        tree_membership   = "results/tree_membership.tsv",
        antismash_fails   = "results/antismash_failures.tsv"
    output:
        tsv     = "results/aggregated_results.tsv",
        # Written by compile_results.R alongside the main table. Previously
        # undeclared, which meant nothing in the DAG knew they existed.
        mibig   = "results/antismash_mibig_summary.tsv",
        product = "results/antismash_product_summary.tsv"
    resources:
        mem_mb=8192,
        node_type="general",
        time="00-01:00:00",
    log:
        "logs/compile_results.log"
    conda:
        "envs/R-main.yaml"
    script:
        "scripts/compile_results.R"

rule plot_tree:
    input:
        tree     = TREE_FILE,
        metadata = "results/aggregated_results.tsv",
        # Produced as side outputs of compile_results; declared so the script
        # reads them from snakemake@input rather than hardcoded paths.
        mibig    = "results/antismash_mibig_summary.tsv",
        product  = "results/antismash_product_summary.tsv"
    output:
        plot = "results/phylogenetic_tree.pdf"
    resources:
        mem_mb=8192,
        node_type="general",
        time="00-01:00:00",
    log:
        "logs/plot_tree.log"
    conda:
        "envs/R-main.yaml"
    script:
        "scripts/plot_tree.R"
