configfile: "config.yaml"

# Input files
(input_np_raw,) = glob_wildcards("data/{filename}.fastq.gz")
(input_np_asm,) = glob_wildcards("all_assemblies/{filename}.flye.fa.gz")

# input files
# raw nanopore data 1 fastq.gz file per sample
# alternative input 1 assembly per sample
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
    return passed


def flye_outputs(wildcards):
    return expand("temp/all_assemblies/{name}.flye.fa.gz", name=passing_samples(wildcards))

def flye_info_outputs(wildcards):
    return expand("temp/all_assemblies/{name}.assembly_info.txt", name=passing_samples(wildcards))

def checkm2_outputs(wildcards):
    return expand("temp/checkm2/{name}", name=passing_samples(wildcards))

def prokka_outputs(wildcards):
    return expand("temp/prokka_out/{name}", name=passing_samples(wildcards))

def bakta_outputs(wildcards):
    return expand("temp/bakta_out/{name}", name=passing_samples(wildcards))

def antismash_outputs(wildcards):
    return expand("temp/antismash_out/{name}", name=passing_samples(wildcards))

def antismash_csv_outputs(wildcards):
    return expand("temp/antismash_csv/{name}_overview.csv", name=passing_samples(wildcards))


rule all:
    input:
        expand("temp/filtered_reads/{name}.filtlong.fastq.gz", name=input_np_raw),
        expand("temp/qc/{name}.qc.tsv", name=input_np_raw),
        "results/qc_summary.tsv",
        flye_outputs,
        checkm2_outputs,
        "temp/gtdbtk_out",
        prokka_outputs,
        bakta_outputs,
        antismash_outputs,
        antismash_csv_outputs,
        "temp/getphylo_out/trees/combined_alignment.tree",
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


rule flye:
    input:
        NPreads="temp/filtered_reads/{sample}.filtlong.fastq.gz"
    output:
        asm="temp/all_assemblies/{sample}.flye.fa.gz",
        asminfo="temp/all_assemblies/{sample}.assembly_info.txt"
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
        
rule gtdb:
    input:
        flye_outputs,
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
    shell:
        """
        export GTDBTK_DATA_PATH="{params.db}"
        gtdbtk classify_wf \
        --genome_dir "temp/all_assemblies/" \
        --out_dir {output} \
        --cpus {threads} \
        --extension fa.gz
        """

rule prokka:
    input:
        asm="temp/all_assemblies/{sample}.flye.fa.gz",
    output:
        directory("temp/prokka_out/{sample}")
    threads: 8
    resources:
        mem_mb=lambda wc, input: max(5 * input.size_mb, 10240),
        node_type="general",
        time="00-05:00:00",
    conda:
        "envs/env_prokka.yml"
    shell:
        """
        # Decompress to a regular file in the output directory's parent folder
        gzip -dc {input.asm} > "temp/prokka_out/{wildcards.sample}_temp.fa"

        # Run Prokka on the temporary file
        prokka "temp/prokka_out/{wildcards.sample}_temp.fa" \
            --kingdom Bacteria \
            --outdir {output} \
            --force \
            --prefix {wildcards.sample} \
            --cpus {threads}

        # Clean up the temporary unzipped file
        rm "temp/prokka_out/{wildcards.sample}_temp.fa"
        """
        
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
    conda:
        "envs/env_antismash.yaml"
    shell:
        """
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
        
rule getphylo:
    input:
        prokka_outputs
    output:
        "temp/getphylo_out/trees/combined_alignment.tree"
    threads: 10
    resources:
        mem_mb=8192,
        node_type="general",
        time="00-05:00:00",
    conda:
        "envs/env_getphylo.yaml"
    shell:
        """
        rm -rf temp/getphylo_out

        WRAPPER=$(mktemp)
        printf '#!/bin/bash\\nif [ $# -eq 0 ]; then exit 0; fi\\nexec fasttree "$@"\\n' > "$WRAPPER"
        chmod +x "$WRAPPER"

        mkdir -p temp/getphylo_in
        for d in {input}; do
            sample=$(basename "$d")
            cp "$d/$sample.gbk" temp/getphylo_in/
        done
        getphylo -g "temp/getphylo_in/*.gbk" -o temp/getphylo_out -c {threads} --fasttree "$WRAPPER"
        rm -rf temp/getphylo_in "$WRAPPER"
        """

rule compile_results:
    input:
        checkm2_reports   = checkm2_outputs,
        flye_info         = flye_info_outputs,
        gtdbtk_summary    = "temp/gtdbtk_out",
        bakta_dirs        = bakta_outputs,
        antismash_overviews = antismash_csv_outputs,
        qc_summary        = "results/qc_summary.tsv"
    output:
        tsv = "results/aggregated_results.tsv"
    resources:
        mem_mb=8192,
        node_type="general",
        time="00-01:00:00",
    conda:
        "envs/R-main.yaml"
    script:
        "scripts/compile_results.R"
        
rule plot_tree:
    input:
        tree     = "temp/getphylo_out/trees/combined_alignment.tree",
        metadata = "results/aggregated_results.tsv"
    output:
        plot = "results/phylogenetic_tree.pdf"
    resources:
        mem_mb=4096,
        node_type="general",
        time="00-01:00:00",
    conda:
        "envs/R-main.yaml"
    script:
        "scripts/plot_tree.R"