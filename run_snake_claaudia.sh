#!/bin/bash
#SBATCH --job-name=tiny_earth_dk
#SBATCH --output=slurm_logs/%x-%j.out
#SBATCH --error=slurm_logs/%x-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1   # Adjust this to the desired number of threads
#SBATCH --mem=16G           # Adjust this to the desired memory allocation
#SBATCH --time=07-00:00:00     # Adjust this to the desired time limit
#SBATCH --mail-user=rhk@bio.aau.dk   # Email address for notifications
#SBATCH --mail-type=BEGIN,END,FAIL    # When to send email notifications

# Exit on first error and if any variables are unset
set -eu

eval "$(conda shell.bash hook)"
conda activate /home/bio.aau.dk/ur36rv/.conda/envs/tiny_earth_denmark

snakemake --use-conda --conda-frontend conda --conda-create-envs-only --cores 1 --latency-wait 30

# Guard against the MOODS/-march=native trap before committing anything to the
# cluster. MOODS compiles for the CPU that builds it, with no runtime fallback,
# so an env built on an AVX-512 machine makes antiSMASH die with SIGILL on every
# node that lacks AVX-512. On 2026-08-05 that wasted 526 of 549 jobs and was
# invisible until the results were inspected days later. The check inspects the
# binaries rather than running them, so it is valid from the submit node.
if ! bash envs/fix_antismash_env.sh --check; then
    echo ""
    echo "Refusing to submit: the antismash environment is not portable." >&2
    echo "Run 'envs/fix_antismash_env.sh --rebuild', then re-run this script." >&2
    exit 1
fi

snakemake --profile profile/

# cluster command not in snakemake 8... downgrading
#snakemake --cluster "sbatch --parsable --output=jobs/{rule}/slurm_%x_%j.out --error=jobs/{rule}/slurm_%x_%j.log --mem={resources.mem_mb} --cpus-per-task={threads} --time 00-10:00:00" -j 400 --use-conda
#snakemake --cluster qsub --jobs 32
#snakemake \
#         --executor cluster-generic \
#         --cluster-generic-submit-cmd 'qsub -N {rule} -q all.q -l h_vmem=8G -pe smp {threads} -V -cwd'
