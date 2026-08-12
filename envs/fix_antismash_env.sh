#!/bin/bash
#
# Make the antiSMASH conda environment portable across CPU generations.
#
# THE PROBLEM
#   antiSMASH pulls in MOODS-python, whose setup.py hardcodes:
#       common_compile_args = ['-march=native', '-O3', '-fPIC', '--std=c++11']
#   -march=native compiles for the CPU that happens to run the install, with no
#   runtime dispatch and no fallback path. Build the env on an AVX-512 machine
#   and antiSMASH dies with SIGILL ("Illegal instruction") the instant it
#   imports MOODS on any node without AVX-512. On 2026-08-05 that killed 526 of
#   549 genomes; the ~4% that survived had landed on zen5 nodes.
#
#   Setting CFLAGS/CXXFLAGS does NOT help: setuptools appends extra_compile_args
#   AFTER the environment flags, and for gcc the last -march wins. The only
#   reliable fixes are to patch setup.py or to append a flag after MOODS's own.
#   This script does the latter, via compiler wrappers.
#
# USAGE
#   envs/fix_antismash_env.sh --check     # is the current env portable?
#   envs/fix_antismash_env.sh --rebuild   # rebuild MOODS portably
#
#   PORTABLE_MARCH=x86-64-v3 envs/fix_antismash_env.sh --rebuild
#
#   x86-64-v2 (default) = SSE4.2/POPCNT, safe on anything modern.
#   x86-64-v3           = AVX2, supported by zen3 and zen5 but not older nodes.
#
set -euo pipefail

MODE="${1:---check}"
MARCH="${PORTABLE_MARCH:-x86-64-v2}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

log() { printf '%s\n' "$*"; }
die() { printf 'ERROR: %s\n' "$*" >&2; exit 1; }

# ── locate the antismash conda env ───────────────────────────────────────────
# Found by content rather than by hash: the hash is derived from the YAML, so it
# changes whenever the environment definition changes.
find_env() {
    local d
    for d in .snakemake/conda/*/; do
        [ -x "${d}bin/antismash" ] && { printf '%s' "${d%/}"; return 0; }
    done
    return 1
}

# ── is any MOODS extension built with AVX-512? ───────────────────────────────
# Host-independent: inspects the binaries rather than trying to run them, so it
# gives the same answer on a login node as on a compute node.
scan_moods() {
    local env_dir="$1" so found=0
    if ! command -v objdump >/dev/null 2>&1; then
        log "  objdump not available; skipping instruction scan"
        return 2
    fi
    while IFS= read -r so; do
        if objdump -d "$so" 2>/dev/null | grep -qE '%zmm[0-9]'; then
            log "  AVX-512 found: ${so#"$env_dir"/}"
            found=1
        fi
    done < <(find "$env_dir" -path '*/MOODS/*.so' 2>/dev/null)
    return $found
}

ENV_DIR="$(find_env)" || die "no conda env with bin/antismash under .snakemake/conda/.
Run 'snakemake --use-conda --conda-create-envs-only --cores 1' first."
log "antismash env: $ENV_DIR"

# ── check mode ───────────────────────────────────────────────────────────────
if [ "$MODE" = "--check" ]; then
    log "Scanning MOODS extensions for AVX-512 instructions..."
    if scan_moods "$ENV_DIR"; then
        log "OK: no AVX-512 in MOODS. Safe to run on any node."
        exit 0
    else
        rc=$?
        [ "$rc" = "2" ] && exit 0
        log ""
        log "MOODS carries AVX-512 instructions and has no runtime fallback."
        log "antiSMASH will die with SIGILL on nodes without AVX-512."
        log "Fix with: envs/fix_antismash_env.sh --rebuild"
        exit 1
    fi
fi

[ "$MODE" = "--rebuild" ] || die "unknown mode '$MODE' (use --check or --rebuild)"

# ── rebuild mode ─────────────────────────────────────────────────────────────
PY="$ENV_DIR/bin/python"
PIP="$ENV_DIR/bin/pip"
[ -x "$PY" ]  || die "no python in $ENV_DIR"
[ -x "$PIP" ] || die "no pip in $ENV_DIR"

MOODS_VER="$("$PIP" show MOODS-python 2>/dev/null | awk '/^Version:/{print $2}')"
[ -n "$MOODS_VER" ] || die "MOODS-python is not installed in $ENV_DIR"
log "MOODS-python version: $MOODS_VER"

# The compiler that setuptools will actually invoke. Taking the name from the
# interpreter's own build config rather than assuming 'gcc' matters: conda
# pythons commonly report something like x86_64-conda-linux-gnu-cc, and a
# wrapper named 'gcc' would then never be consulted.
CC_NAME="$("$PY" -c 'import sysconfig,shlex; print(shlex.split(sysconfig.get_config_var("CC") or "gcc")[0])')"
CXX_NAME="$("$PY" -c 'import sysconfig,shlex; v=sysconfig.get_config_var("CXX") or "g++"; print(shlex.split(v)[0])')"
log "compiler drivers: CC=$CC_NAME CXX=$CXX_NAME"

WRAP_DIR="$(mktemp -d -t portable_cc_XXXXXX)"
trap 'rm -rf "$WRAP_DIR"' EXIT

# Verify the compiler understands the requested baseline before relying on it;
# -march=x86-64-v2 needs gcc >= 11. Fall back to the universal baseline.
REAL_CC="$(command -v "$CC_NAME" || command -v gcc)" || die "no C compiler found"
if ! echo 'int main(void){return 0;}' | "$REAL_CC" -march="$MARCH" -x c - -o /dev/null 2>/dev/null; then
    log "compiler does not support -march=$MARCH; falling back to x86-64"
    MARCH="x86-64"
fi
log "portable baseline: -march=$MARCH"

# The wrapper appends its flag LAST, which is the whole point: MOODS puts
# -march=native in extra_compile_args, setuptools appends those after anything
# from the environment, and for gcc the final -march wins.
make_wrapper() {
    local name="$1" real
    real="$(command -v "$name" 2>/dev/null || true)"
    [ -n "$real" ] || return 0
    cat > "$WRAP_DIR/$name" <<EOF
#!/bin/bash
exec "$real" "\$@" -march=$MARCH -mtune=generic
EOF
    chmod +x "$WRAP_DIR/$name"
    log "  wrapped $name -> $real"
}
for n in "$CC_NAME" "$CXX_NAME" gcc g++ cc c++; do make_wrapper "$n"; done

log "Rebuilding MOODS-python==$MOODS_VER from source..."
PATH="$WRAP_DIR:$PATH" "$PIP" install \
    --force-reinstall --no-binary :all: --no-cache-dir --no-deps \
    "MOODS-python==$MOODS_VER"

# ── verify ───────────────────────────────────────────────────────────────────
log ""
log "Verifying..."
if scan_moods "$ENV_DIR"; then
    log "  no AVX-512 in MOODS extensions"
else
    die "MOODS still contains AVX-512 after rebuild. The wrapper did not take
effect — check that '$CC_NAME' is the driver setuptools actually invoked."
fi

"$PY" -c 'import MOODS.scan, MOODS.tools, MOODS.parsers; print("  MOODS imports cleanly")' \
    || die "MOODS still fails to import"

ASVER="$("$ENV_DIR/bin/antismash" --version 2>&1 | head -1 || true)"
log "  antismash: $ASVER"
printf '%s\n' "$ASVER" > envs/antismash_version.lock
log ""
log "Done. Recorded version in envs/antismash_version.lock"
log "Now clear the stale failure markers and re-run:"
log "    find temp/antismash_out -name ANTISMASH_FAILED -delete"
