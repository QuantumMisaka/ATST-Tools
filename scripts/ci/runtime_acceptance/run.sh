#!/usr/bin/env bash
set -euo pipefail
task_root=$1
exec 3>&1
exec > >(tee "$task_root/runtime.log" >&2) 2>&1
finish() {
    result=$?
    trap - EXIT
    printf '%s\n' "$result" > "$task_root/exit-code.txt"
    set +e
    tar --exclude='venv' --exclude='abacus_basic_gpu' \
        --exclude='*-CHARGE-DENSITY.restart' --exclude='potes.cube' \
        --exclude='ElecStaticPot.cube' -czf - -C "$task_root" \
        exit-code.txt runtime.log environment.txt source-identity.json \
        pins.sh stages.tsv harness >&3
    archive_result=$?
    if [ "$result" -eq 0 ]; then result=$archive_result; fi
    exit "$result"
}
trap finish EXIT
source "$task_root/pins.sh"
export PYTHONDONTWRITEBYTECODE=1
export XDG_CACHE_HOME="$task_root/cache"
export MPLCONFIGDIR="$task_root/cache/matplotlib"
export PIP_CACHE_DIR="$task_root/cache/pip"
export NUMBA_CACHE_DIR="$task_root/cache/numba"
export TMPDIR="$task_root/tmp"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
mkdir -p "$task_root/cache" "$task_root/tmp"
touch "$task_root/stages.tsv"
set +u
source /opt/modules/lmod/9.2.4/init/bash
module use /opt/modules/modulefiles/apps /opt/modules/modulefiles/devtools
module load abacus/v3.11.0-beta9-auto
{
    module list
    command -v mpirun
    mpirun --version
    nvidia-smi --query-gpu=name,driver_version --format=csv
} > "$task_root/environment.txt" 2>&1
set -u
cat "$task_root/environment.txt"
sha256sum "$task_root/abacus_basic_gpu" "$task_root"/wheels/*.whl
export ABACUS_COMMAND="$task_root/abacus_basic_gpu"
export ABACUS_VERSION_COMMAND="$task_root/abacus_basic_gpu --version"
export ATST_WHEEL="$task_root/wheels/$ATST_WHEEL_NAME"
export PLUGIN_WHEEL="$task_root/wheels/$PLUGIN_WHEEL_NAME"
site_python=/home/liuzhaoqing/.conda/envs/atst-dev/bin/python
base_path=$PATH
failures=0
stage() {
    name=$1
    shift
    printf 'BEGIN stage=%s\n' "$name"
    if "$@"; then status=0; else status=$?; failures=1; fi
    printf '%s\t%s\n' "$name" "$status" >> "$task_root/stages.tsv"
    printf 'END stage=%s status=%s\n' "$name" "$status"
    return "$status"
}
if ! stage venv "$site_python" -m venv --system-site-packages "$task_root/venv"; then exit "$failures"; fi
export PATH="$task_root/venv/bin:$base_path"
if ! stage install python -m pip install --no-index --no-deps \
    "$ATST_WHEEL" "$PLUGIN_WHEEL" "$task_root/wheels/$SELLA_WHEEL_NAME"; then exit "$failures"; fi
if ! stage pip_check python -m pip check; then exit "$failures"; fi
python -m pip freeze > "$task_root/harness/pip-freeze.txt"
if ! stage mpi_probe python -c 'from mpi4py import MPI; print(MPI.Get_library_version())'; then exit "$failures"; fi
cd "$task_root/harness"
if ! stage provenance python -u acceptance.py provenance; then exit "$failures"; fi
for case in relax md_initial md_restart md_native neb_serial cp_baseline cp vibration; do
    stage "$case" timeout 40m python -u acceptance.py run-case "$case" || true
done
stage neb_parallel timeout 30m mpirun --oversubscribe --bind-to none -np 2 \
    python -u acceptance.py run-case neb_parallel || true
for case in sella ccqn irc; do
    stage "$case" timeout 40m python -u acceptance.py run-case "$case" || true
done
export PATH="$base_path"
set +u
module purge
module load abacus/LTSv3.10.1-sm70-auto
set -u
export PATH="$task_root/venv/bin:$PATH"
export ABACUS_LTS_COMMAND="$(command -v abacus)"
"$ABACUS_LTS_COMMAND" --version
sha256sum "$ABACUS_LTS_COMMAND"
stage lts_legacy timeout 20m python -u acceptance.py run-case lts_legacy || true
exit "$failures"
