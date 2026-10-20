#!/bin/bash
#SBATCH --job-name=abacus-latest-build
#SBATCH --partition=4V100
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --time=04:00:00
#SBATCH --qos=rush-1o2gpu
set -euo pipefail

ROOT="$HOME/atst-cp-parity-20261019"
SRC="$ROOT/abacus-latest-src"
BUILD="$ROOT/abacus-latest-build"
PREFIX="$ROOT/abacus-latest-install"
CEREAL_SRC="/opt/apps/abacus/abacus-develop-LTSv3.10.1/toolchain/install/cereal-master"
CEREAL_BUILD="$ROOT/cereal-build"
CEREAL_PREFIX="$ROOT/cereal-install"
CASE="$ROOT/case"
LOGS="$ROOT/logs"
mkdir -p "$LOGS"

set +u
source /opt/modules/lmod/9.2.4/init/bash
module use /opt/modules/modulefiles/devtools /opt/modules/modulefiles/apps
module load cmake/3.31.6
module load cuda/12.9.1
module load openmpi/5.0.10-nvhpc26.3-gnu-cuda12-auto
module load elpa/2026.02.001-2603-gnu
module load fftw/3.3.10
set -u
module list 2>&1 | sed -n '1,15p'
export OMP_NUM_THREADS=8
export CUDA_VISIBLE_DEVICES=0

echo "job=${SLURM_JOB_ID:-local} host=$(hostname)"
echo "cereal=$CEREAL_SRC"
echo "elpa=${ELPA_ROOT:-unset} scalapack=${SCALAPACK_ROOT:-unset}"

rm -rf "$CEREAL_BUILD"
cmake -S "$CEREAL_SRC" -B "$CEREAL_BUILD" \
  -DCMAKE_BUILD_TYPE=Release \
  -DJUST_INSTALL_CEREAL=ON \
  -DBUILD_TESTS=OFF \
  -DBUILD_DOC=OFF \
  -DBUILD_SANDBOX=OFF \
  -DWITH_WERROR=OFF \
  -DCMAKE_INSTALL_PREFIX="$CEREAL_PREFIX"
cmake --build "$CEREAL_BUILD" --target install -j8

cmake -S "$SRC" -B "$BUILD" \
  -DCMAKE_BUILD_TYPE=Release \
  -DENABLE_LCAO=ON \
  -DENABLE_ELPA=ON \
  -DENABLE_LIBXC=OFF \
  -DUSE_CUDA=ON \
  -DCMAKE_CUDA_ARCHITECTURES=70 \
  -DCMAKE_C_COMPILER=mpicc \
  -DCMAKE_CXX_COMPILER=mpicxx \
  -DCMAKE_PREFIX_PATH="$CEREAL_PREFIX:$ELPA_ROOT" \
  -DCEREAL_INCLUDE_DIR="$CEREAL_PREFIX/include" \
  -DCMAKE_INSTALL_PREFIX="$PREFIX"
cmake --build "$BUILD" -j8
cmake --install "$BUILD"

echo "=== built version ==="
"$PREFIX/bin/abacus" --version 2>&1 | head -3 || true

export PATH="$PREFIX/bin:$PATH"
export PYTHONPATH="$ROOT/atst-tools/src"
PY="$ROOT/venv-atst/bin/python"
REF_STATUS=0
(cd "$CASE" && "$PY" -m atst_tools.scripts.cli run config_reference_latest.yaml) >"$LOGS/latest_reference.log" 2>&1 || REF_STATUS=$?
COMP_STATUS=0
(cd "$CASE" && "$PY" -m atst_tools.scripts.cli run config_compensated_latest.yaml) >"$LOGS/latest_compensated.log" 2>&1 || COMP_STATUS=$?

"$PY" - <<PY
import json
from pathlib import Path
root = Path("$CASE")
for label, name in (("reference", "constant_potential_results_latest.json"), ("compensated", "constant_potential_results_compensated_latest.json")):
    path = root / name
    if not path.is_file():
        print(label, "missing", path); continue
    data = json.loads(path.read_text()); points = data.get("points") or []
    print(label, "status", data.get("status"), "points", len(points))
    for point in points:
        print(" ", {k: point.get(k) for k in ("nelec","efermi","mu_calc","residual_mu","energy","free_energy","cp_iterations","cp_converged")})
PY
echo "reference_status=$REF_STATUS compensated_status=$COMP_STATUS"
[ "$REF_STATUS" -eq 0 ] && [ "$COMP_STATUS" -eq 0 ]
