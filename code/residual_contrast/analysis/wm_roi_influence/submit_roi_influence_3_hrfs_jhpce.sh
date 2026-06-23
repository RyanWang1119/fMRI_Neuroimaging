#!/bin/bash
#
# Submit the same WM ROI influence contrast for cHRF, cHRFderiv, and sHRF.
#
# Run from JHPCE:
#   cd /users/rwang/fMRI_Neuroimaging
#   bash code/residual_contrast/analysis/wm_roi_influence/submit_roi_influence_3_hrfs_jhpce.sh
#
# Optional overrides:
#   CONTRAST_NAME=Body_LoadDiff_vs_Place_LoadDiff WINDOW_SEC="4.32 8.64" bash ...

set -euo pipefail

REPO_DIR="${REPO_DIR:-/users/rwang/fMRI_Neuroimaging}"
DATA_DIR="${DATA_DIR:-/users/rwang}"
CONTRAST_NAME="${CONTRAST_NAME:-Body_LoadDiff_vs_Face_LoadDiff}"
WINDOW_SEC="${WINDOW_SEC:-4.32 8.64}"
OUTDIR="${OUTDIR:-${REPO_DIR}/data/task_residual/roi_influence_elasticnet}"
NREPEATS="${NREPEATS:-200}"
NPERM="${NPERM:-5000}"
USE_PARFOR="${USE_PARFOR:-1}"
MAKE_PLOTS="${MAKE_PLOTS:-0}"

SBATCH_FILE="${REPO_DIR}/code/residual_contrast/analysis/wm_roi_influence/run_roi_influence_elasticnet_jhpce.sbatch"

if [[ ! -f "${SBATCH_FILE}" ]]; then
    echo "Missing sbatch file: ${SBATCH_FILE}" >&2
    exit 1
fi

for f in WMcHRF.mat WMcHRFderiv.mat WMsHRF.mat; do
    path="${DATA_DIR}/${f}"
    if [[ ! -f "${path}" ]]; then
        echo "Missing data file: ${path}" >&2
        exit 1
    fi
    bytes=$(stat -c%s "${path}")
    if [[ "${bytes}" -lt 1048576 ]]; then
        echo "Data file is too small and may be a Git LFS pointer: ${path} (${bytes} bytes)" >&2
        exit 1
    fi
done

cd "${REPO_DIR}"

echo "Submitting ROI influence jobs:"
echo "  REPO_DIR=${REPO_DIR}"
echo "  DATA_DIR=${DATA_DIR}"
echo "  CONTRAST_NAME=${CONTRAST_NAME}"
echo "  WINDOW_SEC=${WINDOW_SEC}"
echo "  OUTDIR=${OUTDIR}"
echo "  NREPEATS=${NREPEATS}"
echo "  NPERM=${NPERM}"
echo "  USE_PARFOR=${USE_PARFOR}"
echo "  MAKE_PLOTS=${MAKE_PLOTS}"

for model in cHRF cHRFderiv sHRF; do
    echo "Submitting ${model}..."
    sbatch --export=ALL,REPO_DIR="${REPO_DIR}",DATA_DIR="${DATA_DIR}",HRF_MODEL="${model}",CONTRAST_NAME="${CONTRAST_NAME}",WINDOW_SEC="${WINDOW_SEC}",OUTDIR="${OUTDIR}",NREPEATS="${NREPEATS}",NPERM="${NPERM}",USE_PARFOR="${USE_PARFOR}",MAKE_PLOTS="${MAKE_PLOTS}" "${SBATCH_FILE}"
done
