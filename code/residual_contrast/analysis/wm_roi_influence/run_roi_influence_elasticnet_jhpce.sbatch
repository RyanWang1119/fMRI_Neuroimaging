#!/bin/bash
#
# JHPCE submission template for the WM ROI influence elastic-net module.
#
# Example:
#   sbatch --export=ALL,REPO_DIR=/users/rwang/fMRI_Neuroimaging,DATA_DIR=/users/rwang,HRF_MODEL=cHRF,CONTRAST_NAME=Body_LoadDiff_vs_Face_LoadDiff,WINDOW_SEC="4.32 8.64" \
#     code/residual_contrast/analysis/wm_roi_influence/run_roi_influence_elasticnet_jhpce.sbatch
#
# The script writes outputs under:
#   data/task_residual/roi_influence_elasticnet/<contrast>/<HRF_MODEL>/
#
# Reuse the same output root for cHRF, cHRFderiv, and sHRF so
# shared_subject_splits.mat is reused exactly across models.
# Heavy cluster jobs default to MAKE_PLOTS=0. Generate plots locally from
# the saved roi_influence_*.mat with run_roi_influence_local_postprocess.m.

#SBATCH --job-name=wm_roi_enet
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=wm_roi_enet_%x_%j.out
#SBATCH --error=wm_roi_enet_%x_%j.err

set -euo pipefail

REPO_DIR="${REPO_DIR:-/users/rwang/fMRI_Neuroimaging}"
DATA_DIR="${DATA_DIR:-/users/rwang}"
HRF_MODEL="${HRF_MODEL:-cHRF}"
CONTRAST_NAME="${CONTRAST_NAME:-Body_LoadDiff_vs_Face_LoadDiff}"
WINDOW_SEC="${WINDOW_SEC:-4.32 8.64}"
OUTDIR="${OUTDIR:-${REPO_DIR}/data/task_residual/roi_influence_elasticnet}"
NREPEATS="${NREPEATS:-200}"
NPERM="${NPERM:-5000}"
USE_PARFOR="${USE_PARFOR:-1}"
MAKE_PLOTS="${MAKE_PLOTS:-0}"

cd "${REPO_DIR}"

# Adjust this module line to match the MATLAB module available on JHPCE.
# If MATLAB is already on PATH for your account, this block is harmless.
if command -v module >/dev/null 2>&1; then
    module load matlab || module load matlab/R2025a || true
fi

export HRF_MODEL
export DATA_DIR
export CONTRAST_NAME
export WINDOW_SEC
export OUTDIR
export NREPEATS
export NPERM
export USE_PARFOR
export MAKE_PLOTS

echo "Repository: ${REPO_DIR}"
echo "DATA_DIR: ${DATA_DIR}"
echo "HRF_MODEL: ${HRF_MODEL}"
echo "CONTRAST_NAME: ${CONTRAST_NAME}"
echo "WINDOW_SEC: ${WINDOW_SEC}"
echo "OUTDIR: ${OUTDIR}"
echo "NREPEATS: ${NREPEATS}"
echo "NPERM: ${NPERM}"
echo "USE_PARFOR: ${USE_PARFOR}"
echo "MAKE_PLOTS: ${MAKE_PLOTS}"
echo "SLURM_CPUS_PER_TASK: ${SLURM_CPUS_PER_TASK:-unset}"

matlab -nodisplay -nosplash -r "addpath(genpath('code/residual_contrast/analysis/wm_roi_influence')); run_roi_influence_elasticnet_slurm; exit"
