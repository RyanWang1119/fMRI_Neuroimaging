This repository contains the analysis pipeline and code used for the thesis investigating systematic, task-dependent information within fMRI residuals. Using datasets from the Human Connectome Project (HCP), this project tests whether residuals contain structured information that can reliably predict task states.

## WM ROI Influence Elastic-Net Module

The separate module in `code/residual_contrast/analysis/wm_roi_influence` identifies stable, high-confidence influential candidate parcels for WM residual interaction contrasts without modifying the existing linear-SVM pipeline.

The driver `run_roi_influence_elasticnet.m` uses subject-level paired differences:

`D = (2-back category A - 0-back category A) - (2-back category B - 0-back category B)`

and fits nested elastic-net logistic models with grouped subject folds using the symmetric paired representation `[D/2; -D/2]`. Main outputs are kept distinct:

- Haufe pattern: task-linked residual signal map.
- Elastic-net selection frequency: stability of sparse multivariable selection.
- Held-out model reliance: dependence of the fitted decoder on each parcel under held-out paired half-swaps.

The population-level sign-flip maxT regional test is separate from the decoder and uses only unstandardized subject-level paired differences. The `high_confidence_influential_candidate` flag is a reproducibility-based candidate label, not a causal, mechanistic, exclusive, or definitive biological claim.

Example MATLAB command:

```matlab
addpath(genpath('code/residual_contrast/analysis/wm_roi_influence'));
cfg = struct();
cfg.hrfModelName = 'cHRF';
cfg.contrastName = 'Body_LoadDiff_vs_Face_LoadDiff';
cfg.windowSec = [4.32 8.64];
run_roi_influence_elasticnet(cfg);
```
