# amet TODO

## Reconcile `i_total_resid` vs `i_norm`

amet currently has two distinct ways of expressing a "methylation-decoupled" within-cell entropy, and the workflow uses both under different names:

- `i_norm` (in `workflow/scripts/eval_*.R` and `simulations_report.Rmd`): analytical normalization, defined as `i_total / (k_max * H(p_hat))`. Headline score used in the simulations and tool-comparison benchmarks.
- `i_total_resid` (in `workflow/Rmd/crc_windows_sce.Rmd`): empirical residuals from a per-window `lm(i_total ~ mean_meth + I(mean_meth^2))` fit. Used as the input to the SCE-based differential entropy testing and the per-cell embeddings.

These are different quantities computed by different math. Pick one canonical decoupling strategy (or document the regimes where each is preferred) and harmonize naming across the simulations, evals, dataset Rmds, and figure Rmds.

