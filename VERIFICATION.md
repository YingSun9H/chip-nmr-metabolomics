# Numerical verification

## Participant-level reproduction

The documented pipeline was run on the approved primary cohort of 451,743 participants. Checks covered 48 biomarker models, 63 module models, 84 primary any-CHIP stratum-specific module models and all 42 pooled module interactions (237 linear models in total). Sample sizes matched the manuscript results exactly. Maximum absolute coefficient differences were below 1.6e-11 and maximum P/FDR differences below 1.1e-9.

Primary-setting Breslow fits were checked for CHD with any CHIP, HF with large-clone CHIP, cardiovascular-related mortality with any CHIP, and CHD with the signature spline. Participant and event counts matched. Hazard-ratio differences were zero to displayed precision; the CHD spline HR and interval curves differed by less than 3.4e-13 across all 140 points. These checks cover selected survival models, not a new participant-level rerun of every survival analysis.

## Independent numerical checks

Seven primary any-CHIP module models and the sex-by-CHIP LDL-related module interaction were also compared with statsmodels OLS using HC1 covariance. The interaction coefficient and standard error agreed within 8e-13. The four selected survival models were compared with statsmodels PHReg using Breslow ties. With tighter optimisation (`ftol=1e-13`) and Newton refinement, maximum absolute coefficient and standard-error differences were below 8e-12 and 2e-13, respectively.

The primary manuscript setting remains `ftol=1e-8`, `gtol=1e-5` and a maximum of 500 L-BFGS-B iterations. It reproduces the supplied results. Across the four checked models, its maximum coefficient deviation from the stricter numerical solution was 0.157 standard-error units. The optional stricter setting is used for independent engine verification and is reported separately from the primary manuscript estimates.

## Environment and scope

Python 3.12.14, NumPy 2.3.5, pandas 3.0.1, SciPy 1.16.2 and statsmodels 0.14.5 were used. No participant-level records or absolute data paths are distributed. `validate_numerics.py` supports independent checks on a researcher-supplied cohort; The separate submission package includes `verify_aggregate.py` for its supplied result tables and FDR families.

## Paired-model verification and exposure references

An additional 26 adjusted paired models were checked at participant level: ten illustrative biomarker changes in Figure 3B, all 14 module changes in Supplementary Figure S3B and both adjusted signature-change estimates in Supplementary Figure S3A. All analytic sample sizes matched exactly. Figure 3B coefficient differences were below 1e-14; module coefficient differences were below 3.7e-11. Signature-change coefficients and CI limits differed by less than 1.1e-9 and P values by less than 8.1e-9. The full 251-biomarker large-clone FDR family was independently recomputed, with no significant change result. One joint clone-size biomarker model was also independently checked with statsmodels HC1 covariance.

These paired analyses use distinct models. Figure 3B uses clipped concentration change on the baseline SD scale with joint small/large indicators, no-CHIP reference and HC1 covariance. Module changes use no-CHIP controls, excluding small clones from the large-clone contrast, and residual-variance OLS covariance. Adjusted signature change uses indicator contrasts and residual-variance OLS covariance; the large-clone reference includes small clones. Exact references are listed in `metadata/exposure_contrasts.csv`. The supplied estimates were retained; code and reporting were aligned to their generating models.

## Version 1.1.0 release checks (4 October 2026)

The public package was checked in the pinned Python 3.12.14 environment using 1,600 synthetic records, all 251 biomarkers and all five cardiovascular endpoints. The `all` stage completed and generated 21 analysis CSV files. Checks covered output dimensions, endpoint sets, confidence intervals, multiplicity families, categorical coding and the distinct clone-size and driver reference groups. All seven module models and three selected Cox models passed `validate_numerics.py` comparisons with statsmodels. Separate HC1 and Breslow checks included tied event times, censoring, missing covariates and excluded prevalent records.

The corresponding submission package also passed its aggregate-result consistency checks. These release checks are additional software validation; they do not constitute a new full-cohort rerun. The public analysis scripts are identical to those in the checked submission package. Synthetic records, validation outputs, participant-level data and publication aggregate result tables are not included in this release. The endpoint metadata contains definitions only.
