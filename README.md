# CHIP-associated metabolic remodelling and cardiovascular outcomes

Reproducible analysis code and field definitions for the accompanying European Heart Journal manuscript.

## Study and analysis scope

The primary cohort comprises 451,743 UK Biobank participants; 17,842 have repeat NMR measurements. The panel contains 251 NMR parameters. Seven domain modules comprise 52 parameters, listed with UK Biobank field identifiers and parameter types in `metadata/module_membership.csv`.

The cardiovascular endpoints are coronary heart disease, stroke, heart failure, atrial fibrillation and cardiovascular-related mortality, in that order. Cardiovascular-related mortality includes a primary or secondary cardiovascular cause in linked death records. Prevalent diabetes appears only as a covariate or restriction condition.

Individual-level UK Biobank data are available only to approved researchers and are not distributed here. CHIP ascertainment and linked-record endpoint coding must be supplied as curated inputs under the applicable UK Biobank approvals. This code package implements cohort eligibility checks, metabolomic association analyses, metabolic modules, repeat-measurement analyses, survival analyses and sensitivity checks; it does not implement exome variant calling.

## Environment

Python 3.12 with the exact package versions in `requirements.txt`. The computational verification environment used Python 3.12.14, NumPy 2.3.5, pandas 3.0.1 and SciPy 1.16.2.

```sh
python -m venv .venv
```

Activate the environment before installing packages or running the commands below:

```powershell
# Windows PowerShell
.\.venv\Scripts\Activate.ps1
```

```sh
# Linux/macOS
. .venv/bin/activate
```

```sh
python -m pip install -r requirements.txt
```

## Participant data schema

`metadata/data_schema.csv` defines the expected researcher-supplied cohort export. Use a unique `participant_id`, primary covariates, curated CHIP indicators and the 251 raw NMR parameters at each available assessment. Biomarker column names are `baseline__<metabolite_id>` and `repeat__<metabolite_id>`. Endpoint status and time columns follow `<outcome>__<assessment>_status` and `<outcome>__<assessment>_time`.

Endpoint status is 0 for no incident event, 1 for an incident event and 99 for the corresponding prevalent condition. Follow-up times are years from the specified assessment. Nonfatal endpoints exclude the corresponding prevalent condition; mortality includes the full eligible cohort. Follow-up ends at the endpoint, death, loss to follow-up or administrative censoring. The manuscript uses administrative censoring on 8 July 2024.

Driver flags are curated presence indicators and may overlap across gene groups. Each driver contrast includes the corresponding carriers and participants without CHIP, excluding other CHIP carriers. DDR comprises PPM1D/TP53; spliceosome comprises SRSF2/SF3B1/ZRSR2. Large CHIP is maximum VAF ≥10%; small CHIP is CHIP with maximum VAF <10%.

## Run analyses

```sh
python analyse.py all --cohort private_data/cohort.csv --output outputs/analysis
```

The `baseline`, `modules`, `repeat`, `outcomes` and `sensitivity` stages can be run separately. The repeat and outcome stages consume `baseline_biomarker_associations.csv` from the baseline stage in the same output directory.

- Baseline/repeat cross-sectional biomarker models: winsorisation at 0.5/99.5 percentiles, rank inverse-normal transformation, ordinary least squares and HC1 robust standard errors. Transformation occurs within each complete-case biomarker/exposure sample.
- Module models: average available transformed constituents, requiring at least half (rounded up). Ingredients are not sign-reversed or risk-weighted. PCA uses mean imputation only among score-eligible participants, with the first component oriented to the mean score.
- Signature: baseline any-CHIP coefficients for parameters with FDR <0.05; at least half of the selected parameters are required. Baseline survival analyses standardise biomarker transforms before weighting. Repeat analyses retain baseline rank-normal reference distributions and weights, and standardise weighted sums using their baseline cohort mean and SD. Change is follow-up minus baseline on that baseline score scale. Module changes are standardised using their change distribution.
- Paired biomarker change (Figure 3B): baseline-clipped concentration change on the baseline SD scale, with joint small- and large-clone indicators and no-CHIP reference; primary covariates and HC1 covariance. Baseline biomarker value and NMR interval are not covariates in this model. FDR covers 251 large-clone coefficients.
- Paired module change (Supplementary Figure S3B): baseline-anchored rank-normal module scores; change is standardised in the paired score sample. Separate any-CHIP and large-clone contrasts use no-CHIP controls, with small-clone carriers excluded from the latter. Models include baseline module and NMR interval, with residual-variance OLS covariance and 1.96-SE intervals.
- Signature trajectory (Supplementary Figure S3A): raw paired means with 1.96-SE intervals. Adjusted change models include baseline score, NMR interval and primary covariates, using residual-variance OLS covariance. The large-clone indicator contrast retains small-clone participants in its reference group.
- Survival models: Breslow partial likelihood, numerical optimisation and observed-information standard errors. Fits must converge. Splines use baseline score percentiles 5, 35, 65 and 95, with score zero as the HR reference.
- Diagnostics: scaled Schoenfeld residual tests against log follow-up time, including covariates and global tests. Five-year period analyses and a two-year landmark sensitivity analysis are included.
- Sensitivity models: clinical adjustment, disease-free restriction and six subgroup interaction dimensions. Prevalent diabetes is a clinical covariate, not a clinical endpoint.

The score is derived and evaluated in the same cohort. These analyses estimate associations and do not establish mediation, causality or incremental prediction utility.

## Multiplicity families

| Analysis | Benjamini–Hochberg family |
|---|---|
| Biomarker associations | 251 parameters per exposure and assessment/model |
| Modules | Seven modules per exposure and stratum |
| Baseline CHIP–outcome models | Ten tests: two exposures × five endpoints |
| Two-year lag models | Ten tests |
| Five-year period models | Ten tests separately within each period |
| Follow-up signature and score change | Five endpoints separately for each predictor |
| Module-change outcomes | 35 tests: seven modules × five endpoints |
| Paired biomarker changes | 251 large-clone coefficients from joint clone-size models |
| Module interactions | 42 tests: seven modules × six interaction dimensions |

## Exact contrast and coding definitions

Biomarker MWAS and cardiovascular large-clone models compare the large-clone indicator with all participants without a large clone, including small-clone CHIP. Module clone-size analyses instead compare each clone-size group with participants without CHIP; small-clone carriers are excluded from large-clone module contrasts. Driver-specific contrasts use presence indicators and participants without CHIP as controls. Multiple driver indicators may be positive in one participant.

Smoking is coded 0/1/2 with 0 as reference; alcohol is coded 1/2/3 with 1 as reference. Unexpected codes cause an error and missing categories are excluded from complete-case models. `metadata/categorical_coding.csv` defines the exact design. Prevalent CVD is prevalent CHD, stroke or HF before baseline NMR; AF alone does not qualify.

The modules stage exports all 84 primary any-CHIP stratum-specific models (six characteristics, two strata and seven modules) and 42 pooled any-CHIP interactions. Strata are sex, age (<60/≥60), BMI (≤25/>25), prevalent CVD, prevalent diabetes and lipid-lowering medication. Transforms are estimated within the exposure/stratum complete-case sample; pooled interaction models use their own complete sample. Stratum-specific coefficients are descriptive; the pooled interaction tests assess heterogeneity.

The primary Cox engine uses Breslow ties, L-BFGS-B optimisation with `ftol=1e-8`, `gtol=1e-5`, at most 500 iterations and observed-information standard errors. Continuous covariates are standardised, with design columns ordered as in the manuscript analyses. Independent checks use statsmodels PHReg with Breslow ties. The optional `ftol=1e-13, refine=True` setting permits a more precise numerical comparison, without changing the documented primary setting.

```sh
python validate_numerics.py --cohort private_data/cohort.csv --output outputs/numerical_validation.json
```

This independent check requires the optional statsmodels dependency already pinned in `requirements.txt`. No private cohort path or participant records are included in the distribution.


## Distribution scope

This directory contains analysis code and field definitions only. No participant data, manuscript result tables or figure assets are included. `metadata/outcome_definitions.csv` contains endpoint definitions only, without observed sample sizes, event counts or follow-up summaries. The separate submission analysis package contains publication aggregate results and plotting scripts.

Subgroup analyses use only any CHIP and the primary covariate model. No large-clone, repeat-assessment or clinical-adjustment subgroup analyses are included.

`metadata/exposure_contrasts.csv` specifies the reference group separately for each analysis.
