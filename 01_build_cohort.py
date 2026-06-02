# CHIP-NMR metabolomics analysis code

This repository contains analysis templates for reproducing the reported CHIP-NMR metabolomics analyses within an approved UK Biobank workspace. It intentionally does not include individual-level UK Biobank data, derived individual-level data, or project-specific local column names.

The scripts are designed as a transparent analysis framework rather than a standalone public dataset. Users must supply approved UK Biobank extracts and map local fields to the generic schema in `config/column_map_template.yaml`.

## What is included

- `config/column_map_template.yaml`: maps generic analysis variable names to fields available inside an approved workspace.
- `src/01_build_cohort.py`: cohort assembly template.
- `src/02_preprocess_metabolomics.py`: winsorization, rank-inverse-normal transformation, standardization and module-score construction.
- `src/03_run_metabolome_wide_models.py`: metabolome-wide CHIP association models.
- `src/04_run_module_score_models.py`: module, clone-size, driver-gene and subgroup models.
- `src/05_run_outcome_models.py`: Cox outcome and landmark models.
- `src/06_make_source_data_and_figures.py`: figure/source-data assembly templates.
- `tables/module_definitions_template.csv`: editable module definition table.
- `tables/biomarker_derivation_groups_template.csv`: template for the Fig. 1A biomarker-parameter derivation groups.

## Analysis scope

The biomarker-wide NMR results should be interpreted as parameter-level associations rather than independent metabolite discoveries. The UK Biobank/Nightingale NMR panel includes absolute or concentration-based biomarkers, other derived or composite biomarkers, and direct composition-ratio biomarkers.

For the Fig. 1A restricted summaries, each NMR biomarker parameter is assigned to one of three reader-facing derivation groups:

- `Absolute/concentration-based`: absolute biomarker levels, lipid or metabolite concentrations, and lipoprotein particle concentrations.
- `Other derived/composite`: derived or composite readouts that are not direct composition ratios, such as particle-diameter summaries, degree of unsaturation, remnant cholesterol and glucose-lactate.
- `Direct composition-ratio`: percentage or ratio readouts, including to-total-lipids percentages, to-total-fatty-acids percentages, apolipoprotein ratio and fatty-acid ratios.

In the analytic 251-parameter panel used for the manuscript, this classification yielded 164 absolute/concentration-based parameters, 6 other derived/composite parameters and 81 direct composition-ratio parameters. The primary metabolome-wide analyses used all 251 parameters; these groups were used only for interpretation and restricted Fig. 1A summaries.

The seven metabolic modules are prespecified domain summaries based on selected representative NMR biomarkers. They are not intended to be an exhaustive or mutually exclusive classification of all measured NMR parameters. The manuscript reports a principal-component consistency check showing that the mean-based module scores agree closely with the corresponding oriented first principal component.

The CHIP-like metabolic score is provided as a biological and disease-relevance summary of CHIP-associated metabolism, not as a prediction-ready clinical model.

## Data requirements

Users must provide their own approved UK Biobank data extracts, including participant identifier, assessment dates, baseline covariates, CHIP calls, NMR metabolomic biomarkers, linked outcome fields and follow-up times. The scripts expect a configuration file that maps local field names to the generic schema used here.

Individual-level UK Biobank data cannot be redistributed by the authors.

## Reproducible workflow

1. Create an environment using `requirements.txt`.
2. Copy `config/column_map_template.yaml` to `config/column_map.yaml` and map local UK Biobank fields. If generating Fig. 1A source-data summaries, provide biomarker metadata with biomarker names and original measure-type annotations or adapt `tables/biomarker_derivation_groups_template.csv`.
3. Place approved workspace data in a local `data/` directory, which is excluded by `.gitignore`.
4. Run the scripts in numerical order.
5. Compare generated summary tables with the manuscript Supplementary Tables and Source Data files.

The public Source Data workbooks accompanying the manuscript contain derived summary statistics only. They do not contain individual-level UK Biobank data.

## Citation

Please cite the associated manuscript and the archived repository DOI once available.
