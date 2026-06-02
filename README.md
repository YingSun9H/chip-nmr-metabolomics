# Map local approved-workspace fields to generic analysis names.
paths:
  phenotype_table: data/phenotypes.parquet
  metabolomics_table: data/nmr_metabolomics.parquet
  chip_table: data/chip_calls.parquet
  outcome_table: data/outcomes.parquet

columns:
  participant_id: participant_id
  age: age_at_baseline
  sex: sex
  body_mass_index: body_mass_index
  deprivation_index: deprivation_index
  smoking_status: smoking_status
  alcohol_status: alcohol_status
  genetic_principal_components: [genetic_pc_1, genetic_pc_2, genetic_pc_3, genetic_pc_4, genetic_pc_5, genetic_pc_6, genetic_pc_7, genetic_pc_8, genetic_pc_9, genetic_pc_10]
  assessment_date_baseline: nmr_baseline_date
  assessment_date_repeat: nmr_repeat_date
  chip_any: chip_any
  chip_large_clone: chip_large_clone
  chip_driver_group: chip_driver_group
  maximum_vaf: maximum_variant_allele_fraction
  lipid_lowering_medication: lipid_lowering_medication
  antihypertensive_medication: antihypertensive_medication
  kidney_function: estimated_glomerular_filtration_rate

outcomes:
  - cardiometabolic_disease
  - diabetes
  - coronary_heart_disease
  - stroke
  - heart_failure
  - atrial_fibrillation
  - cardiovascular_mortality
