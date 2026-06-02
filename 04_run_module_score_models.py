"""Build the analytic cohort from approved-workspace extracts.

This template uses generic variable names. Users should map local fields in
config/column_map.yaml and keep individual-level data outside the repository.
"""

from __future__ import annotations

import argparse
import yaml
import pandas as pd


def main(config_path: str, output_path: str) -> None:
    with open(config_path, "r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    paths = config["paths"]
    phenotype = pd.read_parquet(paths["phenotype_table"])
    metabolomics = pd.read_parquet(paths["metabolomics_table"])
    chip = pd.read_parquet(paths["chip_table"])
    outcomes = pd.read_parquet(paths["outcome_table"])
    participant_id = config["columns"]["participant_id"]
    cohort = phenotype.merge(metabolomics, on=participant_id, how="inner")
    cohort = cohort.merge(chip, on=participant_id, how="inner")
    cohort = cohort.merge(outcomes, on=participant_id, how="left")
    cohort.to_parquet(output_path, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="config/column_map.yaml")
    parser.add_argument("--out", default="derived_data/analysis_cohort.parquet")
    args = parser.parse_args()
    main(args.config, args.out)
