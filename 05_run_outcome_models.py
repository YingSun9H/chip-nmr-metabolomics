"""Transform NMR biomarkers and construct module scores."""

from __future__ import annotations

import argparse
import pandas as pd
from utils import winsorized_rank_inverse_normal


def build_module_scores(df: pd.DataFrame, module_definitions: pd.DataFrame) -> pd.DataFrame:
    out = pd.DataFrame(index=df.index)
    for module_id, group in module_definitions.groupby("module_id"):
        transformed = []
        for biomarker in group["biomarker_id"]:
            transformed.append(winsorized_rank_inverse_normal(df[biomarker]).rename(biomarker))
        matrix = pd.concat(transformed, axis=1)
        min_required = max(1, int(len(group) * float(group["minimum_nonmissing_fraction"].iloc[0]) + 0.999))
        score = matrix.mean(axis=1, skipna=True)
        score[matrix.notna().sum(axis=1) < min_required] = pd.NA
        out[f"module_score__{module_id}"] = score
    return out


def main(cohort_path: str, module_path: str, output_path: str) -> None:
    df = pd.read_parquet(cohort_path)
    modules = pd.read_csv(module_path)
    scores = build_module_scores(df, modules)
    pd.concat([df, scores], axis=1).to_parquet(output_path, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", default="derived_data/analysis_cohort.parquet")
    parser.add_argument("--modules", default="tables/module_definitions_template.csv")
    parser.add_argument("--out", default="derived_data/analysis_cohort_with_scores.parquet")
    args = parser.parse_args()
    main(args.cohort, args.modules, args.out)
