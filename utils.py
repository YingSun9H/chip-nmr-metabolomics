"""Run module, clone-size, driver-gene and subgroup models."""

from __future__ import annotations

import argparse
import pandas as pd
import statsmodels.api as sm
from utils import benjamini_hochberg


def fit_model(df: pd.DataFrame, outcome: str, exposure: str, covariates: list[str]) -> dict:
    model_df = df[[outcome, exposure, *covariates]].dropna()
    y = model_df[outcome]
    x = sm.add_constant(model_df[[exposure, *covariates]], has_constant="add")
    result = sm.OLS(y, x).fit(cov_type="HC1")
    return {"module_score": outcome, "exposure": exposure, "n": int(result.nobs), "beta": result.params[exposure], "se": result.bse[exposure], "p_value": result.pvalues[exposure]}


def main(input_path: str, output_path: str) -> None:
    df = pd.read_parquet(input_path)
    module_scores = [c for c in df.columns if c.startswith("module_score__")]
    exposures = ["chip_any", "chip_large_clone"]
    covariates = [c for c in df.columns if c.startswith("covariate__")]
    rows = [fit_model(df, score, exposure, covariates) for exposure in exposures for score in module_scores]
    results = pd.DataFrame(rows)
    results["fdr"] = results.groupby("exposure")["p_value"].transform(benjamini_hochberg)
    results.to_csv(output_path, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default="derived_data/analysis_cohort_with_scores.parquet")
    parser.add_argument("--out", default="results/module_score_results.csv")
    args = parser.parse_args()
    main(args.input, args.out)
