"""Independent HC1 linear and Breslow Cox numerical checks on an approved cohort.

Only model-level diagnostics are written. The input cohort is not distributed.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from statsmodels.regression.linear_model import OLS
from statsmodels.duration.hazard_regression import PHReg
import analyse

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cohort',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--checks',choices=['all','linear','cox'],default='all')
    args=parser.parse_args()
    df=pd.read_csv(args.cohort)
    if not df.participant_id.is_unique:raise ValueError('Participant IDs must be unique')
    df=df.loc[analyse.covariates(df).notna().all(axis=1)].copy()
    cov=analyse.covariates(df);x=pd.concat([df.chip_any.rename('exposure'),cov],axis=1)
    mask=df.chip_any.isin([0,1])&x.notna().all(axis=1)
    linear=[]
    for name,y in (analyse.domain_scores(df,mask).items() if args.checks!='cox' else []):
        ok=y.notna();design=x.loc[y.index[ok]]
        result=analyse.fit_linear(y,design)
        independent=OLS(y.loc[ok].to_numpy(),design.to_numpy()).fit(cov_type='HC1',use_t=True)
        diff=abs(result['beta']-float(independent.params[0]))
        sediff=abs(result['se']-float(independent.bse[0]))
        if diff>1e-7 or sediff>1e-7:raise AssertionError((name,diff,sediff))
        linear.append({'module':name,'n':result['n'],'beta_difference':diff,'se_difference':sediff})
    survival=[];cov=analyse.covariates(df,intercept=False)
    for outcome,predictor in ([('CHD','chip_any'),('HF','chip_large10'),('CVD_mortality','chip_any')] if args.checks!='linear' else []):
        design=pd.concat([df[[predictor]],cov],axis=1)
        primary,data=analyse.fit_cox(df,outcome,design)
        strict=analyse.cox_breslow_fit(data.time,data.event,data[primary['columns']],maxiter=500,ftol=1e-13,refine=True)
        if not strict['converged']:raise RuntimeError(strict['message'])
        independent=PHReg(data.time.to_numpy(),data[primary['columns']].to_numpy(),status=data.event.to_numpy(),ties='breslow').fit(disp=0,maxiter=100,tol=1e-8)
        bdiff=float(np.max(np.abs(strict['beta']-independent.params)))
        sediff=float(np.max(np.abs(strict['se']-independent.bse)))
        if bdiff>2e-5 or sediff>1e-6:raise AssertionError((outcome,bdiff,sediff))
        survival.append({'outcome':outcome,'predictor':predictor,'n':primary['n'],'events':primary['events'],'strict_beta_difference':bdiff,'strict_se_difference':sediff,'primary_difference_in_se_units':float(np.max(np.abs(primary['beta']-independent.params)/independent.bse))})
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps({'linear':linear,'cox':survival},indent=2),encoding='utf-8')
    print('Independent numerical checks passed.')

if __name__=='__main__':main()
