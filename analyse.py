"""CHIP–NMR metabolomics and five cardiovascular endpoints.

Input is an approved-researcher cohort export conforming to data_schema.csv.
Individual-level inputs and intermediates are never part of the public release.
"""
from __future__ import annotations
import argparse,json,math
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from engine import rank_inverse_normal,robust_linear_fit,bh_adjust,cox_breslow_fit,rcs_basis,ph_diagnostics

PCS=[f'pc{i}' for i in range(1,11)]
OUTCOMES=['CHD','stroke','HF','AF','CVD_mortality']
EXPOSURES=['chip_any','chip_large10','DNMT3A','TET2','ASXL1','JAK2','chip_DDR','chip_spliceosome']
ROOT=Path(__file__).resolve().parent
MODULES=pd.read_csv(ROOT/'metadata/module_membership.csv')
MANIFEST=pd.read_csv(ROOT/'metadata/biomarker_manifest.csv')

def covariates(df,clinical=False,intercept=True):
 x=pd.DataFrame(index=df.index)
 if intercept:x['intercept']=1.
 for col in ['age','sex',*PCS,'bmi','townsend']:
  x[col]=pd.to_numeric(df[col],errors='coerce')
 x['age_squared']=x.age**2
 for col,levels in [('smoking',[0,1,2]),('alcohol',[1,2,3])]:
  values=pd.to_numeric(df[col],errors='coerce')
  invalid=values.notna()&~values.isin(levels)
  if invalid.any():raise ValueError(f'Unexpected {col} category: expected {levels}')
  for level in levels[1:]:
   x[f'{col}_{level}']=values.eq(level).astype(float).where(values.notna())
 if clinical:
  for col in ['prevalent_cvd','prevalent_diabetes','prevalent_cancer','lipid_lowering','antihypertensive','egfr']:
   x[col]=pd.to_numeric(df[col],errors='coerce')
 if not intercept:
  # Preserve the design-column order used for the manuscript Cox fits.
  x=x[['age','age_squared',*PCS,'bmi','townsend','sex','smoking_1','smoking_2','alcohol_2','alcohol_3',*[c for c in x.columns if c not in ['age','age_squared',*PCS,'bmi','townsend','sex','smoking_1','smoking_2','alcohol_2','alcohol_3']]]]
  for col in x:
   if x[col].nunique()>2:x[col]=(x[col]-x[col].mean())/x[col].std()
 return x

def exposure_data(df,name):
 mask=df.chip_any.isin([0,1])
 if name not in ['chip_any','chip_large10']:
  mask &= (df[name]==1)|(df.chip_any==0)
 return mask,pd.to_numeric(df[name],errors='coerce')

def fit_linear(y,x,index=0):
 data=pd.concat([y.rename('response'),x],axis=1).replace([np.inf,-np.inf],np.nan).dropna()
 if len(data)<=len(x.columns)+2 or data.iloc[:,index+1].nunique()<2:return None
 b,se,p=robust_linear_fit(data.response.to_numpy(float),data[x.columns].to_numpy(float),index)
 ci=stats.t.ppf(.975,len(data)-len(x.columns))*se
 return dict(n=len(data),n_cases=int((data.iloc[:,index+1]==1).sum()),beta=b,se=se,p_value=p,ci95_low=b-ci,ci95_high=b+ci)

def mwas(df,assessment='baseline',clinical=False,restricted=False):
 if restricted:df=df[(df.prevalent_cvd==0)&(df.prevalent_diabetes==0)&(df.prevalent_cancer==0)]
 rows=[]
 cov=covariates(df,clinical)
 if assessment=='repeat':df=df[df.repeat_interval.notna()];cov=cov.loc[df.index]
 for exposure in EXPOSURES:
  mask,v=exposure_data(df,exposure)
  x=pd.concat([v.rename('exposure'),cov],axis=1)
  base=mask&x.notna().all(axis=1)
  for row in MANIFEST.itertuples():
   col=assessment+'__'+row.metabolite_id
   if col not in df:continue
   ok=base&df[col].notna()
   if ok.sum()<50:continue
   raw=df.loc[ok,col]
   z=pd.Series(rank_inverse_normal(raw.to_numpy(float)),index=raw.index)
   result=fit_linear(z,x.loc[ok])
   if result:rows.append(dict(exposure=exposure,assessment=assessment,metabolite_id=row.metabolite_id,metabolite=row.metabolite,**result))
 result=pd.DataFrame(rows)
 result['fdr_bh']=result.groupby(['exposure','assessment']).p_value.transform(bh_adjust)
 return result

def anchored_transform(baseline,values):
 ref=pd.to_numeric(baseline,errors='coerce').dropna()
 if len(ref)<10:return pd.Series(np.nan,index=values.index)
 lo,hi=ref.quantile([.005,.995]);ref=np.sort(ref.clip(lo,hi).to_numpy(float))
 raw=pd.to_numeric(values,errors='coerce');good=raw.notna();out=pd.Series(np.nan,index=values.index)
 v=raw[good].clip(lo,hi).to_numpy(float)
 p=(np.searchsorted(ref,v,'left')+np.searchsorted(ref,v,'right'))/(2*len(ref))
 out.loc[good]=stats.norm.ppf(p)
 return out

def domain_scores(df,mask):
 scores={}
 for module,g in MODULES.groupby('module_id',sort=False):
  transformed={}
  for mid in g.metabolite_id:
   raw=df.loc[mask,'baseline__'+mid];good=raw.notna()
   transformed[mid]=pd.Series(np.nan,index=raw.index)
   if good.sum()>=10:transformed[mid].loc[good]=rank_inverse_normal(raw.loc[good].to_numpy(float))
  z=pd.DataFrame(transformed);v=z.mean(axis=1)
  v[z.notna().sum(axis=1)<int(g.min_nonmissing_metabolites_for_score.iloc[0])]=np.nan
  scores[module]=v
 return scores

def module_exposure_data(df,name):
 # Clone-size module contrasts use participants without CHIP as controls.
 if name=='chip_small':
  value=((df.chip_any==1)&(df.chip_large10==0)).astype(float)
  return (df.chip_any==0)|(value==1),value
 if name=='chip_large10':
  return (df.chip_any==0)|(df.chip_large10==1),df.chip_large10
 return exposure_data(df,name)

def module_models(df):
 cov=covariates(df);rows=[]
 names=[*EXPOSURES,'chip_small']
 for exposure in names:
  mask,v=module_exposure_data(df,exposure)
  x=pd.concat([v.rename('exposure'),cov],axis=1);mask &= x.notna().all(axis=1)
  for module,score in domain_scores(df,mask).items():
   fit=fit_linear(score,x.loc[mask])
   if fit:rows.append(dict(exposure=exposure,module_id=module,**fit))
 out=pd.DataFrame(rows);out['fdr_bh']=out.groupby('exposure').p_value.transform(bh_adjust)
 return out

def signature_scores(df,weights):
 base_raw=pd.DataFrame(index=df.index);repeat_raw=pd.DataFrame(index=df.index);baseline_standardised=pd.DataFrame(index=df.index)
 for h in weights.itertuples():
  mid=h.metabolite_id;raw=df['baseline__'+mid]
  z=anchored_transform(raw,raw)
  base_raw[mid]=z
  # Use the primary survival score's rank transform and array standardisation.
  good=raw.notna();a=rank_inverse_normal(raw.loc[good].to_numpy(float))
  baseline_standardised[mid]=pd.Series(np.nan,index=df.index)
  baseline_standardised.loc[good,mid]=(a-np.nanmean(a))/np.nanstd(a,ddof=1)
  repeat_raw[mid]=anchored_transform(raw,df['repeat__'+mid]) if 'repeat__'+mid in df else np.nan
 w=weights.set_index('metabolite_id').beta
 nmin=math.ceil(len(w)*.5)
 b=base_raw.mul(w,axis=1).sum(axis=1);r=repeat_raw.mul(w,axis=1).sum(axis=1)
 b[base_raw.notna().sum(axis=1)<nmin]=np.nan;r[repeat_raw.notna().sum(axis=1)<nmin]=np.nan
 mean=b.mean();sd=b.std()
 classes=['lipoprotein_lipids_particles','fatty_acids','lipids_other','glycolysis_energy','ketone_bodies','amino_acids','inflammation_protein_renal','other_or_unclassified']
 ordered=weights.merge(MANIFEST[['metabolite_id','metabolic_class']],on='metabolite_id',how='left')
 ordered['class_order']=ordered.metabolic_class.map({c:i for i,c in enumerate(classes)}).fillna(99)
 mids=ordered.sort_values(['class_order','fdr_bh','metabolite_id']).metabolite_id.tolist()
 clinical=baseline_standardised[mids].mul(w.reindex(mids),axis=1).sum(axis=1)
 clinical[baseline_standardised.notna().sum(axis=1)<nmin]=np.nan
 clinical=(clinical-clinical.mean())/clinical.std()
 return clinical,(b-mean)/sd,(r-mean)/sd

def conventional_linear(y,x,index=0):
 # Repeat signature/module models use residual-variance OLS covariance.
 data=pd.concat([y.rename('response'),x],axis=1).replace([np.inf,-np.inf],np.nan).dropna()
 columns=[c for c in x.columns if c=='intercept' or data[c].nunique()>1]
 term=x.columns[index]
 if term not in columns:raise ValueError('Exposure has no variation')
 a=data[columns].to_numpy(float);v=data.response.to_numpy(float)
 coef=np.linalg.lstsq(a,v,rcond=None)[0];resid=v-a@coef;dof=len(v)-a.shape[1]
 vcov=(resid@resid/dof)*np.linalg.pinv(a.T@a);j=columns.index(term)
 b=float(coef[j]);se=float(np.sqrt(max(vcov[j,j],0)));p=float(2*stats.t.sf(abs(b/se),dof))
 return dict(n=len(v),n_cases=int(data[term].eq(1).sum()),n_controls=int(data[term].eq(0).sum()),
             beta=b,se=se,p_value=p,ci95_low=b-1.96*se,ci95_high=b+1.96*se)

def paired_biomarker_models(df):
 # Figure 3B: joint small/large indicators, no-CHIP reference, HC1 covariance.
 # Change uses clipped concentrations on the baseline mean/SD scale.
 cov=covariates(df);rows=[]
 small=((df.chip_any==1)&(df.chip_large10==0)).astype(float)
 x=pd.concat([cov[['intercept']],small.rename('chip_small'),df.chip_large10,
              cov.drop(columns='intercept')],axis=1)
 eligible=(df.chip_any==0)|((df.chip_any==1)&df.chip_large10.isin([0,1]))
 for h in MANIFEST.itertuples():
  mid=h.metabolite_id;bcol='baseline__'+mid;rcol='repeat__'+mid
  if rcol not in df or bcol not in df:continue
  b=df[bcol].where(eligible);r=df[rcol].where(eligible)
  lo,hi=b.quantile([.005,.995]);clipped=b.clip(lo,hi);sd=clipped.std(ddof=0)
  if not np.isfinite(sd) or sd<=0:continue
  delta=(r.clip(lo,hi)-clipped.mean())/sd-(clipped-clipped.mean())/sd
  complete=eligible&delta.notna()&x.notna().all(axis=1)
  if complete.sum()<500:continue
  data=x.loc[complete];values=delta.loc[complete]
  beta,se,p=robust_linear_fit(values.to_numpy(float),data.to_numpy(float),2)
  rows.append(dict(exposure='chip_large10',metabolite_id=mid,n=int(complete.sum()),
                   n_cases=int(df.loc[complete,'chip_large10'].eq(1).sum()),
                   n_small=int(small.loc[complete].eq(1).sum()),beta=beta,se=se,
                   p_value=p,ci95_low=beta-1.96*se,ci95_high=beta+1.96*se))
 out=pd.DataFrame(rows);out['fdr_bh']=bh_adjust(out.p_value)
 return out

def paired_models(df):
 cov=covariates(df);paired=df.repeat_interval.notna();module_rows=[]
 base={};repeat={}
 for mid in MODULES.metabolite_id.unique():
  b=df['baseline__'+mid];base[mid]=anchored_transform(b,b)
  repeat[mid]=anchored_transform(b,df['repeat__'+mid])
 bdf=pd.DataFrame(base);rdf=pd.DataFrame(repeat)
 for module,g in MODULES.groupby('module_id',sort=False):
  mids=list(g.metabolite_id);minimum=int(g.min_nonmissing_metabolites_for_score.iloc[0])
  b=bdf[mids].mean(axis=1);r=rdf[mids].mean(axis=1)
  b[bdf[mids].notna().sum(axis=1)<minimum]=np.nan;r[rdf[mids].notna().sum(axis=1)<minimum]=np.nan
  # Baseline module scale is retained for repeated measurements.
  mean=b.mean();sd=b.std();b=(b-mean)/sd;r=(r-mean)/sd
  delta=r-b;delta=(delta-delta.loc[paired].mean())/delta.loc[paired].std()
  df['module_baseline__'+module]=b;df['module_change__'+module]=delta
  for exposure in ['chip_any','chip_large10']:
   mask,v=module_exposure_data(df,exposure);mask &= paired
   # Match the repeat-module design: exposure, baseline, interval, covariates.
   x=pd.concat([pd.Series(1.,index=df.index,name='intercept'),v.rename(exposure),
                b.rename('baseline_value'),df.repeat_interval,cov.drop(columns='intercept')],axis=1)
   fit=conventional_linear(delta.loc[mask],x.loc[mask],1)
   module_rows.append(dict(exposure=exposure,module_id=module,**fit))
 m=pd.DataFrame(module_rows);m['fdr_bh']=m.groupby('exposure').p_value.transform(bh_adjust)
 return paired_biomarker_models(df),m

def signature_trajectory(df):
 b=df.paired_baseline_metabolic_signature_score;r=df.followup_metabolic_signature_score;delta=r-b
 rows=[]
 for label,mask in [('No CHIP',df.chip_any.eq(0)),('Any CHIP',df.chip_any.eq(1)),('Large CHIP',df.chip_large10.eq(1))]:
  valid=mask&b.notna()&r.notna();row=dict(group=label,n=int(valid.sum()))
  for name,value in [('baseline',b),('repeat',r),('delta',delta)]:
   a=value.loc[valid];mean=a.mean();se=a.std()/np.sqrt(len(a))
   row.update({name+'_mean':mean,name+'_ci95_low':mean-1.96*se,name+'_ci95_high':mean+1.96*se})
  rows.append(row)
 cov=covariates(df)
 for exposure,label in [('chip_any','Any CHIP'),('chip_large10','Large CHIP')]:
  # Large-clone indicator uses non-large participants, including small clones.
  x=pd.concat([cov[['intercept']],df[exposure],b.rename('baseline_score'),df.repeat_interval,
               cov.drop(columns='intercept')],axis=1)
  fit=conventional_linear(delta,x,1)
  row=next(row for row in rows if row['group']==label)
  row.update({'adjusted_change_'+k:fit[k] for k in ['beta','ci95_low','ci95_high','p_value','n']})
 return pd.DataFrame(rows)

def fit_cox(df,outcome,x,assessment='baseline',lag=0):
 status=df[f'{outcome}__{assessment}_status'];duration=df[f'{outcome}__{assessment}_time']
 mask=status.isin([0,1])&(duration>lag)
 data=pd.concat([(duration-lag).rename('time'),status.rename('event'),x],axis=1)
 data=data.loc[mask].replace([np.inf,-np.inf],np.nan).dropna()
 fit=cox_breslow_fit(data.time,data.event,data.drop(columns=['time','event']),maxiter=500)
 if not fit['converged']:raise RuntimeError(f'Cox fit did not converge: {outcome}, {assessment}, {fit["message"]}')
 return fit,data

def cox_row(fit,predictor):
 j=fit['columns'].index(predictor);b=fit['beta'][j];se=fit['se'][j]
 return dict(n=fit['n'],events=fit['events'],beta_log_hr=b,se=se,hr=np.exp(b),ci95_low=np.exp(b-1.96*se),ci95_high=np.exp(b+1.96*se),p_value=fit['p'][j])

def outcome_models(df,output,lag=0):
 cov=covariates(df,intercept=False);rows=[];diagnostics=[]
 for o in OUTCOMES:
  for exposure in ['chip_any','chip_large10']:
   x=pd.concat([df[exposure],cov],axis=1)
   fit,data=fit_cox(df,o,x,lag=lag)
   rows.append(dict(outcome=o,predictor=exposure,**cox_row(fit,exposure)))
   diagnostics.extend(dict(outcome=o,model=exposure,**r) for r in ph_diagnostics(data,fit))
 out=pd.DataFrame(rows);out['fdr_bh']=bh_adjust(out.p_value)
 out.to_csv(output/('chip_outcomes_lag.csv' if lag else 'chip_outcomes.csv'),index=False)
 pd.DataFrame(diagnostics).to_csv(output/('ph_lag.csv' if lag else 'ph_chip.csv'),index=False)
 return out

def spline_models(df,output):
 cov=covariates(df,intercept=False);score=df.metabolic_signature_score
 knots=score.quantile([.05,.35,.65,.95]).to_numpy();grid=np.linspace(-3,3,140)
 basis=rcs_basis(score,knots);basis.index=df.index
 gridbasis=rcs_basis(grid,knots);ref=rcs_basis(np.array([0.]),knots).iloc[0].to_numpy()
 curves=[];summaries=[];diagnostics=[]
 for o in OUTCOMES:
  x=pd.concat([basis,cov],axis=1);fit,data=fit_cox(df,o,x)
  # Null and linear fits use the identical complete-case sample.
  null=cox_breslow_fit(data.time,data.event,data[cov.columns],maxiter=500)
  linear=cox_breslow_fit(data.time,data.event,data[['score_linear',*cov.columns]],maxiter=500)
  if not null['converged'] or not linear['converged']:raise RuntimeError('Spline comparison fit did not converge')
  inds=[fit['columns'].index(c) for c in basis.columns]
  b=fit['beta'][inds];v=fit['cov'][np.ix_(inds,inds)];diff=gridbasis.to_numpy()-ref
  lp=diff@b;se=np.sqrt(np.einsum('ij,jk,ik->i',diff,v,diff))
  curves.extend(dict(outcome=o,score=g,hr=np.exp(l),ci95_low=np.exp(l-1.96*s),ci95_high=np.exp(l+1.96*s)) for g,l,s in zip(grid,lp,se))
  summaries.append(dict(outcome=o,n=fit['n'],events=fit['events'],knots=','.join(map(str,knots)),p_overall=stats.chi2.sf(max(0,2*(fit['loglik']-null['loglik'])),len(inds)),p_nonlinear=stats.chi2.sf(max(0,2*(fit['loglik']-linear['loglik'])),len(inds)-1)))
  diagnostics.extend(dict(outcome=o,model='signature_spline',**r) for r in ph_diagnostics(data,fit))
 pd.DataFrame(curves).to_csv(output/'signature_spline_curves.csv',index=False)
 pd.DataFrame(summaries).to_csv(output/'signature_spline_summary.csv',index=False)
 pd.DataFrame(diagnostics).to_csv(output/'ph_signature.csv',index=False)

def landmark_models(df,output):
 cov=covariates(df,intercept=False);rows=[];module_rows=[]
 for o in OUTCOMES:
  for predictor in ['followup_metabolic_signature_score','metabolic_signature_score_change']:
   x=pd.concat([df[predictor],cov],axis=1)
   if predictor.endswith('_change'):x=pd.concat([x,df[['paired_baseline_metabolic_signature_score','repeat_interval']]],axis=1)
   fit,data=fit_cox(df,o,x,assessment='repeat');rows.append(dict(outcome=o,predictor=predictor,**cox_row(fit,predictor)))
  for module in MODULES.module_id.drop_duplicates():
   pred='module_change__'+module
   x=pd.concat([df[pred],cov,df[['module_baseline__'+module,'repeat_interval']]],axis=1)
   fit,data=fit_cox(df,o,x,assessment='repeat');module_rows.append(dict(outcome=o,module_id=module,**cox_row(fit,pred)))
 a=pd.DataFrame(rows);a['fdr_bh']=a.groupby('predictor').p_value.transform(bh_adjust)
 m=pd.DataFrame(module_rows);m['fdr_bh']=bh_adjust(m.p_value)
 a.to_csv(output/'landmark_signature_outcomes.csv',index=False);m.to_csv(output/'landmark_module_outcomes.csv',index=False)

def module_pca(df,output):
 mask=covariates(df).notna().all(axis=1);means=domain_scores(df,mask);rows=[]
 for module,g in MODULES.groupby('module_id',sort=False):
  z=pd.DataFrame({mid:anchored_transform(df.loc[mask,'baseline__'+mid],df.loc[mask,'baseline__'+mid]) for mid in g.metabolite_id})
  score=means[module];ok=score.notna();z=z.loc[ok];z=z.fillna(z.mean());z-=z.mean()
  values,vec=np.linalg.eigh(np.cov(z.to_numpy(),rowvar=False));pc=z.to_numpy()@vec[:,-1]
  corr=float(np.corrcoef(pc,score.loc[ok])[0,1])
  rows.append(dict(module_id=module,n=len(z),n_biomarkers=len(g),pc1_variance_fraction=float(values[-1]/values.sum()),mean_pc1_correlation=abs(corr)))
 pd.DataFrame(rows).to_csv(output/'module_pca_support.csv',index=False)

SUBGROUPS=[
 ('Sex','Women','Men','sex',0,1),
 ('Age group','Age <60','Age >=60','age',None,60),
 ('BMI group','BMI <=25','BMI >25','bmi',None,25),
 ('Baseline CVD','Baseline CVD no','Baseline CVD yes','prevalent_cvd',0,1),
 ('Baseline diabetes','Baseline diabetes no','Baseline diabetes yes','prevalent_diabetes',0,1),
 ('Lipid-lowering medication','Lipid-lowering medication no','Lipid-lowering medication yes','lipid_lowering',0,1)]

def subgroup_indicator(df,spec):
 name,ref,cmp,col,low,high=spec
 value=pd.to_numeric(df[col],errors='coerce')
 if name=='Age group':out=value.ge(high).astype(float).where(value.notna())
 elif name=='BMI group':out=value.gt(high).astype(float).where(value.notna())
 else:out=value.eq(high).astype(float).where(value.isin([low,high]))
 return out

def subgroup_models(df,output):
 cov=covariates(df);rows=[]
 for spec in SUBGROUPS:
  indicator=subgroup_indicator(df,spec)
  for level,label in [(0,spec[1]),(1,spec[2])]:
   for exposure in ['chip_any']:
    mask,value=module_exposure_data(df,exposure)
    x=pd.concat([value.rename('exposure'),cov],axis=1)
    mask &= indicator.eq(level)&x.notna().all(axis=1)
    # Transformation is fitted within this exposure/stratum sample.
    # Retain the original primary-model design, including constant stratum
    # covariates, for its documented HC1 residual degrees of freedom.
    for module,score in domain_scores(df,mask).items():
     fit=fit_linear(score,x.loc[mask])
     if fit:rows.append(dict(subgroup_dimension=spec[0],subgroup=label,exposure=exposure,module_id=module,**fit))
 out=pd.DataFrame(rows)
 out['fdr_bh']=out.groupby(['subgroup_dimension','subgroup','exposure']).p_value.transform(bh_adjust)
 out.to_csv(output/'module_subgroup_associations.csv',index=False)
 return out

def module_interactions(df,output):
 cov=covariates(df);rows=[]
 for spec in SUBGROUPS:
  indicator=subgroup_indicator(df,spec);x=cov.copy()
  if spec[0]=='Sex':x=x.drop(columns=['sex'])
  x=pd.concat([(df.chip_any*indicator).rename('interaction'),df.chip_any,indicator.rename('subgroup'),x],axis=1)
  mask=df.chip_any.isin([0,1])&x.notna().all(axis=1)
  # Use the complete interaction-model sample to fit ingredient transforms.
  for module,score in domain_scores(df,mask).items():
   fit=fit_linear(score,x.loc[mask])
   if fit:
    eligible=score.notna()
    fit['n_comparison_chip']=fit.pop('n_cases')
    rows.append(dict(subgroup_dimension=spec[0],reference_subgroup=spec[1],comparison_subgroup=spec[2],module_id=module,
      n_chip=int(df.loc[score.index[eligible],'chip_any'].sum()),**fit))
 out=pd.DataFrame(rows);out['fdr_bh']=bh_adjust(out.p_value)
 out.to_csv(output/'module_interactions.csv',index=False)
 return out

def time_period_models(df,output):
 cov=covariates(df,intercept=False);rows=[]
 for o in OUTCOMES:
  status=df[o+'__baseline_status'];duration=df[o+'__baseline_time']
  for period in ['0-5 years','>5 years']:
   eligible=status.isin([0,1])
   if period=='0-5 years':t=duration.clip(upper=5);e=((status==1)&(duration<=5)).astype(float)
   else:eligible &= duration>5;t=duration-5;e=status
   for predictor in ['chip_any','chip_large10']:
    x=pd.concat([df[predictor],cov],axis=1)
    fit=cox_breslow_fit(t.loc[eligible],e.loc[eligible],x.loc[eligible],maxiter=500)
    if not fit['converged']:raise RuntimeError('Time-period model did not converge')
    rows.append(dict(outcome=o,period=period,predictor=predictor,**cox_row(fit,predictor)))
 out=pd.DataFrame(rows);out['fdr_bh']=out.groupby('period').p_value.transform(bh_adjust)
 out.to_csv(output/'time_period_sensitivity.csv',index=False)

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('stage',choices=['baseline','modules','repeat','outcomes','sensitivity','all'])
 p.add_argument('--cohort',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
 args=p.parse_args();args.output.mkdir(parents=True,exist_ok=True)
 df=pd.read_csv(args.cohort);assert df.participant_id.is_unique,'Participant IDs must be unique'
 primary=covariates(df).notna().all(axis=1);df=df.loc[primary].copy()
 audit={'primary_n':len(df),'repeat_n':int(df.repeat_interval.notna().sum()),'endpoints':OUTCOMES}
 (args.output/'cohort_counts.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
 if args.stage in ['baseline','all']:
  baseline=mwas(df);baseline.to_csv(args.output/'baseline_biomarker_associations.csv',index=False)
 elif args.stage in ['repeat','outcomes']:baseline=pd.read_csv(args.output/'baseline_biomarker_associations.csv')
 if args.stage in ['modules','all']:
  module_models(df).to_csv(args.output/'module_associations.csv',index=False);module_pca(df,args.output);subgroup_models(df,args.output);module_interactions(df,args.output)
 if args.stage in ['repeat','all']:mwas(df,'repeat').to_csv(args.output/'repeat_biomarker_associations.csv',index=False)
 if args.stage in ['repeat','outcomes','all']:
  weights=baseline[(baseline.exposure=='chip_any')&(baseline.fdr_bh<.05)]
  clinical,b,r=signature_scores(df,weights);df['metabolic_signature_score']=clinical;df['paired_baseline_metabolic_signature_score']=b;df['followup_metabolic_signature_score']=r
  df['metabolic_signature_score_change']=r-b
  a,m=paired_models(df);a.to_csv(args.output/'paired_biomarker_changes.csv',index=False);m.to_csv(args.output/'paired_module_changes.csv',index=False)
  signature_trajectory(df).to_csv(args.output/'signature_trajectory.csv',index=False)
 if args.stage in ['outcomes','all']:
  outcome_models(df,args.output);spline_models(df,args.output);landmark_models(df,args.output);time_period_models(df,args.output)
 if args.stage in ['sensitivity','all']:
  mwas(df,clinical=True).to_csv(args.output/'clinical_adjusted_biomarker_associations.csv',index=False)
  mwas(df,restricted=True).to_csv(args.output/'disease_free_biomarker_associations.csv',index=False)
  outcome_models(df,args.output,lag=2)
 print('Completed:',args.stage)

if __name__=='__main__':main()
