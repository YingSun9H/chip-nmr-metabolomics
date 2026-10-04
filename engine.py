"""Regression, multiple-testing and survival-model calculations."""
from __future__ import annotations
import math,warnings
import numpy as np
import pandas as pd
from scipy import optimize,stats
MAXITER=500


def rank_inverse_normal(values: np.ndarray) -> np.ndarray:
    s = pd.Series(values)
    lo = s.quantile(0.005)
    hi = s.quantile(0.995)
    x = s.clip(lo, hi).to_numpy(dtype=float)
    ranks = stats.rankdata(x, method="average")
    q = (ranks - 0.5) / len(ranks)
    return stats.norm.ppf(q)

def bh_adjust(p: pd.Series) -> pd.Series:
    pvals = p.to_numpy(dtype=float)
    out = np.full_like(pvals, np.nan, dtype=float)
    valid = np.isfinite(pvals)
    pv = pvals[valid]
    n = len(pv)
    if n == 0:
        return pd.Series(out, index=p.index)
    order = np.argsort(pv)
    ranked = pv[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    tmp = np.empty(n)
    tmp[order] = q
    out[valid] = tmp
    return pd.Series(out, index=p.index)

def robust_linear_fit(y: np.ndarray, x: np.ndarray, exposure_col: int) -> tuple[float, float, float]:
    beta, _, _, _ = np.linalg.lstsq(x, y, rcond=None)
    resid = y - x @ beta
    n, k = x.shape
    xtx_inv = np.linalg.pinv(x.T @ x)
    meat = x.T @ ((resid ** 2)[:, None] * x)
    vcov = (n / max(n - k, 1)) * (xtx_inv @ meat @ xtx_inv)
    se = math.sqrt(max(vcov[exposure_col, exposure_col], 0))
    b = float(beta[exposure_col])
    t = b / se if se > 0 else np.nan
    p = float(2 * stats.t.sf(abs(t), df=max(n - k, 1))) if np.isfinite(t) else np.nan
    return b, se, p

def observed_information(
    x: np.ndarray,
    e: np.ndarray,
    beta: np.ndarray,
    starts: np.ndarray,
    ends_e: np.ndarray,
    d: np.ndarray,
    chunk_size: int = 12000,
) -> np.ndarray:
    eta = np.clip(x @ beta, -50, 50)
    exp_eta = np.exp(eta)
    cum_risk = np.cumsum(exp_eta)
    cum_x = np.cumsum(exp_eta[:, None] * x, axis=0)
    risk = cum_risk[ends_e]
    mean_x = cum_x[ends_e] / risk[:, None]

    p = x.shape[1]
    info = np.zeros((p, p), dtype=float)
    flat_cum = np.zeros(p * p, dtype=float)
    event_ptr = 0
    n = x.shape[0]
    for start in range(0, n, chunk_size):
        stop = min(start + chunk_size, n)
        x_chunk = x[start:stop]
        w_chunk = exp_eta[start:stop]
        xx = np.einsum("ni,nj,n->nij", x_chunk, x_chunk, w_chunk, optimize=True).reshape(stop - start, p * p)
        csum = np.cumsum(xx, axis=0) + flat_cum
        while event_ptr < len(ends_e) and ends_e[event_ptr] < stop:
            if ends_e[event_ptr] >= start:
                s2 = csum[ends_e[event_ptr] - start].reshape(p, p)
                mx = mean_x[event_ptr]
                info += d[event_ptr] * (s2 / risk[event_ptr] - np.outer(mx, mx))
            event_ptr += 1
        flat_cum = csum[-1]
    return info

def cox_breslow_fit(time_s: pd.Series, event_s: pd.Series, x_df: pd.DataFrame, maxiter: int = MAXITER, ftol: float = 1e-8, refine: bool = False) -> dict:
    data = pd.concat(
        [
            pd.to_numeric(time_s, errors="coerce").rename("time"),
            pd.to_numeric(event_s, errors="coerce").rename("event"),
            x_df,
        ],
        axis=1,
    ).replace([np.inf, -np.inf], np.nan)
    data = data[data["time"].gt(0) & data["event"].isin([0, 1])].dropna()
    if data.empty or int(data["event"].sum()) == 0:
        raise ValueError("No valid incident events")
    x_cols = [c for c in data.columns if c not in {"time", "event"} and data[c].nunique(dropna=True) > 1]
    data = data[["time", "event", *x_cols]]
    if not x_cols:
        raise ValueError("No model covariates with variation")
    order = np.argsort(-data["time"].to_numpy(float), kind="mergesort")
    t = data["time"].to_numpy(float)[order]
    e = data["event"].to_numpy(float)[order]
    x = data[x_cols].to_numpy(float)[order]
    starts = np.r_[0, np.flatnonzero(t[1:] != t[:-1]) + 1]
    ends = np.r_[starts[1:] - 1, len(t) - 1]
    d_all = np.add.reduceat(e, starts)
    event_groups = d_all > 0
    starts_e = starts[event_groups]
    ends_e = ends[event_groups]
    d = d_all[event_groups]
    event_x_sum = np.add.reduceat(x * e[:, None], starts, axis=0)[event_groups]

    def objective(beta: np.ndarray) -> tuple[float, np.ndarray]:
        eta = x @ beta
        eta = np.clip(eta, -50, 50)
        exp_eta = np.exp(eta)
        cum_risk = np.cumsum(exp_eta)
        cum_x = np.cumsum(exp_eta[:, None] * x, axis=0)
        risk = cum_risk[ends_e]
        mean_x = cum_x[ends_e] / risk[:, None]
        event_eta_sum = np.add.reduceat(eta * e, starts)[event_groups]
        loglik = event_eta_sum.sum() - np.sum(d * np.log(risk))
        grad = event_x_sum.sum(axis=0) - np.sum(d[:, None] * mean_x, axis=0)
        return -float(loglik), -grad

    beta0 = np.zeros(x.shape[1])
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fit = optimize.minimize(
            lambda b: objective(b),
            beta0,
            method="L-BFGS-B",
            jac=True,
            options={"maxiter": maxiter, "ftol": ftol, "gtol": 1e-5, "maxls": 30},
        )
    if refine:
        # Optional Newton refinement for independent numerical validation.
        for _ in range(3):
            old_value, gradient = objective(fit.x)
            information = observed_information(x, e, fit.x, starts, ends_e, d)
            step = np.linalg.pinv(information) @ gradient
            candidate = fit.x - step
            if objective(candidate)[0] <= old_value + 1e-8:
                fit.x = candidate
            else:
                break
            if np.max(np.abs(step)) < 1e-8:
                break
    info = observed_information(x, e, fit.x, starts, ends_e, d)
    try:
        cov = np.linalg.pinv(info)
        se = np.sqrt(np.clip(np.diag(cov), 0, np.inf))
    except np.linalg.LinAlgError:
        se = np.full(x.shape[1], np.nan)
    neg_loglik, _ = objective(fit.x)
    z = fit.x / se
    p = 2 * stats.norm.sf(np.abs(z))
    return {
        "n": int(len(data)),
        "events": int(e.sum()),
        "columns": x_cols,
        "beta": fit.x,
        "se": se,
        "cov": cov,
        "p": p,
        "loglik": -float(neg_loglik),
        "converged": bool(fit.success),
        "message": str(fit.message),
    }

def rcs_basis(x: pd.Series | np.ndarray, knots: np.ndarray) -> pd.DataFrame:
    x_arr = np.asarray(x, dtype=float)
    k = np.asarray(knots, dtype=float)
    cols = {"score_linear": x_arr}
    last = k[-1]
    prev = k[-2]
    for j in range(len(k) - 2):
        kj = k[j]
        term = np.maximum(x_arr - kj, 0) ** 3
        term -= np.maximum(x_arr - prev, 0) ** 3 * (last - kj) / (last - prev)
        term += np.maximum(x_arr - last, 0) ** 3 * (prev - kj) / (last - prev)
        cols[f"score_rcs{j + 1}"] = term
    return pd.DataFrame(cols)

def ph_diagnostics(data,fit):
 # Breslow risk-set means include all participants tied at an event time.
 cols=fit['columns']; order=np.argsort(-data.time.to_numpy(),kind='mergesort')
 x=data[cols].to_numpy(float)[order]; t=data.time.to_numpy(float)[order]; e=data.event.to_numpy(float)[order]
 starts=np.r_[0,np.flatnonzero(t[1:]!=t[:-1])+1]; ends=np.r_[starts[1:]-1,len(t)-1]
 weights=np.exp(np.clip(x@fit['beta'],-50,50)); risk=np.cumsum(weights)
 rx=np.cumsum(x*weights[:,None],axis=0)
 means=rx[ends]/risk[ends,None]
 residual=x-np.repeat(means,np.diff(np.r_[starts,len(t)]),axis=0)
 event=e==1; residual=residual[event]; et=t[event]
 scaled=len(et)*residual@fit['cov']
 g=np.log(et); g-=g.mean(); u=g@scaled
 denominator=len(et)*np.sum(g*g)
 statistic=u*u/(denominator*np.diag(fit['cov']))
 result=[{'term':c,'chisq':float(q),'df':1,'p_value':float(stats.chi2.sf(q,1))} for c,q in zip(cols,statistic)]
 global_stat=float(u@np.linalg.pinv(denominator*fit['cov'])@u)
 result.append({'term':'GLOBAL','chisq':global_stat,'df':len(cols),'p_value':float(stats.chi2.sf(global_stat,len(cols)))})
 return result
