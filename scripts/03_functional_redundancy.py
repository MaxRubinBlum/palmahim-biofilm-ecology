#!/usr/bin/env python3
"""Functional redundancy/core-guild table for Supplementary Table 5.

Inputs
------
--abundance : samples x MAG relative-abundance matrix (fractions, not percentages)
--traits    : MAG rows; one or more 0/1 functional columns
--taxonomy  : MAG rows; must contain Order

Outputs one row per function with carrier MAG/order counts, min/median/max summed
abundance, number of samples >=1%, core status, and median inverse-Simpson
effective MAG number.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def read(path, index_col=0):
    return pd.read_csv(path, sep="\t" if path.endswith((".tsv", ".txt")) else ",", index_col=index_col)

def effective_mag_number(x):
    x = np.asarray(x, dtype=float); x = x[x > 0]
    if x.size == 0: return np.nan
    p = x / x.sum(); return 1.0 / np.sum(p**2)

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--abundance",required=True); ap.add_argument("--traits",required=True)
    ap.add_argument("--taxonomy",required=True); ap.add_argument("--out",required=True); ap.add_argument("--threshold",type=float,default=0.01)
    args=ap.parse_args()
    A=read(args.abundance); T=read(args.traits); tax=read(args.taxonomy)
    mags=[m for m in A.columns if m in T.index and m in tax.index]; A=A[mags].astype(float); T=T.loc[mags]; tax=tax.loc[mags]
    rows=[]
    for trait in T.columns:
        carriers=T[trait].fillna(0).astype(bool)
        cmags=T.index[carriers].tolist(); X=A[cmags].to_numpy() if cmags else np.zeros((len(A),0))
        summed=X.sum(axis=1); eff=[effective_mag_number(r) for r in X]
        rows.append({"function":trait,"carrier_MAGs":len(cmags),"carrier_orders":tax.loc[cmags,"Order"].nunique() if cmags else 0,
                     "min_abundance":summed.min() if len(summed) else np.nan,"median_abundance":np.median(summed) if len(summed) else np.nan,
                     "max_abundance":summed.max() if len(summed) else np.nan,"samples_ge_1pct":int((summed>=args.threshold).sum()),
                     "core_ge_1pct_all_samples":bool(np.all(summed>=args.threshold)),"median_effective_MAG_number":float(np.nanmedian(eff)) if len(eff) else np.nan})
    Path(args.out).parent.mkdir(parents=True, exist_ok=True); pd.DataFrame(rows).to_csv(args.out,index=False)
if __name__=="__main__": main()
