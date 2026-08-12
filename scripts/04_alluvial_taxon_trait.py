#!/usr/bin/env python3
"""Prepare Fig. 3 abundance-weighted taxon-trait links.

Input table: one row per MAG, containing MAG, Order, abundance, and binary trait
columns. Abundance should be the manuscript's cumulative/mean MAG abundance
weight used for Fig. 3.

The script ranks taxa by cumulative abundance, keeps the top N orders, collapses
all others to 'Other taxa', calculates two-sided Fisher tests and phi on the
MAG-level binary taxon/trait matrix, applies BH correction over all tests, and
exports positive links with sqrt(abundance) ribbon width and phi opacity.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--input",required=True); ap.add_argument("--traits",required=True,nargs="+")
    ap.add_argument("--out",required=True); ap.add_argument("--top-n",type=int,default=15); ap.add_argument("--order-col",default="Order"); ap.add_argument("--abundance-col",default="abundance")
    a=ap.parse_args(); sep="\t" if a.input.endswith((".tsv",".txt")) else ","; df=pd.read_csv(a.input,sep=sep)
    totals=df.groupby(a.order_col)[a.abundance_col].sum().sort_values(ascending=False); top=set(totals.head(a.top_n).index)
    df["display_taxon"]=df[a.order_col].where(df[a.order_col].isin(top),"Other taxa")
    taxa=list(totals.head(a.top_n).index)+(["Other taxa"] if (~df[a.order_col].isin(top)).any() else [])
    rows=[]
    for taxon in taxa:
        in_taxon=df["display_taxon"].eq(taxon)
        for trait in a.traits:
            has=df[trait].fillna(0).astype(bool)
            A=int((in_taxon&has).sum()); B=int((in_taxon&~has).sum()); C=int((~in_taxon&has).sum()); D=int((~in_taxon&~has).sum())
            _,p=fisher_exact([[A,B],[C,D]],alternative="two-sided")
            den=np.sqrt((A+B)*(C+D)*(A+C)*(B+D)); phi=((A*D)-(B*C))/den if den else np.nan
            abundance=df.loc[in_taxon&has,a.abundance_col].sum()
            rows.append({"taxon":taxon,"trait":trait,"a":A,"b":B,"c":C,"d":D,"phi":phi,"p":p,"abundance":abundance})
    out=pd.DataFrame(rows); out["q"]=multipletests(out["p"],method="fdr_bh")[1]
    out["ribbon_width"]=np.sqrt(out["abundance"]); out["ribbon_opacity"]=out["phi"].clip(lower=0)
    out["positive_for_display"]=out["phi"]>0
    Path(a.out).parent.mkdir(parents=True,exist_ok=True); out.to_csv(a.out,index=False)
if __name__=="__main__": main()
