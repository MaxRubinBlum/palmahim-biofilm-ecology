#!/usr/bin/env python3
"""Prepare Fig. 2a top-order abundance table.

Data required: sample x MAG relative abundance matrix and MAG taxonomy with Order.
The final heatmap clustering/ordering can be applied in the plotting program; this
script records the abundance aggregation and top-20 selection used by the figure.
"""
import argparse
from pathlib import Path
import pandas as pd

def read(p,index_col=0): return pd.read_csv(p,sep="\t" if p.endswith((".tsv",".txt")) else ",",index_col=index_col)
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--abundance",required=True); ap.add_argument("--taxonomy",required=True); ap.add_argument("--out",required=True); ap.add_argument("--top-n",type=int,default=20); a=ap.parse_args()
    A=read(a.abundance); T=read(a.taxonomy); mags=[m for m in A.columns if m in T.index]; X=A[mags].T.join(T[["Order"]]); O=X.groupby("Order").sum(numeric_only=True); top=O.sum(axis=1).nlargest(a.top_n).index; Path(a.out).parent.mkdir(parents=True,exist_ok=True); O.loc[top].to_csv(a.out)
if __name__=="__main__": main()
