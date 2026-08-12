#!/usr/bin/env python3
"""Generate gene/MAG TPM and an integrated Figs. 6–7 plotting table.

Inputs:
--genes: gene_id, MAG, length_bp, read_count plus optional functional columns
--dna: MAG-indexed table with DNA relative abundance (one or more sample columns)
--traits: MAG-indexed curated trait table

Outputs gene TPM, MAG TPM, and MAG-level joined matrix. MAGs contributing >=1%
of total transcript abundance can be flagged for the manuscript plotting subset.
"""
import argparse
from pathlib import Path
import pandas as pd

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--genes",required=True); ap.add_argument("--dna",required=True); ap.add_argument("--traits",required=True); ap.add_argument("--outdir",required=True); ap.add_argument("--threshold-pct",type=float,default=1.0)
    a=ap.parse_args(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    sep=lambda p:"\t" if p.endswith((".tsv",".txt")) else ","
    g=pd.read_csv(a.genes,sep=sep(a.genes)); dna=pd.read_csv(a.dna,sep=sep(a.dna),index_col=0); traits=pd.read_csv(a.traits,sep=sep(a.traits),index_col=0)
    g["length_kb"]=g.length_bp/1000; g["RPK"]=g.read_count/g.length_kb; scale=g.RPK.sum()/1e6; g["TPM"]=g.RPK/scale
    g.to_csv(out/"gene_tpm.tsv",sep="\t",index=False)
    m=g.groupby("MAG",as_index=True).TPM.sum().to_frame("MAG_TPM"); m["relative_MAG_TPM_pct"]=100*m.MAG_TPM/m.MAG_TPM.sum(); m["include_ge_threshold"]=m.relative_MAG_TPM_pct>=a.threshold_pct
    integrated=m.join(dna,how="left").join(traits,how="left"); integrated.to_csv(out/"mag_integrated_expression_traits.tsv",sep="\t")
if __name__=="__main__": main()
