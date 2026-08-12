#!/usr/bin/env python3
"""Construct Fig. 4 C1 candidate network and Cytoscape tables.

Inputs
------
--relative  samples x MAG relative abundance matrix
--raw       samples x MAG raw abundance/coverage matrix used in edge score
--metadata  sample, habitat
--traits    MAG-indexed table with methane_oxidizer, methanol_consumer,
            formate_consumer, formaldehyde_processing
--taxonomy  MAG-indexed table containing Order (additional ranks optional)

Biological edge inclusion requires compatible C1 capacity and occurrence of both
MAGs in >=1 shared habitat. CLR correlations are supporting annotations, not
filters. For each source-target-interaction combination, the habitat with the
highest score is retained. The complete table is exported; a simple deterministic
visualization subset retains top-scoring edges globally and per source/target.
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests


def read(path,index_col=0): return pd.read_csv(path,sep="\t" if path.endswith((".tsv",".txt")) else ",",index_col=index_col)
def score(Rs,Rt,Cs,Ct,shared_fraction=None):
    return np.log(1+1e6*np.sqrt(Rs*Rt))*np.log(1+np.sqrt(Cs*Ct))*(0.65 if shared_fraction is None else 1+shared_fraction)
def compatible(s,t):
    out=[]
    if bool(s.get("methane_oxidizer",0)):
        if bool(t.get("methanol_consumer",0)): out.append("methane_to_methanol")
        if bool(t.get("formate_consumer",0)): out.append("methane_to_formate")
    src_down=bool(s.get("methanol_consumer",0) or s.get("formate_consumer",0) or s.get("formaldehyde_processing",0))
    if src_down and bool(t.get("formate_consumer",0)): out.append("methanol_to_formate")
    return out

def clr(X,pseudocount=1e-6):
    Y=np.log(np.asarray(X,float)+pseudocount); return Y-Y.mean(axis=1,keepdims=True)

def main():
    ap=argparse.ArgumentParser();
    for x in ("relative","raw","metadata","traits","taxonomy","outdir"): ap.add_argument(f"--{x}",required=True)
    ap.add_argument("--top-global",type=int,default=600); ap.add_argument("--top-per-source",type=int,default=10); ap.add_argument("--top-per-target",type=int,default=10)
    a=ap.parse_args(); out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    R=read(a.relative); C=read(a.raw); T=read(a.traits); Tax=read(a.taxonomy); meta=pd.read_csv(a.metadata,sep="\t" if a.metadata.endswith((".tsv",".txt")) else ",").set_index("sample")
    samples=[s for s in R.index if s in C.index and s in meta.index]; mags=[m for m in R.columns if m in C.columns and m in T.index and m in Tax.index]
    R=R.loc[samples,mags].astype(float); C=C.loc[samples,mags].astype(float); T=T.loc[mags]; Tax=Tax.loc[mags]; meta=meta.loc[samples]
    rows=[]
    for si,src in enumerate(mags):
        for ti,tgt in enumerate(mags):
            if src==tgt: continue
            kinds=compatible(T.loc[src],T.loc[tgt])
            if not kinds: continue
            for habitat, idx in meta.groupby("habitat").groups.items():
                ss=list(idx); rs=R.loc[ss,src]; rt=R.loc[ss,tgt]
                present_s=rs>0; present_t=rt>0
                if not present_s.any() or not present_t.any(): continue
                shared=(present_s & present_t); sf=float(shared.mean()) if shared.any() else None
                Rs=float(rs.sum()); Rt=float(rt.sum()); Cs=float(C.loc[ss,src].sum()); Ct=float(C.loc[ss,tgt].sum())
                for kind in kinds:
                    rows.append({"source":src,"target":tgt,"interaction":kind,"habitat":habitat,"Rs":Rs,"Rt":Rt,"Cs":Cs,"Ct":Ct,
                                 "shared_samples":int(shared.sum()),"habitat_replicates":len(ss),"shared_fraction":np.nan if sf is None else sf,
                                 "edge_score":score(Rs,Rt,Cs,Ct,sf)})
    edges=pd.DataFrame(rows)
    if edges.empty: raise SystemExit("No compatible shared-habitat edges found")
    edges=edges.sort_values(["source","target","interaction","edge_score","shared_samples"],ascending=[True,True,True,False,False])
    edges=edges.groupby(["source","target","interaction"],as_index=False).head(1).reset_index(drop=True)

    # Global CLR Spearman support; FDR is calculated across all unique source-target tests.
    Z=clr(R.to_numpy()); pos={m:i for i,m in enumerate(mags)}; pairs=edges[["source","target"]].drop_duplicates()
    cor=[]
    for _,r in pairs.iterrows():
        rho,p=spearmanr(Z[:,pos[r.source]],Z[:,pos[r.target]])
        cor.append({"source":r.source,"target":r.target,"clr_rho":rho,"clr_p":p})
    cor=pd.DataFrame(cor); cor["clr_q"]=multipletests(cor["clr_p"].fillna(1),method="fdr_bh")[1]; cor["clr_supported"]=(cor.clr_rho>=0.60)&(cor.clr_q<=0.05)
    edges=edges.merge(cor,on=["source","target"],how="left")
    edges.to_csv(out/"c1_edges_complete.tsv",sep="\t",index=False)

    # Deterministic display subset: high-scoring globally + best edges per source and target.
    keep=set(edges.nlargest(a.top_global,"edge_score").index)
    for _,g in edges.groupby("source"): keep.update(g.nlargest(a.top_per_source,"edge_score").index)
    for _,g in edges.groupby("target"): keep.update(g.nlargest(a.top_per_target,"edge_score").index)
    disp=edges.loc[sorted(keep)].copy().sort_values("edge_score",ascending=False)
    disp.to_csv(out/"c1_edges_cytoscape.tsv",sep="\t",index=False)

    # Node table uses all MAGs present in the displayed network.
    nmags=sorted(set(disp.source)|set(disp.target)); nodes=[]
    for m in nmags:
        dom = R[m].groupby(meta["habitat"]).mean().idxmax() if len(meta) else ""
        row={"MAG":m,"mean_relative_abundance":float(R[m].mean()),"dominant_habitat":dom}
        row.update(Tax.loc[m].to_dict());
        row["C1_guild"] = "methane_oxidizer" if bool(T.loc[m].get("methane_oxidizer",0)) else ("dual_consumer" if bool(T.loc[m].get("methanol_consumer",0)) and bool(T.loc[m].get("formate_consumer",0)) else ("methanol_consumer" if bool(T.loc[m].get("methanol_consumer",0)) else ("formate_consumer" if bool(T.loc[m].get("formate_consumer",0)) else "other")))
        row["formaldehyde_processing"]=bool(T.loc[m].get("formaldehyde_processing",0)); nodes.append(row)
    pd.DataFrame(nodes).to_csv(out/"c1_nodes_cytoscape.tsv",sep="\t",index=False)

if __name__=="__main__": main()
