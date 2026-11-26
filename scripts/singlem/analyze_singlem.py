#!/usr/bin/env python3
"""SingleM marker-based recovery, diversity and comparison with original MAGs."""
import argparse
import re
from pathlib import Path
from itertools import combinations
import numpy as np
import pandas as pd

KEY = ["sample", "gene", "sequence"]
GROUP = ["sample", "domain", "gene"]
CATS = ["MAG_linked", "Assembled_only", "Unaccounted"]
RANKS = [("phylum", "p", "phylum"), ("class", "c", "class"),
         ("order", "o", "Order"), ("family", "f", "family"), ("genus", "g", "genus")]
SAMPLES = sorted("AnemPM22 Bio1 Bio2 Bio3 Bio4 Bio6 Bio7 Bio8 Bio9 CrabPM21 EGGPM35 EGGPM36 PBB PBT PM27 PM28 PM29 PM30 PM39".split())

def substrate(s):
    if s.startswith("Bio"): return "Carbonate"
    if s.startswith("PM"): return "Plastic"
    if s.startswith("EGG"): return "Shark_egg"
    return {"PBB":"ARMS","PBT":"ARMS","CrabPM21":"Crab","AnemPM22":"Anemone"}[s]

def load(path):
    x = pd.read_csv(path, sep="\t", comment="#")
    x["sample"] = x["sample"].str.replace(r"_1$", "", regex=True)
    x["domain"] = x["taxonomy"].str.extract(r"(d__[^;]+)")[0]
    return x

def keys(x):
    return pd.MultiIndex.from_frame(x[KEY].drop_duplicates())

def save(x, path, index=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    x.to_csv(path, sep="\t", index=index)

def norm(x):
    assert (x.sum(axis=1)>0).all()
    return x.div(x.sum(axis=1),axis=0)

def recovery(raw, out):
    dest = out/"analysis_filtered"
    idx = keys(raw)
    for mode in ["exact","imperfect"]:
        source = out/"appraisal"/mode/"pooled"
        b,a,u,m = [idx.isin(keys(load(source/f"{name}.tsv")))
                   for name in ["binned","assembled","unbinned","unaccounted"]]
        assert not (b & ~a).any()
        assert np.array_equal(u,a & ~b)
        assert np.array_equal(m,~a & ~b)
        assert np.all(b.astype(int)+u.astype(int)+m.astype(int)==1)
        df = raw.copy()
        df["Recovery"] = np.select([b,u],CATS[:2],default=CATS[2])
        total = df.groupby(GROUP)["coverage"].transform("sum")
        assert (total>0).all()
        df["marker_fraction"] = df["coverage"]/total
        save(df,dest/f"OTU_recovery_{mode}.tsv.gz")
        marker = df.groupby(GROUP+["Recovery"])["marker_fraction"].sum().unstack(fill_value=0).reindex(columns=CATS,fill_value=0)
        save(marker,dest/f"marker_recovery_{mode}.tsv",True)
        sample = marker.groupby(["sample","domain"]).mean()*100
        sample["N_markers"] = marker.groupby(["sample","domain"]).size()
        save(sample,dest/f"sample_recovery_{mode}.tsv",True)
        n = df.groupby(["sample","domain"])["gene"].nunique()
        for rank,prefix,_ in RANKS:
            df["lineage"] = df["taxonomy"].str.extract(rf"^(.*?{prefix}__[^;]+)")[0].fillna(df["taxonomy"]+"; Unresolved_"+rank)
            tab = df.groupby(["sample","domain","lineage","Recovery"])["marker_fraction"].sum().unstack(fill_value=0).reindex(columns=CATS,fill_value=0).reset_index()
            tab = tab.join(n.rename("N_markers"),on=["sample","domain"])
            tab[CATS] = tab[CATS].div(tab["N_markers"],axis=0)*100
            tab["Community_pct"] = tab[CATS].sum(axis=1)
            tab["MAG_recovery_within_lineage_pct"] = 100*tab["MAG_linked"]/tab["Community_pct"]
            save(tab,dest/f"taxonomy_recovery_{rank}_{mode}.tsv")
        print(mode, "validated",len(df),"nonoverlapping on-target observations")
    p = pd.concat([load(f) for f in sorted((out/"profiles").glob("*_profile.tsv"))],ignore_index=True)
    domain = p.groupby(["sample","domain"],dropna=False)["coverage"].sum().unstack(fill_value=0)
    save(norm(domain)*100,dest/"domain_relative_abundance.tsv",True)
    save(p,dest/"combined_condensed_profiles.tsv")

def alpha(raw,out):
    dest = out/"analysis_filtered"
    records=[]
    for hits in [1,2,5,10]:
        for (sample,domain,gene),g in raw.loc[raw.num_hits>=hits].groupby(GROUP):
            x=g.coverage.to_numpy(float); x=x[x>0]
            if not len(x): continue
            p=x/x.sum(); h=-(p*np.log(p)).sum()
            records.append(dict(sample=sample,domain=domain,gene=gene,minimum_hits=hits,
                                marker_reads=g.num_hits.sum(),richness=len(x),Shannon=h,Effective_Shannon=np.exp(h)))
    d=pd.DataFrame(records)
    save(d,dest/"marker_alpha_diversity.tsv")
    s=d.groupby(["sample","domain","minimum_hits"]).agg(N_markers=("gene","nunique"),Median_Shannon=("Shannon","median"),Mean_Shannon=("Shannon","mean"),Median_richness=("richness","median"),Median_effective_Shannon=("Effective_Shannon","median"),Median_marker_reads=("marker_reads","median")).reset_index()
    save(s,dest/"sample_alpha_diversity.tsv")

def rarefy(bac,out):
    dest=out/"rarefaction"
    rng=np.random.default_rng(42)
    depth=bac.groupby(["gene","sample"]).num_hits.sum().unstack(fill_value=0)
    save(depth,dest/"marker_read_depth.tsv",True)
    minimum=depth.min(axis=1)
    primary=1000 if (minimum>=1000).sum()>=30 else 500
    markers=minimum[minimum>=primary].index
    assert len(markers)>=30
    maximum=int(minimum.loc[markers].min())
    depths=sorted(set(n for n in [50,100,250,500,750,1000,1500,2000,3000,primary,maximum] if n<=maximum))
    save(pd.DataFrame({"Minimum_reads":minimum,"Median_reads":depth.median(axis=1),"Maximum_reads":depth.max(axis=1),"Retained":minimum.index.isin(markers)}),dest/"marker_selection.tsv",True)
    rows=[]
    for (sample,gene),g in bac.loc[bac.gene.isin(markers)].groupby(["sample","gene"]):
        counts=g.num_hits.to_numpy(np.int64); counts=counts[counts>0]
        for n in depths:
            draws=rng.multivariate_hypergeometric(counts,n,size=100)
            for h in [1,2,5]:
                x=draws.copy(); x[x<h]=0
                total=x.sum(axis=1)
                p=np.divide(x,total[:,None],out=np.zeros(x.shape,float),where=total[:,None]>0)
                H=-(p*np.log(np.maximum(p,1e-300))).sum(axis=1)
                H[total==0]=np.nan
                rows.append(dict(Sample=sample,Marker=gene,Depth=n,Minimum_hits=h,Replicates=100,
                                 Richness_mean=(x>0).sum(axis=1).mean(),Shannon_mean=np.nanmean(H),
                                 Shannon_sd=np.nanstd(H),Effective_Shannon_mean=np.nanmean(np.exp(H))))
    r=pd.DataFrame(rows)
    save(r,dest/"marker_rarefaction.tsv")
    s=r.groupby(["Sample","Depth","Minimum_hits"]).agg(N_markers=("Marker","nunique"),Median_Shannon=("Shannon_mean","median"),Mean_Shannon=("Shannon_mean","mean"),Median_richness=("Richness_mean","median"),Median_effective_Shannon=("Effective_Shannon_mean","median")).reset_index()
    s["Substrate"]=s.Sample.map(substrate)
    save(s,dest/"sample_rarefaction.tsv")
    save(s.loc[s.Depth==primary],dest/"primary_rarefied_diversity.tsv")

def beta(bac,out):
    dest=out/"beta_diversity"
    save(pd.DataFrame({"Sample":SAMPLES,"Substrate":[substrate(s) for s in SAMPLES]}),dest/"sample_metadata.tsv")
    for h in [1,2,5]:
        distances=[]
        for gene,g in bac.loc[bac.num_hits>=h].groupby("gene"):
            x=g.pivot_table(index="sample",columns="sequence",values="coverage",aggfunc="sum",fill_value=0).reindex(SAMPLES,fill_value=0)
            x=norm(x).to_numpy(float)
            distances.append(np.stack([0.5*np.abs(x-v).sum(axis=1) for v in x]))
        stack=np.stack(distances)
        assert stack.shape==(35,19,19)
        for label,mat in [("mean",stack.mean(axis=0)),("median",np.median(stack,axis=0))]:
            save(pd.DataFrame(mat,index=SAMPLES,columns=SAMPLES),dest/f"bray_{label}_min{h}.tsv",True)
        np.save(dest/f"marker_distances_min{h}.npy",stack)

def mag_inputs(inputs):
    ab=pd.read_csv(inputs/"mag_abundance.tsv",sep="\t",index_col=0)
    tax=pd.read_csv(inputs/"mag_taxonomy.tsv",sep="\t",index_col=0)
    assert ab.shape==(19,565) and set(ab.columns)==set(tax.index)
    assert np.allclose(ab.sum(axis=1),1)
    tax=tax.loc[tax.Domain.eq("Bacteria")].copy()
    assert len(tax)==559
    return norm(ab[tax.index]),tax

def compare(out,inputs):
    ab,tax=mag_inputs(inputs)
    dest=out/"mag_comparison"
    results=[]; unmatched=[]
    for target in ["phylum","order","genus"]:
        allowed=RANKS[:[r[0] for r in RANKS].index(target)+1]
        prefix=allowed[-1][1]
        def mag_label(row):
            for name,p,col in reversed(allowed):
                v=row[col]
                if pd.notna(v) and str(v).strip():
                    label=f"{p}__{v}"
                    return label if name==target else "Unresolved@"+label
            return "Unresolved@d__Bacteria"
        def singlem_label(lineage):
            parts=[x.strip() for x in lineage.split(";")]
            target_parts=[x for x in parts if x.startswith(prefix+"__") and len(x)>3]
            resolved=[x for x in parts if "__" in x and len(x)>3]
            return target_parts[0] if target_parts else "Unresolved@"+(resolved[-1] if resolved else "d__Bacteria")
        mag=ab.T.groupby(tax.apply(mag_label,axis=1)).sum().T
        df=pd.read_csv(out/"analysis_filtered"/f"taxonomy_recovery_{target}_imperfect.tsv",sep="\t")
        df=df.loc[df.domain.eq("d__Bacteria")].copy()
        df["Taxon"]=df.lineage.map(singlem_label)
        all_sm=df.pivot_table(index="sample",columns="Taxon",values="Community_pct",aggfunc="sum",fill_value=0)
        linked=df.pivot_table(index="sample",columns="Taxon",values="MAG_linked",aggfunc="sum",fill_value=0)
        assert np.allclose(all_sm.sum(axis=1),100)
        assert set(mag.index)==set(all_sm.index)==set(linked.index)
        profiles={"SingleM_all":norm(all_sm),"SingleM_MAG_linked":norm(linked),"Original_MAGs":norm(mag)}
        taxa=sorted(set().union(*(set(x.columns) for x in profiles.values())))
        for name,x in profiles.items():
            profiles[name]=x.reindex(index=SAMPLES,columns=taxa,fill_value=0)
            save(profiles[name],dest/f"{target}_{name}.tsv",True)
        for s in SAMPLES:
            for a,b in combinations(profiles,2):
                x,y=profiles[a].loc[s],profiles[b].loc[s]
                results.append(dict(Rank=target,Sample=s,Profile_A=a,Profile_B=b,Bray_Curtis=np.abs(x-y).sum()/(x+y).sum()))
        for name,other in [("SingleM_all",set(mag.columns)),("Original_MAGs",set(all_sm.columns))]:
            for t in taxa:
                x=profiles[name][t]
                if t not in other and x.sum()>0:
                    unmatched.append(dict(Rank=target,Profile=name,Taxon=t,Mean_abundance_pct=100*x.mean(),Maximum_abundance_pct=100*x.max()))
    save(pd.DataFrame(results),dest/"profile_dissimilarity.tsv")
    save(pd.DataFrame(unmatched),dest/"unmatched_taxonomic_labels.tsv")

def diagnose(bac,out,inputs):
    ab,tax=mag_inputs(inputs)
    dest=out/"mag_comparison"
    uwma=set(tax.index[tax.Order.eq("UWMA-0217")])
    raw=bac.loc[bac["sample"].isin(["EGGPM35","EGGPM36","PM29"])].copy()
    raw["Assigned_order"]=raw.taxonomy.str.extract(r"(o__[^;]+)")[0].fillna("Unresolved")
    raw["marker_fraction"]=raw.coverage/raw.groupby(["sample","gene"]).coverage.transform("sum")
    rows=[]; tables=[]
    for mode in ["exact","imperfect"]:
        matched=load(dest/"targeted_appraisal"/f"{mode}_binned.tsv")
        matched=matched.loc[matched.domain.eq("d__Bacteria")].copy()
        matched["MAGs"]=matched.found_in.map(lambda x:set(re.findall(r"MAG\d+",str(x))))
        assert matched.MAGs.map(bool).all()
        matches=matched.groupby(KEY).MAGs.agg(lambda x:set().union(*x)).reset_index()
        df=raw.merge(matches,on=KEY,how="left",validate="one_to_one")
        df["MAGs"]=df.MAGs.map(lambda x:x if isinstance(x,set) else set())
        df["Any_MAG"]=df.MAGs.map(bool)
        df["MAG124_only"]=df.MAGs.map(lambda x:x=={"MAG124"})
        df["MAG124_possible"]=df.MAGs.map(lambda x:"MAG124" in x)
        df["UWMA_only"]=df.MAGs.map(lambda x:bool(x) and x.issubset(uwma))
        df["UWMA_possible"]=df.MAGs.map(lambda x:bool(x & uwma))
        df["UWMA_taxonomy"]=df.Assigned_order.eq("o__UWMA-0217")
        df["MAG124_other_order"]=df.MAG124_only & ~df.UWMA_taxonomy & df.Assigned_order.ne("Unresolved")
        df["MAG124_unresolved"]=df.MAG124_only & df.Assigned_order.eq("Unresolved")
        metrics=["Any_MAG","MAG124_only","MAG124_possible","UWMA_only","UWMA_possible","UWMA_taxonomy","MAG124_other_order","MAG124_unresolved"]
        for k in metrics: df[k]=df[k].astype(float)*df.marker_fraction
        marker=df.groupby(["sample","gene"])[metrics].sum().reset_index()
        marker[metrics]*=100
        marker["Mode"]=mode; tables.append(marker)
        for s,g in marker.groupby("sample"):
            v=g[metrics].mean(); linked=v.Any_MAG
            row=dict(Mode=mode,Sample=s,Original_MAG124_pct=100*ab.loc[s,"MAG124"],Original_UWMA_pct=100*ab.loc[s,sorted(uwma)].sum(),SingleM_UWMA_taxonomy_pct=v.UWMA_taxonomy,Total_MAG_linked_pct=linked,MAG124_other_order_pct=v.MAG124_other_order,MAG124_unresolved_pct=v.MAG124_unresolved,MAG124_median_marker_pct=g.MAG124_only.median(),MAG124_marker_Q25=g.MAG124_only.quantile(.25),MAG124_marker_Q75=g.MAG124_only.quantile(.75))
            for k in ["MAG124_only","MAG124_possible","UWMA_only","UWMA_possible"]:
                row[k+"_all_pct"]=v[k]
                row[k+"_linked_pct"]=100*v[k]/linked
            rows.append(row)
    save(pd.DataFrame(rows),dest/"UWMA0217_abundance_diagnostic.tsv")
    save(pd.concat(tables),dest/"UWMA0217_marker_abundance.tsv")

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--mag-inputs",type=Path,help="Directory generated from S3/S4 by script 00")
    p.add_argument("--stage",choices=["all","recovery","alpha","rarefaction","beta","compare","diagnose"],default="all")
    a=p.parse_args()
    raw=load(a.out/"analysis_qc/on_target_otu.tsv")
    raw=raw.loc[raw.domain.isin(["d__Bacteria","d__Archaea"])].copy()
    assert not raw.duplicated(KEY).any() and sorted(raw["sample"].unique())==SAMPLES
    assert (raw.coverage>=0).all() and (raw.num_hits>=0).all()
    assert np.all(raw.num_hits==raw.num_hits.astype(int))
    bac=raw.loc[raw.domain.eq("d__Bacteria")].copy()
    assert bac.gene.nunique()==35 and (bac.groupby("sample").gene.nunique()==35).all()
    for name,fn,data in [("recovery",recovery,raw),("alpha",alpha,raw),("rarefaction",rarefy,bac),("beta",beta,bac)]:
        if a.stage in ["all",name]: fn(data,a.out)
    if a.stage in ["compare","diagnose"] and not a.mag_inputs: p.error("--mag-inputs required for this stage")
    if a.mag_inputs and a.stage in ["all","compare"]: compare(a.out,a.mag_inputs)
    if a.mag_inputs and a.stage in ["all","diagnose"]: diagnose(bac,a.out,a.mag_inputs)

if __name__=="__main__":
    main()
