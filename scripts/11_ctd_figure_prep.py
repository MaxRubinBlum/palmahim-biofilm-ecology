#!/usr/bin/env python3
"""Prepare Supplementary Fig. 1 from a processed CTD/ROV timeline table.

Data required: timestamp, dissolved_oxygen, and optional event/event_type columns.
This script does not reconstruct Sea-Bird calibration from raw .hex files; that
instrument-specific conversion is described in Supplementary Methods.
"""
import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--timeline",required=True); ap.add_argument("--out",required=True); a=ap.parse_args(); sep="\t" if a.timeline.endswith((".tsv",".txt")) else ","; d=pd.read_csv(a.timeline,sep=sep); d["timestamp"]=pd.to_datetime(d.timestamp)
    fig,ax=plt.subplots(figsize=(9,4)); ax.plot(d.timestamp,d.dissolved_oxygen,lw=1); ax.set_ylabel("Dissolved oxygen"); ax.set_xlabel("Time")
    if "event" in d.columns:
        ev=d[d.event.notna()]
        for _,r in ev.iterrows(): ax.axvline(r.timestamp,alpha=.2,lw=.7)
    fig.tight_layout(); Path(a.out).parent.mkdir(parents=True,exist_ok=True); fig.savefig(a.out,dpi=300)
if __name__=="__main__": main()
