#!/usr/bin/env python3
"""
Merges GISAID metadata XLS fragments for H1N1 and H3N2.
Output: H1N1_metadata/metadata.csv and H3N2_metadata/metadata.csv
Duplicates removed based on Isolate_Id.
"""

import os
import glob
import pandas as pd

BASE_DIR = os.path.dirname(os.path.abspath(__file__))

for subtype in ["H1N1", "H3N2"]:
    folder = os.path.join(BASE_DIR, f"{subtype}_metadata")
    files  = sorted(glob.glob(os.path.join(folder, "*.xls")))
    out    = os.path.join(folder, "metadata.csv")

    print(f"\n[{subtype}] Merging {len(files)} files...")

    frames = []
    for f in files:
        df = pd.read_excel(f)
        print(f"  {os.path.basename(f):45s}  {len(df):>7,} rows")
        frames.append(df)

    merged = pd.concat(frames, ignore_index=True)
    before = len(merged)
    merged.drop_duplicates(subset="Isolate_Id", keep="first", inplace=True)
    after  = len(merged)

    merged.to_csv(out, index=False)

    print(f"  {'─'*55}")
    print(f"  Total before dedup : {before:>7,}")
    print(f"  Duplicates removed : {before - after:>7,}")
    print(f"  Final rows         : {after:>7,}")
    print(f"  Saved → {os.path.relpath(out, BASE_DIR)}")

print("\nDone.")
