#!/usr/bin/env python3
"""
Splits GISAID Influenza A FASTA files by timeline intervals.

Output structure:
  split_output/
    H1N1/  HA/  1986_A_Singapore_6_1986/
                  sequences.fasta
                  reference.fasta
           NA/  ...
    H3N2/  ...

Reference sequence: the canonical strain listed in the timeline for that epoch,
searched (case-insensitive, partial match) inside the FASTA headers.
Sequences that fall before the first epoch or have unparseable dates are skipped.
"""

import os
import re
from collections import defaultdict
from Bio import SeqIO

BASE_DIR   = os.path.dirname(os.path.abspath(__file__))
TIMELINE   = os.path.join(BASE_DIR, "influenza_timeline_latest.md")
OUTPUT_DIR = os.path.join(BASE_DIR, "split_output")


# ---------------------------------------------------------------------------
# Timeline parsing
# ---------------------------------------------------------------------------

def parse_timeline(md_file: str) -> dict:
    data = {"H1N1": [], "H3N2": []}
    current = None

    with open(md_file) as f:
        for line in f:
            if "H1N1" in line and "Timeline" in line:
                current = "H1N1"
            elif "H3N2" in line and "Timeline" in line:
                current = "H3N2"
            elif line.startswith("|") and current:
                parts = [p.strip() for p in line.split("|")]
                if len(parts) > 3:
                    try:
                        year = int(parts[1])
                        strain = parts[3]
                        data[current].append({"year": year, "strain": strain})
                    except ValueError:
                        pass

    for subtype in data:
        data[subtype] = sorted(data[subtype], key=lambda x: x["year"])

    return data


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def epoch_index(date_str: str, intervals: list):
    """Return which epoch index a collection date belongs to, or None."""
    try:
        year = int(date_str.split("-")[0])
    except Exception:
        return None

    for i in range(len(intervals) - 1):
        if intervals[i]["year"] <= year < intervals[i + 1]["year"]:
            return i

    if year >= intervals[-1]["year"]:
        return len(intervals) - 1

    return None  # before the first epoch


def epoch_name(interval: dict) -> str:
    safe = re.sub(r"[\s/]+", "_", interval["strain"])
    return f"{interval['year']}_{safe}"


def _norm(s: str) -> str:
    """Lowercase and treat underscores as spaces."""
    return s.strip().lower().replace("_", " ")


def _year_match(y_ref: str, y_header: str) -> bool:
    """True if two year strings refer to the same year (handles 2-digit vs 4-digit)."""
    if y_ref == y_header:
        return True
    # 2-digit reference year: "95" matches "1995"
    if len(y_ref) == 2 and y_header.endswith(y_ref) and len(y_header) == 4:
        return True
    # 2-digit header year (unusual but possible)
    if len(y_header) == 2 and y_ref.endswith(y_header) and len(y_ref) == 4:
        return True
    return False


def ref_matches(header: str, strain: str) -> bool:
    """
    Strict component-wise match of a timeline strain name against a GISAID header.

    Rules:
    - Works on the strain-name field only (before first '|').
    - Normalises underscores ↔ spaces.
    - Requires type (A) and location to match exactly.
    - Requires isolate number to match exactly (prevents "10" matching "109").
    - Year comparison handles 2-digit vs 4-digit (e.g. "95" == "1995").
    - 3-part WHO names (A/Location/Year with no isolate number) match if
      type + location match and the year appears as the LAST slash-component
      of the header strain field.
    """
    parts_s = [_norm(p) for p in strain.split("/")]
    parts_h = [_norm(p) for p in header.split("|")[0].split("/")]

    if len(parts_s) < 2 or len(parts_h) < 2:
        return False

    # Type must match (always "a")
    if parts_s[0] != parts_h[0]:
        return False

    # Location must match exactly (after normalisation)
    if parts_s[1] != parts_h[1]:
        return False

    # ---- 4-part strain: A/Location/Isolate/Year ----
    if len(parts_s) == 4:
        isolate_s = parts_s[2]
        year_s    = parts_s[3]

        if len(parts_h) >= 3:
            # Isolate must match exactly (no prefix matching)
            if parts_h[2] != isolate_s:
                return False
            # Year check (if present in header)
            if len(parts_h) >= 4:
                return _year_match(year_s, parts_h[3])
            else:
                return True   # header has no year field — accept on type+loc+isolate

    # ---- 3-part WHO name: A/Location/Year (no isolate) ----
    if len(parts_s) == 3:
        year_s = parts_s[2]
        # Match if the last slash-component of the header is the same year
        if len(parts_h) >= 2:
            return _year_match(year_s, parts_h[-1])

    return False


# ---------------------------------------------------------------------------
# Core processing
# ---------------------------------------------------------------------------

def process_segment(subtype: str, segment: str, fasta_path: str,
                    intervals: list, output_base: str):
    """
    Single-pass through the FASTA file:
      • Streams each record into the correct epoch file.
      • Simultaneously checks EVERY record against ALL epoch reference strains
        (reference strains are often collected before their named epoch begins,
        so we must search the whole file, not just each epoch's own records).
    """
    n_epochs = len(intervals)
    epoch_dirs  = {}
    seq_handles = {}
    seq_counts  = defaultdict(int)
    ref_records = {}   # epoch_idx -> best-matching SeqRecord
    ref_found   = set()

    # Pre-build (lowercased) strain strings for fast per-record checking
    epoch_strains = [iv["strain"] for iv in intervals]

    # Create output dirs and open per-epoch sequence file handles
    for idx, interval in enumerate(intervals):
        ename = epoch_name(interval)
        edir  = os.path.join(output_base, subtype, segment, ename)
        os.makedirs(edir, exist_ok=True)
        epoch_dirs[idx]  = (edir, interval)
        seq_handles[idx] = open(os.path.join(edir, "sequences.fasta"), "w")

    total = 0
    try:
        for record in SeqIO.parse(fasta_path, "fasta"):
            total += 1
            fields = record.description.split("|")
            if len(fields) < 3:
                continue

            # --- epoch assignment (by collection date) ---
            epoch_idx = epoch_index(fields[2].strip(), intervals)
            if epoch_idx is not None:
                SeqIO.write([record], seq_handles[epoch_idx], "fasta")
                seq_counts[epoch_idx] += 1

            # --- reference matching: check against ALL epochs ---
            for eidx, strain in enumerate(epoch_strains):
                if eidx not in ref_found and ref_matches(record.description, strain):
                    ref_records[eidx] = record
                    ref_found.add(eidx)
                    # stop early once all refs are found
                    if len(ref_found) == n_epochs:
                        break

    finally:
        for h in seq_handles.values():
            h.close()

    # Write reference files and print summary
    print(f"    Total in file: {total:,}")
    for idx in range(n_epochs):
        edir, interval = epoch_dirs[idx]
        ename = epoch_name(interval)
        count = seq_counts.get(idx, 0)

        if idx in ref_records:
            ref_path = os.path.join(edir, "reference.fasta")
            SeqIO.write([ref_records[idx]], ref_path, "fasta")
            status = f"ref OK  ({ref_records[idx].id.split('|')[0]})"
        else:
            status = f"ref NOT FOUND for '{interval['strain']}'"

        print(f"      {ename:50s}  {count:6,} seqs  |  {status}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("=" * 65)
    print(" Influenza Segment-Wise Data Splitter")
    print("=" * 65)

    timeline = parse_timeline(TIMELINE)

    for subtype, intervals in timeline.items():
        print(f"\n{'='*65}")
        print(f" {subtype}  —  {len(intervals)} epochs")
        print(f"{'='*65}")

        subtype_dir = os.path.join(BASE_DIR, subtype)
        if not os.path.isdir(subtype_dir):
            print(f"  [!] {subtype_dir} not found — skipping.")
            continue

        fasta_files = sorted(f for f in os.listdir(subtype_dir) if f.endswith(".fasta"))
        for fname in fasta_files:
            segment = fname.replace(".fasta", "")
            fasta_path = os.path.join(subtype_dir, fname)
            print(f"\n  [{segment}]")
            process_segment(subtype, segment, fasta_path, intervals, OUTPUT_DIR)

    print(f"\n{'='*65}")
    print(f" Done. Output in: {OUTPUT_DIR}")
    print(f"{'='*65}")


if __name__ == "__main__":
    main()
