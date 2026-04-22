#!/usr/bin/env python3
"""
Runs Nextclade alignment on clustered.fasta in every epoch folder.

Uses --output-all to produce all Nextclade output files per epoch:
  nextclade.aligned.fasta         — sequences aligned to epoch reference
  nextclade.json                  — full per-sequence results (JSON)
  nextclade.ndjson                — same, newline-delimited (streamable)
  nextclade.tsv                   — mutations + QC (tab-separated)
  nextclade.csv                   — same (comma-separated)
  nextclade.gff                   — per-sequence GFF3 annotations
  nextclade.tbl                   — feature table
  nextclade.cds_translation.*.fasta — amino acid translations per CDS
    (only in codon-aware mode, i.e. when --input-annotation is provided)

Note on .nwk (Newick tree):
  Nextclade tree output requires --input-tree, which must match the
  reference sequence exactly. Since we use epoch-specific references
  (not the dataset root), the dataset tree is incompatible and causes
  a hard crash. Build Newick trees downstream from nextclade.aligned.fasta
  using IQ-TREE or FastTree.

Per-epoch strategy:
  - Uses the epoch's own reference.fasta as --input-ref
  - Tries codon-aware mode (with dataset annotation + pathogen.json)
    when epoch reference length is within 50 bp of the dataset reference
  - Falls back to basic mode (no annotation) otherwise

The reference sequence is prepended to clustered.fasta so it appears
first in nextclade.aligned.fasta.

Idempotent: skips epochs where nextclade.aligned.fasta already exists.
"""

import os
import subprocess
import tempfile
import time

BASE_DIR    = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR  = os.path.join(BASE_DIR, "split_output")
DATASETS    = os.path.join(BASE_DIR, "nextclade_datasets")
NEXTCLADE   = os.path.join(BASE_DIR, "nextclade")

SUBTYPE_MAP   = {"H1N1": "H1N1", "H3N2": "H3N2"}
LEN_TOLERANCE = 50


def seq_len(fasta_path: str) -> int:
    length = 0
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                if length:
                    break
            else:
                length += len(line.strip())
    return length


def count_seqs(fasta_path: str) -> int:
    return sum(1 for l in open(fasta_path) if l.startswith(">"))


def run_nextclade(cmd: list) -> bool:
    try:
        subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        return True
    except (subprocess.CalledProcessError, OSError):
        return False


total = skipped = succeeded = failed = no_ref = 0

for subtype in sorted(os.listdir(OUTPUT_DIR)):
    sdir = os.path.join(OUTPUT_DIR, subtype)
    if not os.path.isdir(sdir) or subtype not in SUBTYPE_MAP:
        continue

    for segment in sorted(os.listdir(sdir)):
        segdir = os.path.join(sdir, segment)
        if not os.path.isdir(segdir):
            continue

        ds_dir        = os.path.join(DATASETS, SUBTYPE_MAP[subtype], segment)
        annotation    = os.path.join(ds_dir, "genome_annotation.gff3")
        pathogen_json = os.path.join(ds_dir, "pathogen.json")
        ds_ref        = os.path.join(ds_dir, "reference.fasta")

        has_dataset = os.path.isdir(ds_dir)
        ds_ref_len  = seq_len(ds_ref) if has_dataset and os.path.exists(ds_ref) else None

        if not has_dataset:
            print(f"\n[WARN] No dataset for {subtype}/{segment} — skipping segment")
            continue

        print(f"\n[{subtype}] {segment}")

        for epoch in sorted(os.listdir(segdir)):
            epoch_dir = os.path.join(segdir, epoch)
            if not os.path.isdir(epoch_dir):
                continue

            clustered  = os.path.join(epoch_dir, "clustered.fasta")
            reference  = os.path.join(epoch_dir, "reference.fasta")
            aligned_out = os.path.join(epoch_dir, "nextclade.aligned.fasta")

            if not os.path.exists(clustered):
                continue

            total += 1

            if not os.path.exists(reference):
                no_ref += 1
                print(f"  [SKIP] {epoch}  — no reference.fasta")
                continue

            if os.path.exists(aligned_out):
                skipped += 1
                print(f"  [--]   {epoch}  (already done)")
                continue

            n_in = count_seqs(clustered)
            epoch_ref_len = seq_len(reference)

            use_annotation = (
                has_dataset
                and ds_ref_len is not None
                and abs(epoch_ref_len - ds_ref_len) <= LEN_TOLERANCE
            )

            # Build combined input: reference first, then clustered sequences
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".fasta", delete=False, dir=epoch_dir
            ) as tmp:
                tmp_path = tmp.name
                with open(reference) as rf:
                    tmp.write(rf.read())
                with open(clustered) as cf:
                    tmp.write(cf.read())

            base_cmd = [
                NEXTCLADE, "run",
                "--input-ref",       reference,
                "--output-all",      epoch_dir + "/",
                "--output-basename", "nextclade",
            ]

            anno_cmd = base_cmd + [
                "--input-annotation",    annotation,
                "--input-pathogen-json", pathogen_json,
                tmp_path,
            ] if use_annotation else None

            no_anno_cmd = base_cmd + [tmp_path]

            t0 = time.time()
            ok = False
            mode = ""

            if anno_cmd:
                ok = run_nextclade(anno_cmd)
                if ok:
                    mode = "codon-aware"
                else:
                    # Clean partial outputs and retry without annotation
                    for fname in os.listdir(epoch_dir):
                        if fname.startswith("nextclade."):
                            os.remove(os.path.join(epoch_dir, fname))

            if not ok:
                ok = run_nextclade(no_anno_cmd)
                mode = "basic" if ok else "—"

            elapsed = time.time() - t0

            if os.path.exists(tmp_path):
                os.remove(tmp_path)

            if ok and os.path.exists(aligned_out):
                n_out = count_seqs(aligned_out)
                succeeded += 1
                print(f"  [OK]   {epoch}  {n_in:,} → {n_out:,} seqs  [{mode}]  ({elapsed:.0f}s)")
            else:
                failed += 1
                for fname in os.listdir(epoch_dir):
                    if fname.startswith("nextclade."):
                        os.remove(os.path.join(epoch_dir, fname))
                print(f"  [FAIL] {epoch}  ({elapsed:.0f}s)")

print(f"\n{'='*65}")
print(f" Done.")
print(f"   Aligned  : {succeeded}")
print(f"   Skipped  : {skipped}  (already done)")
print(f"   No ref   : {no_ref}  (New Caledonia — excluded)")
print(f"   Failed   : {failed}")
print(f"   Total    : {total}")
print(f"{'='*65}")
