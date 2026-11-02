#!/usr/bin/env python3
"""Fetch the TCGA-LUAD teaching subset from the GDC API.

The sample set is pinned in sample_manifest.csv (file_id + md5 + case/sample
barcodes), so every run downloads the same 30 files and rebuilds the same
matrix. Mutation calls are pinned the same way in mutation_manifest.csv.

Usage:
    python3 fetch_tcga_luad.py                     # download + build, output next to this script
    python3 fetch_tcga_luad.py --cache ./gdc_cache  # reuse previously downloaded files
    python3 fetch_tcga_luad.py --workers 4 --out .
    python3 fetch_tcga_luad.py --skip-mutation      # expression only
    python3 fetch_tcga_luad.py --verify-only        # check existing outputs, no download

Requires only the Python standard library.
"""

import argparse
import csv
import gzip
import hashlib
import json
import os
import sys
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor

GDC_DATA = "https://api.gdc.cancer.gov/data/{file_id}"
HERE = os.path.dirname(os.path.abspath(__file__))

# Only these genes are extracted from the MAF files; the list covers the
# actionable drivers discussed in the course. Fusions and copy-number events
# are not in a MAF, so absence of a row is not evidence of wild type.
DRIVER_GENES = [
    "EGFR", "KRAS", "TP53", "STK11", "KEAP1", "BRAF", "PIK3CA",
    "ERBB2", "ALK", "RET", "ROS1", "MET", "NF1", "SMARCA4",
]
NONCODING = {
    "Silent", "Intron", "3'UTR", "5'UTR", "3'Flank", "5'Flank",
    "IGR", "RNA", "Targeted_Region",
}


def md5sum(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256sum(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def download(file_id, expected_md5, dest, attempts=5):
    """Download one GDC file and verify its md5. Returns the local path."""
    if os.path.exists(dest) and md5sum(dest) == expected_md5:
        return dest
    url = GDC_DATA.format(file_id=file_id)
    for i in range(1, attempts + 1):
        try:
            tmp = dest + ".part"
            with urllib.request.urlopen(url, timeout=600) as r, open(tmp, "wb") as f:
                while True:
                    chunk = r.read(1 << 20)
                    if not chunk:
                        break
                    f.write(chunk)
            if md5sum(tmp) != expected_md5:
                raise ValueError("md5 mismatch")
            os.replace(tmp, dest)
            return dest
        except Exception as exc:  # network hiccup or partial file
            print(f"  retry {i}/{attempts} {file_id}: {exc}", file=sys.stderr)
            if os.path.exists(dest + ".part"):
                os.remove(dest + ".part")
            time.sleep(3 * i)
    raise RuntimeError(f"failed to download {file_id}")


def read_manifest(path, required):
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    missing = required - set(rows[0].keys()) if rows else required
    if missing:
        raise SystemExit(f"{path} is missing columns: {sorted(missing)}")
    return rows


def build_expression(rows, cache, out_dir, workers):
    """Merge the per-sample STAR count files into one gene x sample matrix."""
    def fetch(row):
        dest = os.path.join(cache, row["file_id"] + ".tsv")
        download(row["file_id"], row["md5"], dest)
        return row, dest

    with ThreadPoolExecutor(max_workers=workers) as pool:
        fetched = list(pool.map(fetch, rows))

    genes, annotation, matrix = [], {}, {}
    for row, path in fetched:
        sample = row["sample_id"]
        with open(path, encoding="utf-8") as f:
            for line in f:
                if line.startswith("#") or line.startswith("gene_id"):
                    continue
                fields = line.rstrip("\n").split("\t")
                gene_id = fields[0]
                if gene_id.startswith("N_"):        # STAR QC rows
                    continue
                if gene_id not in annotation:
                    annotation[gene_id] = (fields[1], fields[2])
                    genes.append(gene_id)
                    matrix[gene_id] = []
                matrix[gene_id].append(int(fields[3]))  # unstranded counts

    samples = [row["sample_id"] for row, _ in fetched]
    n = len(samples)
    for gene_id, values in matrix.items():
        if len(values) != n:
            raise SystemExit(f"{gene_id}: {len(values)} values for {n} samples")
    # A gene missing from one file would silently shift columns; the check
    # above catches that, so every row can be written as-is.

    expr_path = os.path.join(out_dir, "expression_matrix.csv")
    with open(expr_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["gene_id"] + samples)
        for gene_id in genes:
            w.writerow([gene_id] + matrix[gene_id])

    ann_path = os.path.join(out_dir, "gene_annotation.csv")
    with open(ann_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["gene_id", "symbol", "gene_type"])
        for gene_id in genes:
            w.writerow([gene_id, *annotation[gene_id]])

    meta_path = os.path.join(out_dir, "sample_meta.csv")
    with open(meta_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["sample_id", "case_id", "condition", "sample_type", "tss"])
        for row in rows:
            w.writerow([row["sample_id"], row["case_id"], row["condition"],
                        row["sample_type"], row["sample_id"].split("-")[1]])

    return expr_path, ann_path, meta_path, len(genes), len(samples)


def build_mutations(rows, cache, out_dir, workers):
    """Extract non-silent variants in DRIVER_GENES from the pinned MAF files."""
    def fetch(row):
        dest = os.path.join(cache, row["file_id"] + ".maf.gz")
        download(row["file_id"], row["md5"], dest)
        return row, dest

    with ThreadPoolExecutor(max_workers=workers) as pool:
        fetched = list(pool.map(fetch, rows))

    variants = []
    for row, path in fetched:
        with gzip.open(path, "rt", encoding="utf-8") as f:
            reader = csv.DictReader((l for l in f if not l.startswith("#")),
                                    delimiter="\t")
            for rec in reader:
                gene = rec.get("Hugo_Symbol", "")
                if gene not in DRIVER_GENES:
                    continue
                if rec.get("Variant_Classification", "") in NONCODING:
                    continue
                variants.append({
                    "sample_id": row["sample_id"],
                    "case_id": row["case_id"],
                    "gene": gene,
                    "variant_class": rec["Variant_Classification"],
                    "hgvsp": rec.get("HGVSp_Short", ""),
                    "chromosome": rec.get("Chromosome", ""),
                    "start": rec.get("Start_Position", ""),
                })

    variants.sort(key=lambda v: (v["sample_id"], v["gene"], v["hgvsp"]))
    long_path = os.path.join(out_dir, "mutation_status.csv")
    with open(long_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["sample_id", "case_id", "gene",
                                         "variant_class", "hgvsp",
                                         "chromosome", "start"])
        w.writeheader()
        w.writerows(variants)

    samples = [row["sample_id"] for row in rows]
    hits = {(v["sample_id"], v["gene"]) for v in variants}
    matrix_path = os.path.join(out_dir, "mutation_matrix.csv")
    with open(matrix_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["sample_id"] + DRIVER_GENES)
        for sample in samples:
            w.writerow([sample] + [int((sample, g) in hits) for g in DRIVER_GENES])

    return long_path, matrix_path, len(variants)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", default=os.path.join(HERE, "sample_manifest.csv"))
    ap.add_argument("--mutation-manifest",
                    default=os.path.join(HERE, "mutation_manifest.csv"))
    ap.add_argument("--out", default=HERE)
    ap.add_argument("--cache", default=os.path.join(HERE, "gdc_cache"))
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--skip-mutation", action="store_true")
    ap.add_argument("--verify-only", action="store_true",
                    help="re-hash the built CSV files in --out and compare with "
                         "the values printed by a rebuild; downloads nothing")
    args = ap.parse_args()

    os.makedirs(args.out, exist_ok=True)
    os.makedirs(args.cache, exist_ok=True)

    if args.verify_only:
        # The published values live in data/README.md; this mode only reports
        # what is on disk so a teacher can diff it against those values.
        missing = []
        for name in ("expression_matrix.csv", "gene_annotation.csv", "sample_meta.csv",
                     "mutation_status.csv", "mutation_matrix.csv"):
            path = os.path.join(args.out, name)
            if os.path.exists(path):
                print(f"  {sha256sum(path)}  {name}")
            else:
                missing.append(name)
        if missing:
            raise SystemExit("Missing output files: " + ", ".join(missing))
        return

    rows = read_manifest(args.manifest,
                         {"file_id", "file_name", "md5", "sample_id", "case_id",
                          "condition", "sample_type"})
    print(f"Manifest: {len(rows)} expression files "
          f"({sum(r['condition'] == 'Tumor' for r in rows)} tumor / "
          f"{sum(r['condition'] == 'Normal' for r in rows)} normal)")

    expr, ann, meta, n_genes, n_samples = build_expression(
        rows, args.cache, args.out, args.workers)
    print(f"Wrote {expr}: {n_genes} genes x {n_samples} samples")
    print(f"Wrote {ann}")
    print(f"Wrote {meta}")

    if not args.skip_mutation:
        mrows = read_manifest(args.mutation_manifest,
                              {"file_id", "file_name", "md5", "sample_id", "case_id"})
        long_path, matrix_path, n_variants = build_mutations(
            mrows, args.cache, args.out, args.workers)
        print(f"Wrote {long_path}: {n_variants} non-silent variants in {len(DRIVER_GENES)} genes")
        print(f"Wrote {matrix_path}")

    print("\nSHA256 (record these in data/README.md after a rebuild):")
    for name in ("expression_matrix.csv", "gene_annotation.csv", "sample_meta.csv",
                 "mutation_status.csv", "mutation_matrix.csv"):
        path = os.path.join(args.out, name)
        if os.path.exists(path):
            print(f"  {sha256sum(path)}  {name}")


if __name__ == "__main__":
    main()
