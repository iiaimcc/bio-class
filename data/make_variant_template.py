#!/usr/bin/env python3
"""Build the per-variant functional-annotation table from the MAF-derived status file.

Every non-silent variant in mutation_status.csv gets one row. Truncating variants
(nonsense, frameshift) are filled in automatically: a destroyed protein cannot be
switched on, so the call is loss_of_function and the evidence comes from the MAF
itself. Everything else - missense, in-frame indels, splice-site, start-codon,
nonstop - is left as needs_review, because the amino-acid change alone says
nothing about what it does to the protein.

Splice-site, start-codon and nonstop variants are deliberately NOT auto-called:
splice-site effects depend on which transcript and frame is used, a lost start
codon does not necessarily abolish the protein, and a nonstop variant extends
the C-terminus rather than truncating it.

Usage:  python3 make_variant_template.py [data_dir] [out_file]
"""

import csv
import os
import sys

AUTO_LOSS = {
    "Nonsense_Mutation": "stop gained, protein truncated",
    "Frame_Shift_Del": "frameshift, protein truncated",
    "Frame_Shift_Ins": "frameshift, protein truncated",
}
NEEDS_REVIEW = {
    "Missense_Mutation": "amino-acid substitution; effect not implied by the MAF",
    "In_Frame_Del": "in-frame deletion; effect not implied by the MAF",
    "In_Frame_Ins": "in-frame insertion; effect not implied by the MAF",
    "Splice_Site": "splice effect depends on transcript and frame; check it",
    "Translation_Start_Site": "lost start codon does not always abolish the protein",
    "Nonstop_Mutation": "stop lost, C-terminus extended; not a truncation",
}
FIELDS = ["symbol", "sample_id", "hgvsp", "variant_class",
          "functional_consequence", "evidence_level", "evidence_source"]


def main():
    data_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
    out_file = sys.argv[2] if len(sys.argv) > 2 else "variant_calls.csv"

    src = os.path.join(data_dir, "mutation_status.csv")
    if not os.path.exists(src):
        sys.exit(f"Missing {src} - run data/fetch_tcga_luad.py or extract the fallback copy")

    rows = []
    unknown = set()
    with open(src) as fh:
        for rec in csv.DictReader(fh):
            variant_class = rec["variant_class"]
            if variant_class in AUTO_LOSS:
                call, level, source = "loss_of_function", "high", f"truncating variant ({variant_class})"
            elif variant_class in NEEDS_REVIEW:
                call, level, source = "needs_review", "", NEEDS_REVIEW[variant_class]
            else:
                unknown.add(variant_class)
                call, level, source = "needs_review", "", "variant class not handled by the template"
            rows.append({
                "symbol": rec["gene"],
                "sample_id": rec["sample_id"],
                "hgvsp": rec["hgvsp"],
                "variant_class": variant_class,
                "functional_consequence": call,
                "evidence_level": level,
                "evidence_source": source,
            })

    if unknown:
        print("Note: unhandled variant classes left as needs_review:", ", ".join(sorted(unknown)))

    rows.sort(key=lambda r: (r["functional_consequence"] != "needs_review", r["symbol"], r["hgvsp"], r["sample_id"]))
    with open(out_file, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    auto = sum(1 for r in rows if r["functional_consequence"] != "needs_review")
    print(f"Wrote {out_file}: {len(rows)} variant rows, {auto} auto-called "
          f"(truncating), {len(rows) - auto} left as needs_review")
    print("Fill the needs_review rows with data/lookup_variant.py, then run 03_target_shortlist.R")


if __name__ == "__main__":
    main()
