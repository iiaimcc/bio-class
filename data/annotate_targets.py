#!/usr/bin/env python3
"""Resolve the driver genes to UniProt entries and record where each protein sits.

The third question of the course-1 target funnel is "can the large molecule we
design in course 3 reach it?". That is a property of the protein, so it is
fetched from UniProt rather than maintained by hand: an extracellular
topological domain (or secretion) is what makes an antibody-like binder
possible, while a cytoplasmic or nuclear protein is out of reach for it.

Output: <data_dir>/target_context.csv

Usage:  python3 annotate_targets.py [data_dir]
"""

import csv
import json
import os
import sys
import urllib.parse
import urllib.request

UNIPROT = "https://rest.uniprot.org/uniprotkb/search"
FIELDS = "accession,gene_primary,protein_name,cc_subcellular_location,ft_topo_dom"
EXTRACELLULAR_DOMAINS = {"Extracellular", "Extracellular space", "Lumenal"}


def driver_genes(data_dir):
    with open(os.path.join(data_dir, "mutation_matrix.csv")) as fh:
        header = next(csv.reader(fh))
    return [g for g in header if g != "sample_id"]


def fetch(genes):
    query = " OR ".join(f"(gene_exact:{g})" for g in genes)
    url = (
        f"{UNIPROT}?query={urllib.parse.quote(query, safe='=&:+()')}"
        f" AND organism_id:9606 AND reviewed:true"
        f"&fields={FIELDS}&format=json&size=100"
    )
    url = url.replace(" ", "%20")
    with urllib.request.urlopen(url, timeout=120) as resp:
        release = resp.headers.get("x-uniprot-release", "unknown")
        return json.load(resp)["results"], release


def topology_of(entry):
    for feature in entry.get("features", []):
        if feature["type"] == "Topological domain":
            desc = feature["description"]
            if desc in EXTRACELLULAR_DOMAINS:
                return desc
    return ""


def locations_of(entry):
    out = []
    for comment in entry.get("comments", []):
        if comment["commentType"] != "SUBCELLULAR LOCATION":
            continue
        for loc in comment.get("subcellularLocations", []):
            value = loc["location"]["value"]
            if "topology" in loc:
                value += f" ({loc['topology']['value']})"
            out.append(value)
    return sorted(set(out))


def main():
    data_dir = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
    genes = driver_genes(data_dir)
    results, release = fetch(genes)

    by_gene = {}
    for entry in results:
        primary = entry["genes"][0]["geneName"]["value"]
        # gene_exact also matches synonyms (e.g. MET is a synonym of SLTM)
        if primary in genes and primary not in by_gene:
            by_gene[primary] = entry

    missing = [g for g in genes if g not in by_gene]
    if missing:
        sys.exit(f"No reviewed human UniProt entry for: {', '.join(missing)}")

    rows = []
    for gene in genes:
        entry = by_gene[gene]
        topo = topology_of(entry)
        locations = locations_of(entry)
        reachable = bool(topo) or any("Secreted" in loc for loc in locations)
        rows.append({
            "symbol": gene,
            "uniprot_accession": entry["primaryAccession"],
            "protein_name": entry["proteinDescription"]["recommendedName"]["fullName"]["value"],
            "subcellular_location": "; ".join(locations) or "not annotated in UniProt",
            "extracellular_domain": topo or "none",
            "has_extracellular_domain": "yes" if reachable else "no",
            "uniprot_release": release,
            "source": f"https://rest.uniprot.org/uniprotkb/{entry['primaryAccession']}/entry",
        })

    out_path = os.path.join(data_dir, "target_context.csv")
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    surface = [r["symbol"] for r in rows if r["has_extracellular_domain"] == "yes"]
    print(f"Wrote {out_path}: {len(rows)} genes, UniProt release {release}")
    print(f"Reachable by a large molecule (extracellular domain): {', '.join(surface)}")


if __name__ == "__main__":
    main()
