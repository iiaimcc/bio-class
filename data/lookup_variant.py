#!/usr/bin/env python3
"""Look a variant up in ClinVar and CIViC and print what those databases say.

The script only fetches: it converts the short protein change from the MAF
(e.g. p.L858R) into the HGVS form ClinVar indexes (p.Leu858Arg), prints every
matching ClinVar record with its classification and review status, and prints
the accepted CIViC evidence for the same variant. Deciding what the evidence
means - activating, loss of function, or uncertain - stays with the reader.

Usage:  python3 lookup_variant.py GENE VARIANT
        python3 lookup_variant.py EGFR p.L858R
        python3 lookup_variant.py RET p.T350N
"""

import json
import sys
import time
import urllib.parse
import urllib.request

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
CIVIC = "https://civicdb.org/api/graphql"
AA = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln",
    "E": "Glu", "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys",
    "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp",
    "Y": "Tyr", "V": "Val",
}


def get_json(url, tries=5):
    last = None
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=90) as resp:
                return json.load(resp)
        except Exception as exc:  # NCBI returns 500 under load; retry with backoff
            last = exc
            time.sleep(3 * (attempt + 1))
    raise RuntimeError(f"request failed after {tries} attempts: {last}")


def to_hgvs3(short):
    """p.L858R -> p.Leu858Arg ; p.E746_A750del -> p.Glu746_Ala750del"""
    body = short[2:] if short.startswith("p.") else short
    if "_" in body:
        left, rest = body.split("_")
        right = rest.replace("del", "").replace("ins", "").replace(">", "")
        try:
            return f"p.{AA[left[0]]}{left[1:]}_{AA[right[0]]}{right[1:]}del"
        except (KeyError, IndexError):
            return short
    try:
        return f"p.{AA[body[0]]}{body[1:-1]}{AA[body[-1]]}"
    except (KeyError, IndexError):
        return short


def clinvar(gene, hgvs3):
    term = f'{gene}[gene] AND "{hgvs3}"'
    ids = get_json(EUTILS + "esearch.fcgi?db=clinvar&retmode=json&retmax=5&term="
                   + urllib.parse.quote(term))["esearchresult"]["idlist"]
    if not ids:
        return []
    result = get_json(EUTILS + "esummary.fcgi?db=clinvar&retmode=json&id=" + ",".join(ids))["result"]
    out = []
    for key in result.get("uids", []):
        rec = result[key]
        germ = rec.get("germline_classification") or {}
        conditions = ", ".join(t.get("trait_name", "") for t in (germ.get("trait_set") or []))
        out.append({
            "accession": rec.get("accession"),
            "title": rec.get("title"),
            "classification": germ.get("description"),
            "review_status": germ.get("review_status"),
            "conditions": conditions,
            "url": f"https://www.ncbi.nlm.nih.gov/clinvar/variation/{key}/",
        })
    return out


def civic_query(query):
    req = urllib.request.Request(
        CIVIC, data=json.dumps({"query": query}).encode(),
        headers={"Content-Type": "application/json", "User-Agent": "bio-class-course/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=90) as resp:
            return json.load(resp)
    except Exception as exc:
        return {"error": str(exc)}


def civic(gene, short):
    name = short[2:] if short.startswith("p.") else short
    found = civic_query(
        '{browseVariants(featureName:"%s", variantName:"%s", first:5){nodes{id name}}}' % (gene, name))
    nodes = ((found.get("data") or {}).get("browseVariants") or {}).get("nodes") or []
    out = []
    for node in nodes:
        evidence = civic_query(
            "{evidenceItems(variantId:%d, status:ACCEPTED, first:25)"
            "{nodes{evidenceType significance evidenceLevel evidenceDirection}}}"
            % node["id"])
        items = ((evidence.get("data") or {}).get("evidenceItems") or {}).get("nodes") or []
        seen = []
        for item in items:
            line = (f"{item['evidenceType']} / {item['significance']} "
                    f"(level {item['evidenceLevel']}, {item['evidenceDirection']})")
            if line not in seen:
                seen.append(line)
        out.append({
            "id": node["id"], "name": node["name"], "evidence": seen,
            "url": f"https://civicdb.org/variants/{node['id']}/summary",
        })
    return out


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    gene, short = sys.argv[1], sys.argv[2]
    hgvs3 = to_hgvs3(short)
    print(f"# {gene} {short}  (ClinVar index form: {hgvs3})")

    print("\n## ClinVar")
    try:
        records = clinvar(gene, hgvs3)
    except Exception as exc:
        records = None
        print(f"  lookup failed ({exc}); try again in a minute")
    if records is not None and not records:
        print("  no record for this protein change")
    for rec in records:
        print(f"  {rec['accession']} | {rec['title']}")
        print(f"    classification: {rec['classification']} ({rec['review_status']})")
        if rec["conditions"]:
            print(f"    conditions: {rec['conditions']}")
        print(f"    {rec['url']}")

    print("\n## CIViC (accepted evidence only)")
    variants = civic(gene, short)
    if not variants:
        print("  no variant with this name")
    for var in variants:
        print(f"  {var['name']} (id {var['id']}) {var['url']}")
        for line in var["evidence"]:
            print(f"    {line}")

    print("\nDecide the call yourself (gain_of_function / loss_of_function / uncertain),")
    print("record the evidence level and the source, and note what you checked.")


if __name__ == "__main__":
    main()
