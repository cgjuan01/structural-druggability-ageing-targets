"""
Step 1 - Resolve gene symbols to UniProt accessions and structure availability.

For each gene we record its reviewed UniProt accession and the number of
experimental PDB structures. Genes with experimental structures are used for the
primary (trustworthy) fpocket analysis; AlphaFold-only genes are handled as a
separate, lower-confidence exploratory set.
"""
import urllib.request
import urllib.parse
import json
import time
import sys


def resolve(gene):
    q = f"gene_exact:{gene} AND organism_id:9606 AND reviewed:true"
    url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(
        {"query": q, "fields": "accession,xref_pdb", "format": "json", "size": 1})
    with urllib.request.urlopen(url, timeout=30) as r:
        d = json.load(r)
    if not d.get("results"):
        return None, 0
    res = d["results"][0]
    acc = res["primaryAccession"]
    npdb = sum(1 for x in res.get("uniProtKBCrossReferences", [])
               if x["database"] == "PDB")
    return acc, npdb


def main(gene_file, out_file):
    genes = [l.strip() for l in open(gene_file) if l.strip()]
    resolved = []
    for i, g in enumerate(genes):
        try:
            acc, npdb = resolve(g)
            resolved.append({"gene": g, "uniprot": acc, "n_pdb": npdb})
        except Exception:
            resolved.append({"gene": g, "uniprot": None, "n_pdb": 0})
        if (i + 1) % 100 == 0:
            print(f"  ...{i+1}/{len(genes)}")
        time.sleep(0.05)
    json.dump(resolved, open(out_file, "w"), indent=2)
    n_ok = sum(1 for r in resolved if r["uniprot"])
    n_xtal = sum(1 for r in resolved if r["n_pdb"] > 0)
    print(f"Resolved {n_ok}/{len(genes)} | crystal {n_xtal} | "
          f"AlphaFold-only {n_ok - n_xtal} | unresolved {len(genes) - n_ok}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "top1000_clock.txt",
         sys.argv[2] if len(sys.argv) > 2 else "resolved.json")
