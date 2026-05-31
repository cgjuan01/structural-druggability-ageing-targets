"""
Step 2 - Structure-based druggability scoring with fpocket.

For each gene with an experimental structure, the first associated PDB entry is
downloaded and fpocket is run to detect ligand-binding pockets. The highest
fpocket Druggability Score across detected pockets is retained (0-1; >= 0.5 is
conventionally considered druggable).

Resumable: results are written line-by-line to a JSONL file, and genes already
present are skipped on re-run, so an interrupted session can be resumed by simply
re-running the script.

LIMITATION (stated honestly): a single representative structure is used per
target. A given crystal form may not contain the biologically relevant pocket
(e.g. a catalytic site absent from a particular construct), so individual scores
are hypotheses; the aggregate distribution is robust to this, individual calls
are not. Proper per-target analysis would curate multiple structures.
"""
import urllib.request
import json
import os
import subprocess
import re
import sys


def get_first_pdb_id(uniprot):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot}.json"
    with urllib.request.urlopen(url, timeout=30) as r:
        d = json.load(r)
    for x in d.get("uniProtKBCrossReferences", []):
        if x["database"] == "PDB":
            return x["id"]
    return None


def run_fpocket(pdb_path):
    subprocess.run(["fpocket", "-f", pdb_path], capture_output=True,
                   text=True, timeout=300)
    base = pdb_path.replace(".pdb", "")
    info = f"{base}_out/{os.path.basename(base)}_info.txt"
    if not os.path.exists(info):
        return None, None
    text = open(info).read()
    best_d, best_v = None, None
    for p in text.split("Pocket ")[1:]:
        dm = re.search(r"Druggability Score\s*:\s*([\d.]+)", p)
        vm = re.search(r"Volume\s*:\s*([\d.]+)", p)
        if dm:
            d = float(dm.group(1))
            if best_d is None or d > best_d:
                best_d, best_v = d, (float(vm.group(1)) if vm else 0)
    return best_d, best_v


def main(resolved_file="resolved.json", results_file="fpocket_results.jsonl"):
    resolved = json.load(open(resolved_file))
    crystals = [r for r in resolved if r["n_pdb"] > 0]
    os.makedirs("pdb", exist_ok=True)

    done = set()
    if os.path.exists(results_file):
        for line in open(results_file):
            try:
                done.add(json.loads(line)["gene"])
            except Exception:
                pass
    print(f"{len(crystals)} crystal-structure genes | {len(done)} already done")

    out = open(results_file, "a")
    for i, r in enumerate(crystals):
        g = r["gene"]
        if g in done:
            continue
        try:
            pdb_id = get_first_pdb_id(r["uniprot"])
            if not pdb_id:
                out.write(json.dumps({"gene": g, "uniprot": r["uniprot"],
                                      "drug_score": None, "note": "no_pdb_id"}) + "\n")
                out.flush()
                continue
            dest = f"pdb/{g}_{pdb_id}.pdb"
            if not os.path.exists(dest):
                urllib.request.urlretrieve(
                    f"https://files.rcsb.org/download/{pdb_id}.pdb", dest)
            d, v = run_fpocket(dest)
            out.write(json.dumps({"gene": g, "uniprot": r["uniprot"],
                                  "pdb": pdb_id, "drug_score": d, "volume": v}) + "\n")
            out.flush()
            os.system(f"rm -rf pdb/{g}_{pdb_id}_out {dest}")
        except Exception as e:
            out.write(json.dumps({"gene": g, "uniprot": r["uniprot"],
                                  "drug_score": None, "note": str(e)[:80]}) + "\n")
            out.flush()
        if (i + 1) % 25 == 0:
            print(f"  ...{i+1}/{len(crystals)}")
    out.close()
    print("done (re-run to resume if interrupted)")


if __name__ == "__main__":
    main(*sys.argv[1:])
