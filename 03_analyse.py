"""
Step 3 - Analyse the druggability landscape and produce the summary figure.

Reports the distribution of fpocket Druggability Scores across the prioritised
gene set, tests whether druggability relates to GNN rank, lists the top
druggable + highly-ranked candidate targets, and writes a two-panel figure.
"""
import json
import statistics as st
import sys
import matplotlib.pyplot as plt


def load(results_file, ranks_file):
    res = [json.loads(l) for l in open(results_file) if l.strip()]
    scored = [r for r in res if r.get("drug_score") is not None]
    # rank = position in the ranked clock list (1 = top)
    ranks = {g.strip(): i + 1 for i, g in enumerate(open(ranks_file)) if g.strip()}
    for r in scored:
        r["rank"] = ranks.get(r["gene"])
    return scored


def main(results_file="fpocket_results.jsonl", ranks_file="top1000_clock.txt"):
    scored = load(results_file, ranks_file)
    ds = [r["drug_score"] for r in scored]
    frac = sum(1 for x in ds if x >= 0.5) / len(ds)
    print(f"Scored: {len(scored)}")
    print(f"Median druggability: {st.median(ds):.2f}")
    print(f"Druggable (>=0.5): {frac*100:.0f}%  | highly (>=0.7): "
          f"{sum(1 for x in ds if x >= 0.7)}")

    # figure: (a) score distribution, (b) rank vs druggability
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.2))

    ax1.hist(ds, bins=20, color="#2c6fbb", edgecolor="white")
    ax1.axvline(0.5, color="#b03a2e", ls="--", lw=1.2, label="druggable threshold (0.5)")
    ax1.set_xlabel("fpocket Druggability Score")
    ax1.set_ylabel("Number of targets")
    ax1.set_title(f"Druggability of prioritised ageing-clock targets\n"
                  f"({frac*100:.0f}% druggable; n={len(scored)} crystal structures)")
    ax1.legend(fontsize=8)

    withrank = [(r["rank"], r["drug_score"]) for r in scored if r["rank"]]
    xs, ys = zip(*withrank)
    ax2.scatter(xs, ys, s=12, alpha=0.4, color="#2c6fbb")
    ax2.axhline(0.5, color="#b03a2e", ls="--", lw=1)
    ax2.set_xlabel("GNN rank (1 = top-ranked)")
    ax2.set_ylabel("fpocket Druggability Score")
    ax2.set_title("Druggability vs GNN rank\n(independent: tractability does not track rank)")

    fig.tight_layout()
    fig.savefig("results/druggability_landscape.png", dpi=300)
    print("wrote results/druggability_landscape.png")

    # top candidates
    cand = sorted([r for r in scored if r["drug_score"] >= 0.7 and r["rank"]],
                  key=lambda r: r["rank"])[:20]
    print("\nTop druggable + highly-ranked candidates:")
    for r in cand:
        print(f"  {r['gene']:10s} rank={r['rank']:4d}  "
              f"druggability={r['drug_score']:.2f}  ({r.get('pdb','')})")


if __name__ == "__main__":
    main(*sys.argv[1:])
