
import csv, gzip, json, os, re, collections

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATASET = os.path.join(ROOT, "dataset")
ANALYSIS = os.path.join(ROOT, "analysis")
VS = os.path.join(os.path.dirname(HERE), "variant_summary.txt.gz")
PTC = re.compile(r'p\.[A-Za-z]{3}(\d+)(Ter|\*)', re.I)
FS = re.compile(r'p\.[A-Za-z]{3}(\d+)[A-Za-z]{0,3}fs', re.I)

ns = collections.Counter(); fs = collections.Counter()
with gzip.open(VS, "rt", encoding="utf-8", errors="replace") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        if row.get("Assembly") != "GRCh38":
            continue
        sig = (row.get("ClinicalSignificance") or "").lower()
        if "pathogenic" not in sig or "conflicting" in sig:
            continue
        name = row.get("Name") or ""
        if FS.search(name):
            cls = "fs"
        elif PTC.search(name):
            cls = "ns"
        else:
            continue
        for s in (row.get("GeneSymbol") or "").split(";"):
            s = s.strip()
            if not s or s == "-":
                continue
            (fs if cls == "fs" else ns)[s] += 1

genes = sorted(set(ns) | set(fs))
rows = []
for g in genes:
    a, b = fs[g], ns[g]
    t = a + b
    if t:
        rows.append((g, b, a, t, a / t))

def tally(minboth, mintot):
    return [r for r in rows if min(r[1], r[2]) >= minboth and r[3] >= mintot]

print(f"genes with at least one pathogenic truncating variant: {len(rows)}")
for mb, mt in ((1, 10), (3, 10), (5, 20), (5, 30), (10, 40)):
    s = tally(mb, mt)
    print(f"  both classes >= {mb} and total >= {mt}: {len(s)} genes")

pool = [r for r in rows if r[3] >= 20]
dspp = next((r for r in rows if r[0] == "DSPP"), None)
print(f"\ngenes with >= 20 pathogenic truncating variants: {len(pool)}")
if dspp:
    print(f"DSPP: nonsense {dspp[1]}, frameshift {dspp[2]}, share {dspp[4]:.3f}")
    hi = [r for r in pool if r[4] >= dspp[4]]
    lo = [r for r in pool if r[4] <= 0.25]
    print(f"  genes at or above DSPP's frameshift share: {len(hi)} "
          f"({round(100*len(hi)/len(pool),1)} percent of the pool)")
    print(f"  genes at or below a 0.25 frameshift share: {len(lo)} "
          f"({round(100*len(lo)/len(pool),1)} percent)")
    print("  top frameshift-skewed, with counts:")
    for r in sorted(hi, key=lambda r: -r[3])[:12]:
        print(f"    {r[0]:12s} ns={r[1]:4d} fs={r[2]:4d} share={r[4]:.3f}")
    print("  top nonsense-skewed, with counts:")
    for r in sorted(lo, key=lambda r: -r[3])[:12]:
        print(f"    {r[0]:12s} ns={r[1]:4d} fs={r[2]:4d} share={r[4]:.3f}")
shares = sorted(r[4] for r in pool)
q = lambda p: round(shares[int(p * (len(shares) - 1))], 3)
print(f"\nframeshift share across the pool: min {q(0)}, 10th {q(.1)}, median {q(.5)}, "
      f"90th {q(.9)}, max {q(1)}")
json.dump({"n_genes_with_ptv": len(rows),
           "qualifying": {f"both_ge{mb}_total_ge{mt}": len(tally(mb, mt))
                          for mb, mt in ((1,10),(3,10),(5,20),(5,30),(10,40))},
           "pool_ge20": len(pool),
           "dspp": {"nonsense": dspp[1], "frameshift": dspp[2], "share": round(dspp[4],3)} if dspp else None,
           "n_at_or_above_dspp_share": len([r for r in pool if dspp and r[4] >= dspp[4]]),
           "share_quantiles": {"min": q(0), "p10": q(.1), "median": q(.5), "p90": q(.9), "max": q(1)}},
          open(os.path.join(ANALYSIS, "derived_genomewide_power.json"), "w"), indent=1, sort_keys=True)
