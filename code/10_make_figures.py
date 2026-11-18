#!/usr/bin/env python3
import csv, json, os, re, sys, collections
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D
import numpy as np
from scipy.stats import fisher_exact

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATASET = os.path.join(ROOT, "dataset")
ANALYSIS = os.path.join(ROOT, "analysis")
DATA = ANALYSIS
FIGS = os.path.join(ROOT, "figures")

S1, S2, S3, S4 = "#7a2d8f", "#167c51", "#c0196b", "#6f6a16"
GREY, WHITE = "#b9b8b3", "#ffffff"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8a8984"
TITLE_SIZE, TITLE_PAD = 9.5, 7

plt.rcParams.update({
    "svg.hashsalt": "dental-truncation-frame",
    "figure.facecolor": WHITE, "axes.facecolor": WHITE, "savefig.facecolor": WHITE,
    "font.family": "DejaVu Sans", "font.size": 8.5,
    "axes.edgecolor": MUTED, "axes.linewidth": 0.7,
    "axes.labelcolor": INK, "axes.titlesize": TITLE_SIZE, "axes.titleweight": "bold",
    "xtick.color": INK2, "ytick.color": INK2,
    "xtick.labelsize": 8, "ytick.labelsize": 8,
    "legend.frameon": False, "legend.fontsize": 8,
    "grid.color": "#e9e8e4", "grid.linewidth": 0.5,
    "lines.linewidth": 2.0, "lines.markersize": 4.5,
})

def load(name):
    return json.load(open(os.path.join(DATA, name)))

V2, PH = load("derived_v2.json"), load("derived_phenotypes.json")
FR, LS = load("derived_frame.json"), load("derived_labels_submission.json")
FS2, EX = load("derived_fsd2.json"), load("derived_fsd_exclusion.json")
GW = load("derived_genomewide_power.json")

with open(os.path.join(DATASET, "variant_table.csv")) as fh:
    VT = list(csv.DictReader(fh))
for r in VT:
    for k in ("codon", "protein_aa", "fsd", "allele_count", "homozygote_count",
              "nmd_escape_from_codon", "in_escape_window"):
        r[k] = int(r[k]) if r[k] not in ("", None) else None
    for k in ("rel_pos", "nmd_escape_rel"):
        r[k] = float(r[k]) if r[k] not in ("", None) else None

LIT = {r["gene"]: r["literature_symbol"] for r in VT}
POOL = sorted((g for g, x in PH["dual_phenotype_genes"].items() if x["positioned_ptv"] >= 10),
              key=lambda g: -(LS["genes"][g]["n_dental"] + LS["genes"][g]["n_non_dental"]))

def sel(**kw):
    out = VT
    for k, v in kw.items():
        out = [r for r in out if r[k] == v] if not isinstance(v, (list, tuple, set)) \
              else [r for r in out if r[k] in v]
    return out

def tidy(ax, ylab=None, xlab=None, title=None, grid="y"):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    if grid:
        ax.set_axisbelow(True); ax.grid(axis=grid, which="major")
    if ylab: ax.set_ylabel(ylab)
    if xlab: ax.set_xlabel(xlab)
    if title: ax.set_title(title, loc="left", color=INK, pad=TITLE_PAD)

def letter(fig, ax, ch, dx=-0.085, dy=1.035):
    ax.text(dx, dy, ch, transform=ax.transAxes, fontsize=11, weight="bold",
            color=INK, ha="left", va="bottom")

def save(fig, name):
    os.makedirs(FIGS, exist_ok=True)
    for ext, md, kw in (("png", {}, {"dpi": 400}), ("svg", {"Date": None}, {})):
        fig.savefig(os.path.join(FIGS, f"{name}.{ext}"), bbox_inches="tight",
                    facecolor=WHITE, metadata=md, **kw)
    plt.close(fig)
    print("wrote", name + ".png and .svg")

def stars(p):
    if p is None: return ""
    return "p < 0.001" if p < 1e-3 else f"p = {p:.3f}"

def figure1():
    fig = plt.figure(figsize=(7.4, 9.1))
    gs = fig.add_gridspec(3, 1, height_ratios=[1.45, 0.95, 1.15], hspace=0.18)
    axt = fig.add_subplot(gs[0])
    tooth_panel(axt)
    axt.text(0.0, 1.0, "a", fontsize=11, weight="bold", color=INK, va="top",
             ha="left", transform=axt.transAxes)
    axt.text(0.038, 1.0, "The panel by tissue", fontsize=TITLE_SIZE, weight="bold",
             color=INK, va="top", ha="left", transform=axt.transAxes)
    axa = fig.add_subplot(gs[1]); axb = fig.add_subplot(gs[2])

    axa.set_xlim(0, 10); axa.set_ylim(-1.15, 5.15); axa.axis("off")
    axa.text(0.0, 4.86, "Three ways a truncating variant can act", fontsize=TITLE_SIZE,
             weight="bold", color=INK, va="bottom")
    axa.text(-0.55, 4.86, "b", fontsize=11, weight="bold", color=INK, va="bottom")
    X0, X1, XE = 1.7, 8.6, 7.2
    for k, (lab, col, note) in enumerate([
            ("Uniform loss of function", GREY, "disease anywhere"),
            ("Position restricted", S1, "disease only in the escape window"),
            ("Frame dependent", S2, "frameshift only")]):
        y = 3.30 - k * 1.46
        axa.add_patch(Rectangle((X0, y - 0.15), X1 - X0, 0.30, facecolor="#f3f2ef",
                                edgecolor=MUTED, linewidth=0.7))
        axa.add_patch(Rectangle((XE, y - 0.15), X1 - XE, 0.30, facecolor="#ece2f2",
                                edgecolor=MUTED, linewidth=0.7))
        axa.text(X0, y + 0.62, lab, fontsize=8.8, weight="bold", color=col, va="bottom")
        axa.text(X0, y + 0.29, note, fontsize=8, color=INK2, va="bottom")
        tri = {0: [2.3, 3.3, 4.3, 5.4, 6.4, 7.9], 1: [7.5, 7.9, 8.3],
               2: [3.6, 5.0, 6.2, 7.8]}[k]
        for x in tri:
            axa.plot([x], [y], marker="v", color=col, markersize=6,
                     markeredgecolor=WHITE, markeredgewidth=0.6, clip_on=False)
        if k == 2:
            for x in (2.5, 4.3, 8.3):
                axa.plot([x], [y], marker="o", color=GREY, markersize=5.2,
                         markeredgecolor=WHITE, markeredgewidth=0.6, clip_on=False)
    axa.annotate("escape window", xy=((XE + X1) / 2, 3.48), xytext=((XE + X1) / 2, 4.26),
                 fontsize=7.8, color=INK2, ha="center",
                 arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=0.7,
                                 shrinkA=0, shrinkB=2))
    axa.plot([X0 + 0.15], [-0.86], marker="v", color=INK2, markersize=6, clip_on=False)
    axa.text(X0 + 0.35, -0.86, "pathogenic", fontsize=8, color=INK2, va="center")
    axa.plot([X0 + 2.6], [-0.86], marker="o", color=GREY, markersize=5.2, clip_on=False)
    axa.text(X0 + 2.8, -0.86, "tolerated in gnomAD", fontsize=8, color=INK2, va="center")

    genes = sorted((g for g in LIT if V2["genes"][g]["plp_ptv"] > 0),
                   key=lambda g: -V2["genes"][g]["plp_ptv"])
    x = np.arange(len(genes))
    ptv = [V2["genes"][g]["plp_ptv"] for g in genes]
    vus = [V2["genes"][g]["vus_total"] for g in genes]
    axb.bar(x - 0.2, ptv, width=0.4, color=S1, label="pathogenic truncating")
    axb.bar(x + 0.2, vus, width=0.4, color=GREY, label="uncertain significance")
    axb.set_yscale("log")
    axb.set_xticks(x)
    axb.set_xticklabels([LIT[g] for g in genes], rotation=90, fontsize=7.2)
    for i, g in enumerate(genes):
        if g in POOL:
            axb.text(i, 0.55, "●", ha="center", va="center", fontsize=6.5, color=S2)
    tidy(axb, ylab="variants (log scale)",
         title="ClinVar entries per gene")
    letter(fig, axb, "c", dx=-0.062)
    axb.legend(handles=[Patch(facecolor=S1, label="pathogenic truncating"),
                        Patch(facecolor=GREY, label="uncertain significance"),
                        Line2D([], [], marker="o", ls="none", color=S2, markersize=4.5,
                               label="also causes a non-dental disease")],
               loc="upper right", ncol=1, bbox_to_anchor=(1.0, 1.03), fontsize=7.6)
    axb.set_ylim(0.4, max(vus) * 4.2)
    save(fig, "Figure1")

def figure2():
    fig, axes = plt.subplots(len(POOL), 1, figsize=(7.4, 1.05 * len(POOL)), sharex=True)
    for ax, g in zip(axes, POOL):
        rec = LS["genes"][g]
        esc = next((r["nmd_escape_rel"] for r in sel(gene=g) if r["nmd_escape_rel"]), None)
        if esc:
            ax.add_patch(Rectangle((esc, -0.6), 1.02 - esc, 2.2, facecolor="#f4edf8",
                                   edgecolor="none", zorder=0))
            ax.axvline(esc, color=MUTED, lw=0.8, ls=(0, (3, 2)), zorder=1)
        groups = [("dental", S2, 1.0), ("non_dental", S1, 0.45), ("population", GREY, -0.1)]
        for side, col, y in groups:
            src = "gnomad_plof" if side == "population" else "clinvar_pathogenic"
            xs = [r["rel_pos"] for r in VT if r["gene"] == g and r["source"] == src
                  and r["phenotype_side"] == side and r["rel_pos"] is not None]
            ax.plot(xs, [y] * len(xs), marker="|", ls="none", color=col,
                    markersize=9, markeredgewidth=1.1, alpha=0.85)
            ax.text(1.03, y, str(len(xs)), fontsize=7.4, color=col, va="center")
        ax.set_ylim(-0.6, 1.6); ax.set_yticks([])
        ax.set_xlim(0, 1.0)
        for sp in ("top", "right", "left"):
            ax.spines[sp].set_visible(False)
        if g != POOL[-1]:
            ax.spines["bottom"].set_visible(False)
            ax.tick_params(axis="x", length=0)
        nm = rec.get("literature_symbol") or g
        w = rec["windows"].get("nmd_50nt_rule") or {}
        p = w.get("p_dental_vs_non_dental")
        lab = f"{nm}" + (f"   {stars(p)}" if p is not None else "")
        ax.text(-0.075, 0.5, lab, transform=ax.transAxes, fontsize=8.4,
                weight="bold" if (p is not None and p < 0.05) else "normal",
                color=INK if (p is not None and p < 0.05) else INK2,
                ha="right", va="center")
    axes[-1].set_xlabel("position along the protein, as a fraction of its length")
    axes[0].set_title("Truncation position by reported phenotype", loc="left",
                      color=INK, pad=26)
    axes[0].legend(handles=[Line2D([], [], color=S2, lw=2, label="dental phenotype"),
                            Line2D([], [], color=S1, lw=2, label="other phenotype"),
                            Line2D([], [], color=GREY, lw=2, label="gnomAD")],
                   loc="lower left", bbox_to_anchor=(0.0, 1.04), ncol=3)
    save(fig, "Figure2")

def figure3():
    defs = ["nmd_50nt_rule", "last_exon", "last_10pct", "last_25pct"]
    labs = ["50 nt rule", "final exon", "last 10%", "last 25%"]
    fig = plt.figure(figsize=(7.4, 5.0))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.06, 1.06, 1.06], wspace=0.82)
    for j, (key, ttl) in enumerate([("p_dental_vs_non_dental", "Dental against other"),
                                    ("p_dental_vs_gnomad", "Dental against gnomAD")]):
        ax = fig.add_subplot(gs[j])
        M = np.full((len(POOL), len(defs)), np.nan)
        for i, g in enumerate(POOL):
            for k, d in enumerate(defs):
                v = (LS["genes"][g]["windows"].get(d) or {}).get(key)
                if v is not None:
                    M[i, k] = -np.log10(max(v, 1e-20))
        im = ax.imshow(M, cmap="Purples", vmin=0, vmax=6, aspect="auto")
        ax.set_xticks(range(len(defs))); ax.set_xticklabels(labs, rotation=45, ha="right")
        ax.set_yticks(range(len(POOL)))
        ax.set_yticklabels([LS["genes"][g].get("literature_symbol") or g for g in POOL], fontsize=8)
        for i in range(len(POOL)):
            for k in range(len(defs)):
                if not np.isnan(M[i, k]):
                    v = 10 ** (-M[i, k])
                    ax.text(k, i, "*" if v < 0.05 else "", ha="center", va="center",
                            fontsize=11, color=S2 if v < 0.05 else INK2, weight="bold")
        ax.set_title(ttl, loc="left", color=INK, pad=TITLE_PAD)
        ax.grid(False)
        letter(fig, ax, "ab"[j], dx=-0.42 if j == 0 else -0.3)
        if j == 1:
            cb = fig.colorbar(im, ax=ax, fraction=0.050, pad=0.06)
            cb.ax.set_title("$-\\log_{10}$ p", fontsize=7.6, color=INK2, pad=5)
            cb.ax.tick_params(labelsize=7.5)
            cb.outline.set_linewidth(0.6)
    ax = fig.add_subplot(gs[2])
    g = next(g for g in POOL if (LS["genes"][g].get("literature_symbol") or g) == "LAMB3")
    w = LS["genes"][g]["windows"]["nmd_50nt_rule"]
    cats = ["dental", "other", "gnomAD"]
    inside = [w["dental_in"], w["non_dental_in"], w["gnomad_in"]]
    outside = [w["dental_out"], w["non_dental_out"], w["gnomad_out"]]
    frac = [a / (a + b) * 100 for a, b in zip(inside, outside)]
    bars = ax.bar(cats, frac, color=[S2, S1, GREY], width=0.58)
    ax.tick_params(axis="x", labelsize=7.6)
    for b, a, t in zip(bars, inside, [a + b for a, b in zip(inside, outside)]):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.6, f"{a}/{t}",
                ha="center", fontsize=7.8, color=INK)
    tidy(ax, ylab="inside the escape window (%)",
         title="LAMB3")
    ax.set_ylim(0, max(frac) * 1.35)
    ax.text(0.98, 0.95, stars(w["p_dental_vs_non_dental"]), transform=ax.transAxes,
            fontsize=8, color=S2, weight="bold", ha="right", va="top")
    letter(fig, ax, "c", dx=-0.14)
    save(fig, "FigureS1")

def figure4():
    fig = plt.figure(figsize=(7.4, 3.9))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.1], wspace=0.34)
    ax = fig.add_subplot(gs[0])
    ax.plot([0, 1], [0, 1], color=MUTED, lw=0.8, ls=(0, (3, 2)), zorder=1)
    for g, r in FR["genes"].items():
        if not r.get("path_fs_share") or not r.get("gnomad_fs_share"):
            continue
        n = r["n_path_nonsense"] + r["n_path_frameshift"]
        if n < 5:
            continue
        q = r.get("T1_q")
        hit = q is not None and q < 0.05
        ax.scatter([r["gnomad_fs_share"]], [r["path_fs_share"]], s=16 + 2.2 * n,
                   color=S2 if hit else S1, alpha=0.9 if hit else 0.5,
                   edgecolor=WHITE, linewidth=0.7, zorder=3)
        nm = r.get("literature_symbol") or g
        if hit or nm in ("FAM83H", "LAMB3", "WDR72", "EDA"):
            ax.annotate(nm, (r["gnomad_fs_share"], r["path_fs_share"]),
                        textcoords="offset points", xytext=(7, -3), fontsize=7.6,
                        color=S2 if hit else INK2, weight="bold" if hit else "normal")
    tidy(ax, xlab="frameshift share in gnomAD", ylab="frameshift share, pathogenic",
         title="Which class causes disease", grid="both")
    ax.set_xlim(0, 1.02); ax.set_ylim(0, 1.05)
    letter(fig, ax, "a", dx=-0.17)

    ax = fig.add_subplot(gs[1])
    q = GW["share_quantiles"]
    ax.axvspan(q["p10"], q["p90"], color="#f1f0ec", zorder=0)
    for key, col, nm in (("dspp", S2, "DSPP"),):
        ax.axvline(GW["dspp"]["share"], color=col, lw=1.6, zorder=3)
        ax.text(GW["dspp"]["share"] - 0.02, 0.86, nm, color=col, fontsize=8.4,
                weight="bold", ha="right", transform=ax.get_xaxis_transform())
    f83 = next(r for g, r in FR["genes"].items()
               if (r.get("literature_symbol") or g) == "FAM83H")
    ax.axvline(f83["path_fs_share"], color=S1, lw=1.6, zorder=3)
    ax.text(f83["path_fs_share"] + 0.02, 0.86, "FAM83H", color=S1, fontsize=8.4,
            weight="bold", transform=ax.get_xaxis_transform())
    for v, lab in ((q["median"], "median"),):
        ax.axvline(v, color=MUTED, lw=0.9, ls=(0, (3, 2)), zorder=2)
        ax.text(v, 0.965, lab, color=INK2, fontsize=7.6, ha="center", va="top",
                transform=ax.get_xaxis_transform(),
                bbox=dict(facecolor=WHITE, edgecolor="none", pad=1.2))
    ax.set_xlim(0.15, 1.02); ax.set_yticks([])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.set_xlabel("frameshift share, pathogenic")
    ax.set_title("Genome-wide context", loc="left", color=INK, pad=TITLE_PAD)
    ax.text(0.5, 0.42, "10th to 90th centile", transform=ax.transAxes,
            fontsize=7.8, color=INK2, ha="center")
    letter(fig, ax, "b", dx=-0.06)
    save(fig, "Figure3")

def ecdf(xs):
    xs = np.sort(np.asarray(xs, dtype=float))
    return xs, np.arange(1, len(xs) + 1) / len(xs)

def figure5():
    fig = plt.figure(figsize=(7.4, 5.3))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.88], hspace=0.62, wspace=0.26)
    dspp = next(g for g in FS2["genes"]
                if (FS2["genes"][g].get("literature_symbol") or g) == "DSPP")

    def fsds(gene, source):
        return [r["fsd"] for r in VT if r["gene"] == gene and r["source"] == source
                and r["variant_class"] == "frameshift" and r["fsd"] is not None]

    allp = [r["fsd"] for r in VT if r["source"] == "clinvar_pathogenic"
            and r["variant_class"] == "frameshift" and r["fsd"] is not None]
    allq = [r["fsd"] for r in VT if r["source"] == "gnomad_plof"
            and r["variant_class"] == "frameshift" and r["fsd"] is not None]
    exp = [r["fsd"] for r in VT if r["source"] == "clinvar_pathogenic" and r["gene"] != dspp
           and r["variant_class"] == "frameshift" and r["fsd"] is not None]
    exq = [r["fsd"] for r in VT if r["source"] == "gnomad_plof" and r["gene"] != dspp
           and r["variant_class"] == "frameshift" and r["fsd"] is not None]

    ax = fig.add_subplot(gs[0, 0])
    curves = ((fsds(dspp, "clinvar_pathogenic"), S2, "pathogenic"),
              (fsds(dspp, "gnomad_plof"), GREY, "gnomAD"))
    for xs, col, lab in curves:
        x, y = ecdf([max(v, 1) for v in xs])
        ax.step(x, y, where="post", color=col, lw=2.0)
    ax.set_xscale("log")
    tidy(ax, xlab="novel residues before the new stop", ylab="cumulative fraction",
         title="DSPP", grid="both")
    ax.set_ylim(-0.03, 1.10)
    ax.text(780, 0.21, "pathogenic\n(n = %d)" % len(curves[0][0]), fontsize=7.6,
            color=S2, ha="right", va="center", weight="bold")
    ax.text(1.05, 0.82, "gnomAD\n(n = %d)" % len(curves[1][0]), fontsize=7.6,
            color=INK2, ha="left", va="center", weight="bold")
    r = FS2["genes"][dspp]
    ax.text(0.03, 0.99, f"median {r['median_fsd_pathogenic']} against "
                        f"{r['median_fsd_population']}\n{stars(r['mwu_p'])}",
            transform=ax.transAxes, fontsize=7.6, color=INK, va="top")
    letter(fig, ax, "a", dx=-0.17)

    ax = fig.add_subplot(gs[0, 1])
    for xs, col, ls in ((allp, S2, "-"), (exp, S2, (0, (3, 2))),
                        (allq, GREY, "-"), (exq, GREY, (0, (3, 2)))):
        x, y = ecdf([max(v, 1) for v in xs])
        ax.step(x, y, where="post", color=col, lw=1.8, ls=ls)
    ax.set_xscale("log")
    tidy(ax, xlab="novel residues before the new stop", title="Whole panel", grid="both")
    ax.set_ylim(-0.03, 1.10)
    ax.text(0.03, 0.96, f"all {stars(EX['pooled_all']['mwu_p'])}\n"
                        f"minus DSPP {stars(EX['pooled_excluding_DSPP']['mwu_p'])}",
            transform=ax.transAxes, fontsize=7.6, color=INK, va="top")
    ax.legend(handles=[Line2D([], [], color=S2, lw=1.8, label="pathogenic"),
                       Line2D([], [], color=GREY, lw=1.8, label="gnomAD"),
                       Line2D([], [], color=INK2, lw=1.2, ls=(0, (3, 2)),
                              label="dashed: minus DSPP")],
              loc="lower right", fontsize=7.0, handlelength=1.8,
              borderpad=0.2, labelspacing=0.25)
    letter(fig, ax, "b", dx=-0.13)

    ax = fig.add_subplot(gs[1, :])
    ts = [5, 10, 20, 40, 80, 160]
    for xs, col, ls, lab in ((fsds(dspp, "clinvar_pathogenic"), S2, "-", "DSPP, pathogenic"),
                             (fsds(dspp, "gnomad_plof"), GREY, "-", "DSPP, gnomAD"),
                             (exp, S2, (0, (3, 2)), "panel minus DSPP, pathogenic"),
                             (exq, GREY, (0, (3, 2)), "panel minus DSPP, gnomAD")):
        ys = [100 * sum(1 for v in xs if v >= t) / len(xs) for t in ts]
        ax.plot(ts, ys, color=col, ls=ls, lw=1.8, marker="o", markersize=3.4, label=lab)
    ax.set_xscale("log"); ax.set_xticks(ts); ax.set_xticklabels([str(t) for t in ts])
    ax.set_xlim(4.4, 230)
    tidy(ax, xlab="threshold, novel residues", ylab="at or above the threshold (%)",
         title="Threshold swept")
    ax.legend(loc="center right", fontsize=7.4, handlelength=2.2, borderpad=0.4,
              labelspacing=0.35)
    letter(fig, ax, "c", dx=-0.075, dy=1.11)
    save(fig, "Figure4")

def figure6_grid():
    genes = sorted((g for g in LIT if V2["genes"][g]["plp_ptv"] > 0),
                   key=lambda g: -V2["genes"][g]["plp_ptv"])
    tests = ["phenotype separation", "pathogenic against population",
             "composition, whole gene", "composition, in window", "neo-tail length"]
    fig, ax = plt.subplots(figsize=(7.4, 0.3 * len(genes) + 1.7))
    for i, g in enumerate(genes):
        ls = LS["genes"].get(g); fr = FR["genes"].get(g, {}); f2 = FS2["genes"].get(g, {})
        vals = []
        w = (ls or {}).get("windows", {}).get("nmd_50nt_rule", {}) if ls else {}
        vals.append(w.get("p_dental_vs_non_dental"))
        vals.append(w.get("p_dental_vs_gnomad"))
        vals.append(fr.get("T1_q"))
        vals.append(fr.get("T3_q"))
        vals.append(f2.get("mwu_p"))
        for k, v in enumerate(vals):
            if v is None:
                ax.plot([k], [i], marker="x", color="#d9d8d4", markersize=5,
                        markeredgewidth=1.0)
            elif v < 0.05:
                ax.plot([k], [i], marker="o", color=S2, markersize=8,
                        markeredgecolor=WHITE, markeredgewidth=0.7)
            else:
                ax.plot([k], [i], marker="o", color="#e7e6e2", markersize=6.5,
                        markeredgecolor=MUTED, markeredgewidth=0.5)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels([LIT[g] for g in genes], fontsize=7.6)
    ax.set_xticks(range(len(tests)))
    ax.set_xticklabels(tests, rotation=28, ha="right", fontsize=8)
    ax.set_xlim(-0.6, len(tests) - 0.4); ax.set_ylim(-0.8, len(genes) - 0.2)
    ax.invert_yaxis()
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.grid(axis="y", color="#f2f1ed", lw=0.5); ax.set_axisbelow(True)
    ax.set_title("Where each test reaches significance, by gene",
                 loc="left", color=INK, pad=28)
    ax.legend(handles=[Line2D([], [], marker="o", ls="none", markerfacecolor=S2,
                              markeredgecolor="none", label="p or q below 0.05", markersize=8),
                       Line2D([], [], marker="o", ls="none", markerfacecolor="#e7e6e2",
                              markeredgecolor=MUTED, label="tested, not significant", markersize=6.5),
                       Line2D([], [], marker="x", ls="none", color="#d9d8d4",
                              label="counts too low to test", markersize=5)],
              loc="lower left", bbox_to_anchor=(0.0, 1.005), ncol=3)
    save(fig, "FigureS2")


TISSUE = [
    ("Enamel", re.compile(r'amelogenesis|enamel', re.I), S1),
    ("Dentine", re.compile(r'dentinogenesis|dentin dysplasia|dental radicular', re.I), S2),
    ("Tooth number", re.compile(r'tooth agenesis|hypodontia|oligodontia|anodontia', re.I), S3),
    ("Periodontium", re.compile(r'periodont|tooth eruption', re.I), S4),
]

def tissue_groups():
    out = collections.OrderedDict((t[0], []) for t in TISSUE)
    for g, e in V2["genes"].items():
        if e["tier"] != "isolated_dental":
            continue
        conds = " ; ".join(e["isolated"])
        for name, pat, _ in TISSUE:
            if pat.search(conds):
                out[name].append(LIT.get(g, g))
                break
    for k in out:
        out[k] = sorted(out[k])
    return out

ENAMEL_F, DENTINE_F, PULP_F, SOCKET_F = "#efedf4", "#eae6dd", "#e8bfcd", "#dedad3"

CROWN_OUT = [(0.00, 0.66), (0.12, 0.86), (0.26, 0.97), (0.39, 0.87), (0.49, 0.66),
             (0.54, 0.40), (0.51, 0.16), (0.44, 0.00)]
ROOT_OUT = [(0.42, -0.23), (0.38, -0.62), (0.30, -0.99), (0.23, -1.22), (0.165, -1.11),
            (0.125, -0.77), (0.08, -0.40), (0.03, -0.24), (0.00, -0.19)]
CROWN_IN = [(0.00, 0.54), (0.11, 0.70), (0.23, 0.79), (0.35, 0.70), (0.43, 0.54),
            (0.47, 0.35), (0.44, 0.15), (0.42, 0.00)]
PULP_PTS = [(0.00, 0.40), (0.07, 0.49), (0.12, 0.53), (0.16, 0.38), (0.175, 0.18),
            (0.155, 0.00), (0.13, -0.44), (0.095, -0.88), (0.07, -1.02), (0.045, -0.88),
            (0.035, -0.50), (0.03, -0.16), (0.00, -0.08)]

def _mirror(half):
    half = list(half)
    return half + [(-x, y) for x, y in reversed(half[1:-1])]

def _spline(pts, n, closed):
    A = np.asarray(pts, float)
    if closed and (A[0] != A[-1]).any():
        A = np.vstack([A, A[:1]])
    try:
        from scipy.interpolate import splprep, splev
        tck, _u = splprep([A[:, 0], A[:, 1]], s=0, per=bool(closed), k=3)
        xs, ys = splev(np.linspace(0, 1, n), tck)
        return list(zip(xs, ys))
    except Exception:
        return [tuple(p) for p in A]

def _tooth_paths():
    body = _spline(_mirror(CROWN_OUT + ROOT_OUT), 480, True)
    outer = _spline([(-x, y) for x, y in reversed(CROWN_OUT)] + CROWN_OUT[1:], 220, False)
    inner = _spline(list(reversed(CROWN_IN)) + [(-x, y) for x, y in CROWN_IN[1:]], 220, False)
    pulp = _spline(_mirror(PULP_PTS), 320, True)
    return body, outer + inner, pulp

def tooth_panel(ax):
    from matplotlib.patches import Polygon
    body, shell, pulp = _tooth_paths()
    rx = [0.455, 0.405, 0.325, 0.255]
    ry = [-0.04, -0.58, -1.00, -1.24]
    for sgn in (-1, 1):
        ax.plot([sgn * v for v in rx], ry, color=SOCKET_F, lw=5.2,
                solid_capstyle="round", zorder=0)
    ax.add_patch(Polygon(body, closed=True, facecolor=DENTINE_F, edgecolor=MUTED,
                         lw=0.9, zorder=2))
    ax.add_patch(Polygon(shell, closed=True, facecolor=ENAMEL_F, edgecolor=MUTED,
                         lw=0.9, zorder=3))
    ax.add_patch(Polygon(pulp, closed=True, facecolor=PULP_F, edgecolor=MUTED,
                         lw=0.7, zorder=4))
    for sgn in (-1, 1):
        ax.plot([sgn * v for v in rx], ry, color=S4, lw=0.9, alpha=0.75,
                solid_capstyle="round", zorder=5)
    ax.plot([-0.47, 0.47], [0.00, 0.00], color=MUTED, lw=0.6, ls=(0, (2, 2)), zorder=6)

    groups = tissue_groups()
    gx, gy, gsc = 1.06, -0.40, 0.30
    ax.add_patch(Polygon([(gx + x * gsc, gy + y * gsc) for x, y in body], closed=True,
                         facecolor="none", edgecolor=MUTED, lw=0.8, ls=(0, (3, 2)),
                         alpha=0.85, zorder=2))
    ax.text(gx, gy - 0.52, "Tooth number (%d)" % len(groups.get("Tooth number", [])),
            fontsize=9.2, weight="bold", color=dict((t[0], t[2]) for t in TISSUE)["Tooth number"],
            ha="center", va="center", zorder=7)
    ax.text(gx, gy - 0.74, "a tooth that never forms", fontsize=7.2, color=INK2,
            style="italic", ha="center", va="center", zorder=7)

    anchors = {"Enamel": (0.47, 0.62, 0.88, 0.82),
               "Dentine": (-0.33, 0.42, -0.88, 0.64),
               "Periodontium": (-0.43, -0.66, -0.88, -0.54)}
    for name, _pat, col in TISSUE:
        genes = groups.get(name, [])
        if not genes or name not in anchors:
            continue
        xa, ya, xt, yt = anchors[name]
        ha = "left" if xt > 0 else "right"
        ax.annotate("", xy=(xa, ya), xytext=(xt + (0.04 if xt > 0 else -0.04), yt),
                    arrowprops=dict(arrowstyle="-", color=col, lw=0.9, shrinkA=0,
                                    shrinkB=2), zorder=6)
        ax.text(xt + (0.06 if xt > 0 else -0.06), yt, "%s (%d)" % (name, len(genes)),
                fontsize=9.2, weight="bold", color=col, ha=ha, va="center", zorder=7)
    ax.annotate("pulp", xy=(0.13, 0.42), xytext=(0.62, 0.28), fontsize=7.4,
                color=INK2, style="italic", ha="left", va="center",
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.7, shrinkA=1,
                                shrinkB=2), zorder=7)
    ax.set_xlim(-1.70, 1.72); ax.set_ylim(-1.38, 1.12); ax.axis("off")
    ax.set_aspect("equal", adjustable="box")

def figure6():
    feats = json.load(open(os.path.join(DATASET, "dspp_features.json")))
    Lp = feats["length"]
    dom = {f["desc"] or f["type"]: (f["start"], f["end"]) for f in feats["features"]}
    sig = next((f for f in feats["features"] if f["type"] == "Signal"), None)
    dsp = next((f for f in feats["features"] if f["desc"] == "Dentin sialoprotein"), None)
    dpp = next((f for f in feats["features"] if f["desc"] == "Dentin phosphoprotein"), None)
    dis = max((f for f in feats["features"] if f["desc"] == "Disordered"),
              key=lambda f: f["end"] - f["start"], default=None)

    rows = [r for r in VT if r["gene"] == "DSPP"]
    esc = next((r["nmd_escape_from_codon"] for r in rows if r["nmd_escape_from_codon"]), None)
    path_fs = sorted([r for r in rows if r["source"] == "clinvar_pathogenic"
                      and r["variant_class"] == "frameshift" and r["codon"]],
                     key=lambda r: r["codon"])
    path_ns = sorted([r for r in rows if r["source"] == "clinvar_pathogenic"
                      and r["variant_class"] == "nonsense" and r["codon"]],
                     key=lambda r: r["codon"])
    pop = sorted([r for r in rows if r["source"] == "gnomad_plof" and r["codon"]],
                 key=lambda r: r["codon"])

    fig, ax = plt.subplots(figsize=(7.4, 5.2))
    top = len(path_fs) + 3.0
    if esc:
        ax.add_patch(Rectangle((esc, -len(pop) - 2.2), Lp - esc, top + len(pop) + 2.2,
                               facecolor="#f4edf8", edgecolor="none", zorder=0))
        ax.axvline(esc, color=MUTED, lw=0.9, ls=(0, (3, 2)), zorder=1)
        ax.text(esc + 12, top + 0.45, "escape window", fontsize=7.6, color=INK2, va="top")
    yb = top - 1.0
    for f, col, lab in ((sig, GREY, "signal"), (dsp, S1, "dentin sialoprotein"),
                        (dpp, S2, "dentin phosphoprotein")):
        if not f:
            continue
        ax.add_patch(Rectangle((f["start"], yb), f["end"] - f["start"], 0.9,
                               facecolor=col, alpha=0.30, edgecolor=col, lw=0.9, zorder=3))
        if f["end"] - f["start"] > 120:
            ax.text((f["start"] + f["end"]) / 2, yb + 0.45, lab, fontsize=7.6, ha="center",
                    va="center", color=INK, weight="bold", zorder=4)
        else:
            ax.annotate(lab, xy=((f["start"] + f["end"]) / 2, yb + 0.9),
                        xytext=(26, yb + 1.55), textcoords=("offset points", "data"),
                        fontsize=7.0, color=INK2, ha="left",
                        arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.6))
    if dis:
        ax.plot([dis["start"], dis["end"]], [yb - 0.42, yb - 0.42], color=MUTED, lw=2.2,
                solid_capstyle="butt", zorder=3)
        ax.text((dis["start"] + dis["end"]) / 2, yb - 0.72, "disordered",
                fontsize=7.0, ha="center", va="top", color=INK2)

    for i, r in enumerate(path_fs):
        y = top - 2.6 - i * 0.30
        tail = r["fsd"] or 0
        ax.plot([r["codon"], r["codon"] + tail], [y, y], color=S2, lw=1.5,
                solid_capstyle="butt", alpha=0.85, zorder=3)
        ax.plot([r["codon"]], [y], marker="|", color=INK, markersize=4.5,
                markeredgewidth=1.0, zorder=4)
    for r in path_ns:
        y = top - 2.6 - len(path_fs) * 0.30 - 0.45
        ax.plot([r["codon"]], [y], marker="v", color=S1, markersize=6.5,
                markeredgecolor=WHITE, markeredgewidth=0.6, zorder=4)
    ax.text(Lp + 60, top - 2.6 - (len(path_fs) - 1) * 0.15,
            f"{len(path_fs)} frameshifts", fontsize=7.8, color=S2,
            va="center", ha="left")
    ax.text(Lp + 60, top - 2.6 - len(path_fs) * 0.30 - 0.45,
            f"{len(path_ns)} nonsense", fontsize=7.8, color=S1, va="center", ha="left")

    ybase = top - 2.6 - len(path_fs) * 0.30 - 1.5
    for i, r in enumerate(pop):
        y = ybase - (i % 18) * 0.17
        tail = r["fsd"] or 0
        ax.plot([r["codon"], r["codon"] + tail], [y, y], color=GREY, lw=1.1,
                solid_capstyle="butt", alpha=0.85, zorder=2)
        ax.plot([r["codon"]], [y], marker="|", color=MUTED, markersize=3.6,
                markeredgewidth=0.8, zorder=3)
    ax.text(Lp + 60, ybase - 1.5, f"{len(pop)} gnomAD", fontsize=7.8,
            color=INK2, va="center", ha="left")

    ax.set_xlim(-60, Lp + 450)
    ax.set_ylim(ybase - 4.3, top + 1.1)
    ax.set_yticks([])
    ax.set_xlabel("residue")
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.set_title("DSPP truncations and their novel tails", loc="left", color=INK,
                 pad=TITLE_PAD)
    ax.plot([Lp, Lp], [ybase - 3.35, top - 1.35], color=MUTED, lw=0.8, ls=(0, (2, 2)),
            zorder=1)
    ax.text(Lp, ybase - 3.55, "original stop", fontsize=7.2, color=INK2, ha="center",
            va="top")
    save(fig, "Figure5")

PLAN = {"1": figure1, "2": figure2, "3": figure4, "4": figure5, "5": figure6,
        "S1": figure3, "S2": figure6_grid}


if __name__ == "__main__":
    want = sys.argv[1:] or ["1", "2", "3", "4", "5", "S1", "S2"]
    for w in want:
        PLAN[w]()
