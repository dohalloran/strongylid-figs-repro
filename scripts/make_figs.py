from __future__ import annotations

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from common import ensure_unpacked, load_manifest, default_manifest_path, load_genus_table, nonhost_reads_from_kreport, FOCAL_GENERA

# Match manuscript figure colors (approx.)
GENUS_COLORS = {
    "Necator": "#56B4E9",          # light blue
    "Ancylostoma": "#E69F00",      # orange
    "Oesophagostomum": "#009E73",  # green
    "Trichostrongylus": "#F0E442", # yellow
}

def fig2_stacked(genus_df: pd.DataFrame, manifest: pd.DataFrame, out_png: Path) -> None:
    # Use genus-level Bracken table: fractions already sum to 1.0 per SRR for the 4 focal genera
    g = genus_df[genus_df["name"].isin(FOCAL_GENERA)].copy()

    srr_order = manifest["SRR"].tolist()
    # matrix rows = genus, cols = SRR
    mat = (
        g.pivot_table(index="name", columns="SRR", values="fraction_total_reads", aggfunc="first")
          .reindex(index=FOCAL_GENERA, columns=srr_order)
          .fillna(0.0)
    )

    fig, ax = plt.subplots(figsize=(18, 5))
    x = np.arange(len(srr_order))
    bottom = np.zeros(len(srr_order), dtype=float)

    # stack order to match the PDF (bottom: Ancylostoma, then Necator, Oeso, Tricho)
    stack_order = ["Ancylostoma","Necator","Oesophagostomum","Trichostrongylus"]
    for genus in stack_order:
        vals = mat.loc[genus].to_numpy(dtype=float)
        ax.bar(x, vals, bottom=bottom, label=genus, color=GENUS_COLORS.get(genus, None), linewidth=0)
        bottom += vals

    # host blocks + divider
    chimp_n = int((manifest["host"] == "chimpanzee").sum())
    ax.axvline(chimp_n - 0.5, color="black", linewidth=1)
    ax.text((chimp_n - 1) / 2, 1.02, "Chimpanzee", ha="center", va="bottom", transform=ax.get_xaxis_transform())
    ax.text(chimp_n + (len(srr_order) - chimp_n - 1) / 2, 1.02, "Gorilla", ha="center", va="bottom", transform=ax.get_xaxis_transform())

    ax.set_ylabel("Proportion of strongylid reads")
    ax.set_xlabel("Sample (SRR), grouped by host")
    ax.set_ylim(0, 1.05)
    ax.set_xticks([])

    leg = ax.legend(title="Genus", frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    # tight layout with room for legend
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

def fig4_richness(genus_df: pd.DataFrame, manifest: pd.DataFrame, out_png: Path) -> None:
    g = genus_df[genus_df["name"].isin(FOCAL_GENERA)].copy()
    # thresholds and RA cutoff
    ra_cut = 0.001  # 0.1%
    thresholds = list(range(0, 2501, 50))

    rows = []
    for thr in thresholds:
        for srr in manifest["SRR"]:
            sub = g[g["SRR"] == srr].copy()
            keep = sub[(sub["fraction_total_reads"] >= ra_cut) & (sub["new_est_reads"] >= thr)]
            rows.append({
                "SRR": srr,
                "host": manifest.loc[manifest["SRR"] == srr, "host"].values[0],
                "threshold": thr,
                "richness": int(keep["name"].nunique()),
            })
    df = pd.DataFrame(rows)

    mean_rich = df.groupby(["host","threshold"])["richness"].mean().reset_index()

    fig, ax = plt.subplots(figsize=(12, 5))

    # chimp: dashed + circles; gorilla: solid + triangles
    chimp = mean_rich[mean_rich["host"]=="chimpanzee"]
    gor = mean_rich[mean_rich["host"]=="gorilla"]

    ax.plot(chimp["threshold"], chimp["richness"], linestyle="--", marker="o", markersize=4, linewidth=1.8, label="chimpanzee mean")
    ax.plot(gor["threshold"], gor["richness"], linestyle="-", marker="^", markersize=5, linewidth=1.8, label="gorilla mean")

    ax.axvline(600, linestyle=":", linewidth=1.5, color="black")
    ax.text(610, 3.02, "First drop (chimp): 600 reads", va="bottom", ha="left")

    ax.set_xlabel("Minimum reads per genus (beyond RA≥0.1%)")
    ax.set_ylabel("Mean strongylid genus richness")
    ax.set_xlim(0, 2500)
    ax.set_ylim(2.8, 4.05)
    ax.set_xticks([0,500,1000,1500,2000,2500])
    ax.legend(frameon=False, loc="lower left")

    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

def fig5_pri(genus_df: pd.DataFrame, manifest: pd.DataFrame, kreports_dir: Path, out_png: Path) -> None:
    g = genus_df[genus_df["name"].isin(FOCAL_GENERA)].copy()
    ra_cut = 0.001
    min_reads = 100

    vals = []
    for _, row in manifest.iterrows():
        srr = row["SRR"]
        host = row["host"]
        sub = g[g["SRR"] == srr]
        keep = sub[(sub["fraction_total_reads"] >= ra_cut) & (sub["new_est_reads"] >= min_reads)]
        numerator = float(keep["new_est_reads"].sum())

        kreport = kreports_dir / f"{srr}.kreport"
        denom = float(nonhost_reads_from_kreport(kreport))
        pri_rpm = (numerator / denom) * 1e6
        vals.append((host, pri_rpm))
    df = pd.DataFrame(vals, columns=["host","PRI_RPM"])

    chimp = df[df["host"]=="chimpanzee"]["PRI_RPM"].to_numpy()
    gor = df[df["host"]=="gorilla"]["PRI_RPM"].to_numpy()

    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.set_title("Relative strongylid read burden")
    ax.boxplot([chimp, gor], labels=["Chimpanzee","Gorilla"])

    # overlay points (match PDF markers/colors)
    rng = np.random.default_rng(0)
    ax.scatter(1 + rng.normal(0, 0.05, size=len(chimp)), chimp, s=40, color="#1f77b4", marker="o")
    ax.scatter(2 + rng.normal(0, 0.05, size=len(gor)), gor, s=45, color="#ff7f0e", marker="^")

    ax.set_ylabel("Strongylid PRI_RPM")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close(fig)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--share-zip", required=True, type=Path, help="Path to share.zip")
    ap.add_argument("--workdir", default=Path("data"), type=Path, help="Working directory for unpacking (default: data/)")
    ap.add_argument("--manifest", default=None, type=Path, help="Path to srr_manifest.txt (tab-delimited). If omitted, tries inside share bundle.")
    ap.add_argument("--outdir", default=Path("outputs/figures"), type=Path, help="Output directory for figures")
    args = ap.parse_args()

    paths = ensure_unpacked(args.share_zip, args.workdir)
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)

    mp = args.manifest if args.manifest else default_manifest_path(paths)
    manifest = load_manifest(mp)
    genus_df = load_genus_table(paths)

    fig2_stacked(genus_df, manifest, outdir/"Figure2.png")
    fig4_richness(genus_df, manifest, outdir/"Figure4.png")
    fig5_pri(genus_df, manifest, paths.kreports_dir, outdir/"Figure5.png")

    print("Wrote Figure2.png, Figure4.png, Figure5.png to", outdir)

if __name__ == "__main__":
    main()
