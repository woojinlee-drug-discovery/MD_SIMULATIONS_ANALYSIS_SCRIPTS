#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re, glob, os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns

# ------------------ Config ------------------
YIELD_FILE = "experimental_yields.txt"
OUTPUT_DIR = "plots_per_file"
VMAX_YIELD = 200.0  # clamp yields at this value for color scale

# ------------------ Pair sets (as numeric tuples) ------------------
ICL_TUPLES = {
    (264, 63), (266, 63), (273, 63),
    (266, 65), (267, 65), (268, 65), (270, 65), (271, 65), (272, 65), (273, 65),
    (266, 67), (268, 67), (270, 67), (271, 67), (272, 67), (273, 67),
    (268, 73), (269, 73), (270, 73), (271, 73), (272, 73), (273, 73),
    (268, 74), (269, 74), (270, 74), (271, 74), (272, 74), (273, 74),
    (264, 75), (267, 75), (268, 75), (269, 75), (270, 75), (272, 75), (273, 75),
    (268, 76), (269, 76),
    (268, 77), (270, 77), (271, 77), (272, 77),
    (267, 78), (268, 78), (269, 78), (271, 78), (273, 78),
    (267, 81), (268, 81), (269, 81), (270, 81),
    (267, 83), (268, 83), (269, 83), (270, 83),
    (101, 135), (173, 249), (176, 249), (263, 313), (263, 315)
}

CTERM_TUPLES = {
    (349, 13), (350, 13), (358, 13), (359, 13), (360, 13), (362, 13), (363, 13), (364, 13),
    (365, 13), (366, 13), (367, 13), (368, 13), (369, 13), (370, 13), (371, 13), (372, 13),
    (373, 13), (358, 14), (359, 14), (360, 14), (361, 14), (362, 14), (363, 14), (364, 14),
    (365, 14), (366, 14), (367, 14), (368, 14), (369, 14), (370, 14), (371, 14), (372, 14),
    (373, 14), (358, 15), (359, 15), (360, 15), (361, 15), (362, 15), (363, 15), (364, 15),
    (365, 15), (366, 15), (367, 15), (368, 15), (369, 15), (370, 15), (371, 15), (372, 15),
    (367, 16), (369, 16), (370, 16), (371, 16), (372, 16), (374, 16),
    (381, 99), (382, 99), (383, 99), (384, 99), (385, 99), (386, 99), (388, 99), (391, 99),
    (375, 103), (376, 103), (379, 103), (380, 103), (381, 103), (382, 103), (383, 103),
    (384, 103), (385, 103), (386, 103), (388, 103), (389, 103), (391, 103),
    (373, 104), (375, 104), (376, 104), (379, 104), (380, 104), (381, 104), (382, 104),
    (383, 104), (384, 104), (385, 104), (386, 104), (388, 104),
    (383, 106), (384, 106), (385, 106), (386, 106), (387, 106), (388, 106),
    (358, 159), (359, 159), (360, 159), (361, 159), (362, 159), (363, 159), (364, 159),
    (358, 160), (359, 160), (360, 160), (361, 160), (362, 160), (363, 160), (364, 160),
    (365, 160), (366, 160), (367, 160), (368, 160), (369, 160), (370, 160), (371, 160),
    (372, 160), (373, 160), (374, 160),
    (343, 165), (347, 165), (348, 165), (366, 165), (365, 165), (367, 165), (368, 165),
    (369, 165), (370, 165), (371, 165), (372, 165), (373, 165),
    (365, 295), (369, 295), (370, 295), (371, 295), (372, 295)
}

# ------------------ Regex (Chain A/B, CA/CB, Å or A) ------------------
YIELD_RE = re.compile(
    r"Resid (\d+) \(([A-Z?])-C[AB], Chain [AB]\) -> Resid (\d+) \(([A-Z?])-C[AB], Chain [AB]\) : ([\d.]+)%",
    re.UNICODE
)
DIST_RE = re.compile(
    r"Resid (\d+) \(([A-Z?])-C[AB], Chain [AB]\) -> Resid (\d+) \(([A-Z?])-C[AB], Chain [AB]\) : ([\d.]+)\s*(?:Å|A)",
    re.UNICODE
)

# ------------------ Parsers ------------------
def parse_yields(path: str) -> dict:
    y = {}
    if not os.path.exists(path):
        print(f"[info] No yields file: {path} (using 0%)")
        return y
    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            m = YIELD_RE.search(line)
            if m:
                r1, aa1, r2, aa2, val = m.groups()
                y[f"{aa1}{r1}-{aa2}{r2}"] = float(val)
    print(f"[ok] Loaded {len(y)} yield pairs")
    return y

def parse_distance_file(file_path: str, ymap: dict) -> pd.DataFrame:
    pairs, dists, ylds, tuples = [], [], [], []
    with open(file_path, "r", encoding="utf-8") as df:
        for line in df:
            m = DIST_RE.search(line)
            if not m:
                continue
            r1, aa1, r2, aa2, d = m.groups()
            pair = f"{aa1}{r1}-{aa2}{r2}"
            try:
                dval = float(d)
            except ValueError:
                continue
            pairs.append(pair)
            dists.append(dval)
            ylds.append(float(ymap.get(pair, 0.0)))
            tuples.append((int(r1), int(r2)))

    out = pd.DataFrame({
        "Residue Pair": pairs,
        "PairTuple": tuples,         # numeric (left_resid, right_resid)
        "Distance (Å)": dists,
        "Yield (%)": ylds
    })
    return out

# ------------------ Plotting helpers ------------------
def build_pair_colors(df: pd.DataFrame, vmax: float = VMAX_YIELD):
    sns.set(style="whitegrid")
    cmap = sns.color_palette("YlOrRd", as_cmap=True)
    norm = mcolors.Normalize(vmin=0, vmax=vmax)

    def yield_to_color(y):
        y = max(0.0, min(vmax, float(y)))
        return cmap(y / vmax)

    pair_to_color = {}
    for pair, sub in df.groupby("Residue Pair"):
        pair_to_color[pair] = yield_to_color(sub["Yield (%)"].iloc[0])

    return pair_to_color, cmap, norm

def left_resid_from_key(pair_key: str):
    m = re.match(r"[A-Za-z?](\d+)-[A-Za-z?]\d+$", pair_key)
    return int(m.group(1)) if m else None

def plot_subset(df: pd.DataFrame, title_suffix: str, out_png: str):
    if df.empty:
        print(f"[warn] No rows for {title_suffix} → skip {out_png}")
        return

    df = df.copy()
    df["LeftResid"] = df["Residue Pair"].apply(left_resid_from_key)
    order = (df.drop_duplicates(subset=["Residue Pair"])
               .sort_values(by="LeftResid", ascending=True)["Residue Pair"]
               .tolist())

    pair_to_color, cmap, norm = build_pair_colors(df)

    plt.figure(figsize=(10, 14))
    ax = sns.stripplot(
        x="Distance (Å)",
        y="Residue Pairs",  # temporary label for pretty axis name
        data=df.rename(columns={"Residue Pair": "Residue Pairs"}),
        order=order,
        hue="Residue Pairs",
        palette=pair_to_color,
        orient="h",
        size=8,
        jitter=False
    )
    if ax.legend_ is not None:
        ax.legend_.remove()

    plt.xlabel("Distance (Å)", fontsize=12, fontweight="bold", labelpad=10)
    plt.ylabel(f"Residue Pairs • {title_suffix}", fontsize=12, fontweight="bold")
    for ref in [0, 5, 10, 15, 20, 25, 30]:
        plt.axvline(ref, linestyle="--", color="gray", alpha=0.5)

    # Colorbar for yields
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cax = inset_axes(ax, width="30%", height="1%", loc="lower right", borderpad=2)
    cbar = plt.colorbar(sm, cax=cax, orientation='horizontal')
    cbar.set_label("Exp. Yield (%)", fontsize=10)
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.tick_params(labelsize=8)

    ax.xaxis.set_label_position('top')
    ax.xaxis.set_ticks_position('top')
    plt.grid(axis='y', linestyle="--", alpha=0.3)

    os.makedirs(os.path.dirname(out_png), exist_ok=True)
    plt.tight_layout(pad=0.5)
    plt.savefig(out_png, dpi=300)
    plt.close()
    print(f"[ok] Saved: {out_png}")

# ------------------ Main ------------------
def main():
    files = sorted([f for f in glob.glob("*.txt") if os.path.basename(f) != YIELD_FILE])
    if not files:
        print("[err] No distance *.txt files found.")
        return

    ymap = parse_yields(YIELD_FILE)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    for f in files:
        df = parse_distance_file(f, ymap)
        base = os.path.splitext(os.path.basename(f))[0]

        # Filter to subsets
        df_icl   = df[df["PairTuple"].isin(ICL_TUPLES)]
        df_cterm = df[df["PairTuple"].isin(CTERM_TUPLES)]

        # Plot each subset separately
        out_icl   = os.path.join(OUTPUT_DIR, f"{base}_ICL.png")
        out_cterm = os.path.join(OUTPUT_DIR, f"{base}_CTerm.png")
        plot_subset(df_icl,   "ICL",   out_icl)
        plot_subset(df_cterm, "C-terminus", out_cterm)

if __name__ == "__main__":
    main()
