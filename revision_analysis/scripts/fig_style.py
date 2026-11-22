"""Figure style: Arial, editable PDF text, no top and right spines, black text on white.

Female and male are drawn in the colours used in the figures of the paper.
"""
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import font_manager

FONT_DIR = os.environ.get("FONT_DIR", "")   # optional folder with Arial font files

SEX_COL = {"Female": "#FFBF00", "Male": "#000080"}
SEXES = ["Female", "Male"]
ACCENT = "#CC3311"        # fitted lines, highlighted modules
ACCENT_2 = "#EE7733"
POINT = "#4477AA"         # highlighted genes
LIGHT = "#9ECAE1"         # distributions
GREY = "#7F7F7F"
VAC_COL = {"T1": "#9ECAE1", "T2": "#4292C6", "T3": "#08306B"}


def setup(base=7):
    if os.path.isdir(FONT_DIR):
        for f in font_manager.findSystemFonts(fontpaths=FONT_DIR):
            try:
                font_manager.fontManager.addfont(f)
            except Exception:
                pass
    sns.set_theme(style="ticks")
    mpl.rcParams.update({
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "font.family": "Arial", "font.sans-serif": ["Arial", "DejaVu Sans"],
        "font.size": base, "axes.titlesize": base + 1, "axes.labelsize": base,
        "xtick.labelsize": base - 0.5, "ytick.labelsize": base - 0.5,
        "legend.fontsize": base - 0.5, "legend.frameon": False, "legend.handletextpad": 0.5,
        "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5, "axes.titlepad": 4,
        "axes.unicode_minus": False,
        # all text pure black and fully opaque (seaborn's themes set text to ".15" dark grey).
        # Tick marks and spines keep the seaborn colour; only text is changed.
        "text.color": "black", "axes.labelcolor": "black", "axes.titlecolor": "black",
        "xtick.labelcolor": "black", "ytick.labelcolor": "black", "legend.labelcolor": "black",
        "savefig.transparent": False, "figure.facecolor": "white", "savefig.facecolor": "white",
        "figure.dpi": 150, "savefig.dpi": 400, "lines.linewidth": 1.0, "patch.linewidth": 0.5,
    })


def save(fig, outdir, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"{name}.{ext}"), bbox_inches="tight", dpi=400)
    plt.close(fig)
    print(f"  saved {name}.pdf + .png")


def sex_handles(**kw):
    return [plt.Line2D([], [], ls="", marker="o", ms=4, color=c, label=s, **kw)
            for s, c in SEX_COL.items()]
