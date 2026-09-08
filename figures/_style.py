"""House style shared by the three preprint figures.

Colours: one per reference engine plus two greys. Fonts: DejaVu Sans, 8 pt
body, 7 pt ticks. PDF text stays text (Type 42).
"""
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = pathlib.Path(__file__).resolve().parent

SANDER = "#2a78d6"
OPENMM = "#eb6834"
GROMACS = "#1baf7a"
BLUE = SANDER
TEXT = "#52514e"
GRID = "#d9d8d4"
LIGHT = "#a8a7a3"          # reference lines and secondary marks

ENGINE_COLOUR = {"sander": SANDER, "openmm": OPENMM, "gromacs": GROMACS}
ENGINE_MARKER = {"sander": "o", "openmm": "s", "gromacs": "^"}
ENGINE_NAME = {"sander": "sander", "openmm": "OpenMM", "gromacs": "GROMACS"}


def apply():
    plt.rcParams.update({
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
        "font.size": 8,
        "axes.labelsize": 8, "axes.titlesize": 8,
        "xtick.labelsize": 7, "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "xtick.minor.size": 1.5, "ytick.minor.size": 1.5,
        "figure.facecolor": "white", "axes.facecolor": "white",
        "savefig.facecolor": "white",
        "text.color": TEXT, "axes.labelcolor": TEXT,
        "xtick.color": TEXT, "ytick.color": TEXT,
        "axes.edgecolor": TEXT,
        "grid.color": GRID, "grid.linewidth": 0.5,
        "axes.spines.top": False, "axes.spines.right": False,
        "mathtext.default": "regular",
    })


def ygrid(ax):
    ax.grid(True, axis="y", which="major", color=GRID, lw=0.5, zorder=0)
    ax.grid(False, axis="x")
    ax.set_axisbelow(True)


def letter(ax, text, dx=-0.02, dy=0.02, fig=None):
    """Bold panel letter at the top-left corner, outside the axes.

    dx, dy are offsets in figure inches from the axes' top-left corner.
    """
    fig = fig or ax.figure
    ax.annotate(text, xy=(0, 1), xycoords="axes fraction",
                xytext=(dx * 72, dy * 72), textcoords="offset points",
                ha="right", va="bottom", fontsize=8, fontweight="bold",
                color="#222222", annotation_clip=False)


def save(fig, base):
    base = pathlib.Path(base)
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.02)
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight",
                pad_inches=0.02)
    return base.with_suffix(".pdf"), base.with_suffix(".png")
