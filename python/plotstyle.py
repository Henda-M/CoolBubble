# CoolBubble -- Copyright (C) 2026 Henda Mansour, Yann Gouttenoire
# SPDX-License-Identifier: GPL-3.0-or-later (see the LICENSE file)
"""Common look of all CoolBubble plots. Edit this file to change the style everywhere.

The preferred fonts (Palatino-like text, TeX Gyre Pagella for math) are used
when they are installed; otherwise matplotlib's default fonts are kept.
"""

import matplotlib.pyplot as plt
from matplotlib import font_manager

TEXT_FONTS = ["URW Palladio L", "P052"]  # first one found is used
MATH_FONT = "TeX Gyre Pagella"

SIZES = {
    "font.size": 16,
    "axes.labelsize": 17,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 13,
    "figure.titlesize": 16,
    "grid.linewidth": 0.2,
    "grid.alpha": 0.5,
    "savefig.bbox": "tight",   # no labels cut off, no tight_layout() needed
    "savefig.dpi": 150,
}


def apply():
    """Set the CoolBubble style for all following plots."""
    plt.rcParams.update(SIZES)

    available = {f.name for f in font_manager.fontManager.ttflist}
    text_font = next((f for f in TEXT_FONTS if f in available), None)
    if text_font:
        plt.rcParams["font.family"] = text_font

    if MATH_FONT in available:
        def pick(font, fallback):
            return font if font in available else fallback
        plt.rcParams.update({
            "mathtext.fontset": "custom",
            "mathtext.rm": MATH_FONT,
            "mathtext.it": f"{MATH_FONT}:italic",
            "mathtext.bf": f"{MATH_FONT}:bold",
            "mathtext.sf": pick("TeX Gyre Heros", "DejaVu Sans"),
            "mathtext.tt": "DejaVu Sans Mono",
            "mathtext.cal": pick("TeX Gyre Chorus", "DejaVu Sans"),
        })
    elif text_font:
        plt.rcParams["mathtext.fontset"] = "stix"  # serif math to match a serif text font
