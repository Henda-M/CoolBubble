# CoolBubble -- Copyright (C) 2026 Henda Mansour, Yann Gouttenoire
# SPDX-License-Identifier: GPL-3.0-or-later (see the LICENSE file)
"""Plots of one field snapshot and its power spectrum.

Called by the simulation after every snapshot, or by hand / by replot.py:
    python3 python/plot_FieldAndPowerSpectrum.py --folder <run folder> --filename _phi_t=0003.h5

Plots (written into the run folder, NNNN = snapshot number):
    _phi_t=NNNN.png      the field: phi(x) in 1D, phi(x, y) in 2D, a slice phi(x, z) in 3D
    Radial__phi_t=NNNN.png 2D/3D, single bubble (mode -2): the field along a line from the bubble centre
    Wall__phi_t=NNNN.png 1D with wall tracking: close-up of the wall
    _phi_k_t=NNNN.png    the power spectrum

With wall tracking the simulation also passes the wall radii --rmin, --rmax and
the Lorentz factor --gamma (in simulation units), which are marked in the plots.
Lengths and times are shown in units of 1/m_phi, momenta in units of m_phi.
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # write files only, no window
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy.signal import savgol_filter

import plotstyle
from runfiles import mass, read_parameters, read_snapshot, read_spectrum

plotstyle.apply()
palette = sns.color_palette("viridis", n_colors=10)
FIELD_CMAP = "magma"


def title():
    if np.isnan(t):  # time unknown (no power spectrum file, FFT = 0)
        return ""
    text = rf"$t = {t * m:.2f}\ m_\phi^{{-1}}$"
    if p.get("expansion", 0) == 1 and D == 3:
        text += rf"$,\ a = {a:.2f}$"
    return text


def save(name):
    plt.title(title())
    plt.savefig(folder / f"{name}.png")
    plt.close()


def bubbles():
    """True if the run started from bubbles (not from thermal fluctuations or the test wave)."""
    return p["readxns"] != -3 and p.get("TEST", 0) == 0


def wall_markers(ax, y=0.5):
    """Vertical lines at the wall radii r_min and r_max (wall tracking only)."""
    if not args.rmin:
        return
    for r, label in [(args.rmin, r"$r_\mathrm{min}$"), (args.rmax, r"$r_\mathrm{max}$")]:
        ax.axvline(r * m, color=palette[6], ls="--", lw=1)
        ax.text(r * m, y, label, fontsize=11, color=palette[6], rotation=90,
                transform=ax.get_xaxis_transform())


#######################################################################
# Field plots
#######################################################################

def plot_field_1d():
    x = axes[0] * m
    plt.figure(figsize=(6, 4))
    plt.plot(x, phi, color="black")
    wall_markers(plt.gca())
    plt.ylim(-0.2, 2)
    plt.xlim(0, p["L"] * m)
    plt.xlabel(r"$x\ [m_\phi^{-1}]$")
    plt.ylabel(r"$\phi(x)$")
    plt.grid(True)
    save(stem)


def plot_wall_1d():
    """Close-up of the wall around r_min .. r_max (1D wall tracking, before the bubbles collide)."""
    if not args.rmax or t > p["L"] / 2:
        return
    x = axes[0] * m
    width = (args.rmax - args.rmin) * m
    plt.figure(figsize=(6, 4))
    plt.plot(x, phi, color="black", lw=2)
    wall_markers(plt.gca(), y=0.8)
    plt.xlim(args.rmin * m - 8 * width, args.rmax * m + 5 * width)
    plt.ylim(-0.1, 1.5)
    plt.xlabel(r"$x\ [m_\phi^{-1}]$")
    plt.ylabel(r"$\phi(x)$")
    plt.grid(True)
    save(f"Wall_{stem}")


def image(data2d, ax0, ax1, label, vmin, vmax):
    """Colour map of data2d[i, j] with axis 0 horizontal and axis 1 vertical."""
    plt.figure(figsize=(6, 5))
    im = plt.imshow(data2d.T, origin="lower", cmap=FIELD_CMAP, vmin=vmin, vmax=vmax,
                    extent=(0, p["L"] * m, 0, p["L"] * m), aspect="equal")
    cbar = plt.colorbar(im, label=label)
    plt.xlabel(rf"${ax0}\ [m_\phi^{{-1}}]$")
    plt.ylabel(rf"${ax1}\ [m_\phi^{{-1}}]$")
    return cbar


def field_range(data):
    """Fixed colour range around the two vacua for bubbles, the data range otherwise."""
    return (0.0, 1.5) if bubbles() else (data.min(), data.max())


def plot_field_2d():
    vmin, vmax = field_range(phi)
    cbar = image(phi, "x", "y", r"$\phi(x,y)$", vmin, vmax)
    if bubbles():
        cbar.set_ticks([0, 0.5, 1, 1.5])
        cbar.set_ticklabels([r"$\phi_f = 0$", "0.5", r"$\phi_t = 1$", "1.5"])

    # wall tracking: circles at r_min and r_max around the tracked bubble
    centre = {-1: (p["L"] / 4, p["L"] / 4), -2: (p["L"] / 2, p["L"] / 2)}.get(p["readxns"])
    if args.rmin and centre:
        for r, color in [(args.rmin, "purple"), (args.rmax, "blue")]:
            plt.gca().add_patch(plt.Circle((centre[0] * m, centre[1] * m), r * m,
                                           color=color, fill=False, lw=2))
    save(stem)


def plot_field_3d():
    """Slice phi(x, z) at fixed y: through the centre for a single bubble, at y = L/4 otherwise
    (this cuts through the first bubble for two diagonal bubbles)."""
    n = phi.shape[1]
    j = n // 2 if p["readxns"] == -2 else n // 4
    y = axes[1][j] * m
    vmin, vmax = (-0.5, 1.5) if bubbles() else (phi.min(), phi.max())
    image(phi[:, j, :], "x", "z", rf"$\phi(x,\ y = {y:.1f},\ z)$", vmin, vmax)
    save(stem)


def plot_radial_profile():
    """2D/3D, single bubble (mode -2): the field along a line from the bubble centre (the middle of the box)
    outwards along the first axis, until the bubble reaches the box edge. Same line as the wall tracking."""
    if not (bubbles() and p["readxns"] == -2) or t > p["L"] / 2:
        return
    half = p["L"] / 2
    n = phi.shape[0]
    mid = n // 2                                    # saved point next to the centre on the other axes
    line = phi[(slice(mid, None),) + (mid,) * (D - 1)]
    # distance to the centre, including the offset of the line from the centre on the other axes
    offset2 = sum((ax[mid] - half) ** 2 for ax in axes[1:])
    r = np.sqrt((axes[0][mid:] - half) ** 2 + offset2)

    plt.figure(figsize=(6, 4))
    plt.plot(r * m, line, color="black")
    wall_markers(plt.gca())
    plt.xlim(0, half * m)
    plt.ylim(-0.2, 2)
    plt.xlabel(r"$r\ [m_\phi^{-1}]$")
    plt.ylabel(r"$\phi(r)$")
    plt.grid(True)
    save(f"Radial_{stem}")


#######################################################################
# Power spectrum
#######################################################################

def plot_power_spectrum():
    kfile = folder / f"_phi_k_t={number}.dat"
    if not kfile.exists():
        return
    _, a_k, _, k, delta = read_spectrum(kfile)   # k: physical momentum
    if D == 1:
        delta = savgol_filter(delta, 11, 3)

    dx = p["L"] / p["nGrids"]
    # largest and smallest momentum on the lattice, physical (comoving / a) like the k axis
    kUV = np.sqrt(D) * np.pi / dx / a_k
    kIR = 2 * np.pi / p["L"] / a_k
    ymax = delta.max()

    _, ax = plt.subplots(figsize=(6, 4))
    ax.loglog(k / m, delta, color=palette[4])

    def marker(k_over_m, label, color="grey", ls="--", y=0.05):
        ax.axvline(k_over_m, color=color, ls=ls, lw=1)
        ax.text(k_over_m * 1.05, y, label, fontsize=11, color=color, rotation=90,
                transform=ax.get_xaxis_transform())

    marker(kIR / m, r"$k_\mathrm{IR}$")
    marker(kUV / m, r"$k_\mathrm{UV}$")
    marker(1, r"$m_\phi$", color=palette[1], ls="-", y=0.75)
    if bubbles():
        # inverse radius of the nucleated bubble, and of a bubble that has grown at the speed of light
        R_ini = 0.823 * p["c"] ** 0.203 * (1.3 if D == 3 else 1.0)  # same fit as CRadius() in src/bubble.cpp
        marker(1 / R_ini / m, r"$R_\mathrm{ini}^{-1}$", y=0.75)
        if t < p["L"] / 2:
            marker(1 / (R_ini + t) / m, r"$R(t)^{-1}$", y=0.75)
    if args.gamma:
        marker(args.gamma, r"$\gamma\, m_\phi$", color=palette[2], ls="-", y=0.75)

    ax.set_xlim(0.5 * kIR / m, 2 * kUV / m)
    ax.set_ylim(ymax * 1e-10, ymax * 1e2)
    ax.set_xlabel(r"$k/m_\phi$")
    ax.set_ylabel(r"$\Delta_\phi(k)$")
    ax.grid(True)
    save(f"_phi_k_t={number}")


#######################################################################
# Main
#######################################################################

parser = argparse.ArgumentParser(description="Plot one field snapshot and its power spectrum.")
parser.add_argument("--folder", required=True, type=Path, help="run folder (contains parameters.dat)")
parser.add_argument("--filename", required=True, help="snapshot file, e.g. _phi_t=0003.h5")
# passed by the simulation only when WallTracking = 1 (simulation units)
parser.add_argument("--rmin", type=float)
parser.add_argument("--rmax", type=float)
parser.add_argument("--gamma", type=float)
args = parser.parse_args()

folder = args.folder
stem = Path(args.filename).stem                      # e.g. _phi_t=0003
number = stem.split("=")[-1]                         # e.g. 0003
p = read_parameters(folder)
D = p["D"]
m = mass(p)
fields, axes, info = read_snapshot(folder / args.filename, p)
# time and scale factor: stored in the snapshot; older snapshots only have them in the power spectrum file
kfile = folder / f"_phi_k_t={number}.dat"
if "t" in info:
    t, a = info["t"], info["a"]
else:
    t, a = read_spectrum(kfile)[:2] if kfile.exists() else (np.nan, 1.0)
phi = fields["phi"]

if D == 1:
    plot_field_1d()
    if p.get("WallTracking", 0) == 1:
        plot_wall_1d()
elif D == 2:
    plot_field_2d()
    plot_radial_profile()
else:
    plot_field_3d()
    plot_radial_profile()

if p.get("FFT", 1) == 1:
    plot_power_spectrum()
