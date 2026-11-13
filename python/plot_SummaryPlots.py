# CoolBubble -- Copyright (C) 2026 Henda Mansour, Yann Gouttenoire
# SPDX-License-Identifier: GPL-3.0-or-later (see the LICENSE file)
"""Summary plots of a finished CoolBubble run.

Called at the end of every run, or by hand:
    python3 python/plot_SummaryPlots.py --folder <run folder>

Plots (written into the run folder):
    Energies           energy densities over time, and the drift of the total energy
    EoS                equation of state w = P/rho over time (and over a, if the box expands)
    Averages           box averages of phi and phi'
    PowerSpectrum      power spectrum at several times (and in comoving units, if the box expands)
    Expansion          scale factor and Hubble constraint (expanding runs only)
    Gamma              Lorentz factor of the wall (wall tracking only)
    FieldEvolution     phi(x, t) (1D only)
    Ekin__phi_t=NNNN   2D: kinetic energy density at every snapshot         } one log colour scale for all
    Egrad__phi_t=NNNN  2D: gradient energy density (if memorysave = 0)      } snapshots, for videos

Times and momenta are shown in units of the true-vacuum mass m_phi = sqrt(2 (c + 6)).
"""

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # write files only, no window
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import seaborn as sns
from scipy.signal import savgol_filter

import plotstyle
from runfiles import mass, read_parameters, read_snapshot, read_spectrum

plotstyle.apply()
palette = sns.color_palette("rocket", n_colors=10)

# decades below the maximum shown by the logarithmic colour scale of the 2D energy maps
ENERGY_LOG_DECADES = 4

# columns of EnergyDensities.dat (see the file header)
T, A, ETOT, EKIN, EGRAD, EPOT = 0, 1, 2, 3, 4, 5
AVPHI, AVPHIP, HUBBLE, CONSTRAINT = 10, 12, 14, 15


def smooth(y):
    """Smooth fast oscillations (1D runs and long runs), keeping short series unchanged."""
    if not (p["D"] == 1 or len(y) > 5000):
        return y
    window = len(y) // 60
    window += 1 - window % 2  # odd window length
    if window <= 3:
        return y
    return savgol_filter(y, window, 3)


def save(name):
    plt.savefig(folder / f"{name}.png")
    plt.savefig(folder / f"{name}.pdf")
    plt.close()


#######################################################################
# Plots
#######################################################################

def plot_energies():
    t = E[:, T] * m
    fig, (ax, ax_err) = plt.subplots(2, 1, figsize=(6, 6), sharex=True,
                                     gridspec_kw={"height_ratios": [3, 1]})
    ax.plot(t, smooth(E[:, ETOT]), color=palette[0], lw=2.5, label=r"$\rho_\mathrm{tot}$")
    ax.plot(t, smooth(E[:, EKIN]), color=palette[2], lw=2.5, label=r"$\rho_\mathrm{kin}$")
    ax.plot(t, smooth(E[:, EGRAD]), color=palette[6], lw=2.5, label=r"$\rho_\mathrm{grad}$")
    ax.plot(t, smooth(E[:, EPOT]), color=palette[9], lw=2.5, label=r"$\rho_\mathrm{pot}$")
    ax.set_ylabel(r"$\rho/\Delta V$")
    ax.legend(ncol=2, handlelength=1)
    ax.grid(True)

    # energy conservation check (only meaningful without expansion)
    drift = E[:, ETOT] / E[0, ETOT] - 1
    ax_err.plot(t, drift, color=palette[0])
    ax_err.set_ylabel(r"$\Delta\rho_\mathrm{tot}/\rho_\mathrm{tot}$", fontsize=13)
    ax_err.set_xlabel(tlabel)
    ax_err.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
    ax_err.grid(True)
    ax_err.set_xlim(0, t.max())
    save("Energies")


def equation_of_state():
    """w = P/rho with P = rho_kin - (D-2)/D rho_grad - rho_pot."""
    D = p["D"]
    w = (E[:, EKIN] - (D - 2) / D * E[:, EGRAD] - E[:, EPOT]) / E[:, ETOT]
    return smooth(w)


def plot_eos():
    D = p["D"]
    w = equation_of_state()
    for x, xlabel, name in [(E[:, T] * m, tlabel, "EoS"), (E[:, A], r"$a$", "EoS_of_a")]:
        if name == "EoS_of_a" and not expanding:
            continue
        plt.figure(figsize=(6, 4))
        plt.plot(x, w, color=palette[0], lw=2.5)
        plt.axhline(1 / D, color=palette[2], lw=1, alpha=0.5)
        plt.axhline(0, color=palette[2], lw=1, alpha=0.5)
        plt.text(0.6, 1 / D + 0.03, "radiation-like", color=palette[2],
                 transform=plt.gca().get_yaxis_transform())
        plt.text(0.6, 0.03, "matter-like", color=palette[2],
                 transform=plt.gca().get_yaxis_transform())
        plt.ylim(-1, 1)
        plt.xlim(x.min(), x.max())
        plt.xlabel(xlabel)
        plt.ylabel(r"$w$")
        plt.grid(True)
        save(name)


def plot_averages():
    t = E[:, T] * m
    plt.figure(figsize=(6, 4))
    plt.plot(t, E[:, AVPHI], color=palette[1], label=r"$\langle\phi\rangle$")
    plt.plot(t, E[:, AVPHIP], color=palette[6], label=r"$\langle\dot\phi\rangle$")
    plt.xlabel(tlabel)
    plt.xlim(0, t.max())
    plt.grid(True)
    plt.legend()
    save("Averages")


def plot_power_spectra(comoving):
    """Power spectrum at up to 8 snapshot times, from early (light) to late (dark).
    comoving=False: against the physical momentum, comoving=True: against the comoving one
    (Delta itself is the same in both coordinates)."""
    D, dx = p["D"], p["L"] / p["nGrids"]
    files = sorted(folder.glob("_phi_k_t=*.dat"))[1:]  # the first one (t = 0) repeats the next
    if not files:
        return
    files = [files[i] for i in np.unique(np.linspace(0, len(files) - 1, min(8, len(files))).astype(int))]
    colors = sns.color_palette("Purples", n_colors=len(files) + 1)[1:]

    plt.figure(figsize=(6, 4))
    ymax, a_max = 0, 1.0
    for f, color in zip(files, colors):
        t, a, k_com, k_phys, Delta = read_spectrum(f)
        if comoving:
            x = k_com / m
            label = rf"$t = {t * m:.0f}\,m_\phi^{{-1}},\ a = {a:.2f}$"
        else:
            x = k_phys / m
            label = rf"$t = {t * m:.0f}\,m_\phi^{{-1}}$"
        plt.loglog(x, Delta, color=color, label=label)
        ymax, a_max = max(ymax, Delta.max()), max(a_max, a)

    # comoving lattice momenta; on the physical axis the smallest one moves down to kIR / a_max
    kIR, kUV = 2 * np.pi / p["L"], np.sqrt(D) * np.pi / dx
    plt.xlim(kIR / m / (1.0 if comoving else a_max), 0.5 * kUV / m)
    plt.ylim(ymax * 1e-5, ymax * 2)
    plt.xlabel(r"$k/m_\phi$" + (" (comoving)" if comoving else ""))
    plt.ylabel(r"$\Delta_\phi(k)$")
    plt.grid(True)
    plt.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=11)
    save("PowerSpectrum_comoving" if comoving else "PowerSpectrum")


def plot_expansion():
    t = E[:, T] * m
    fig, (ax_a, ax_c) = plt.subplots(2, 1, figsize=(6, 6), sharex=True)
    ax_a.plot(t, E[:, A], color=palette[1])
    ax_a.set_ylabel(r"$a$")
    ax_a.grid(True)
    H = E[:, HUBBLE]
    with np.errstate(divide="ignore", invalid="ignore"):
        ax_c.plot(t, E[:, CONSTRAINT] / H**2, color=palette[1])
    ax_c.set_ylabel(r"$C/H^2$")  # Hubble constraint, should stay small
    ax_c.set_xlabel(tlabel)
    ax_c.grid(True)
    ax_c.set_xlim(0, t.max())
    save("Expansion")


def plot_gamma():
    data = np.loadtxt(folder / "gamma.dat", ndmin=2)
    if data.size == 0:
        return
    plt.figure(figsize=(6, 4))
    plt.plot(data[:, 0] * m, data[:, 1], color=palette[3])
    plt.xlabel(tlabel)
    plt.ylabel(r"Lorentz factor $\gamma$")
    plt.ylim(0, data[:, 1].max() + 1)
    plt.grid(True)
    save("Gamma")


def plot_1d_field_evolution():
    """phi(x, t) for 1D runs, from the full field saved every Full1DFieldEvery steps."""
    data = np.loadtxt(folder / "FullField.dat", ndmin=2)
    x = np.linspace(0, p["L"] * m, data.shape[1])
    t = np.linspace(0, p["nt"] * p["TIME_STEP"] * m, data.shape[0])
    # at most 400 x 400 points, so the file stays small
    ix = np.linspace(0, len(x) - 1, min(400, len(x))).astype(int)
    it = np.linspace(0, len(t) - 1, min(400, len(t))).astype(int)

    plt.figure(figsize=(6, 4))
    mesh = plt.pcolormesh(x[ix], t[it], data[np.ix_(it, ix)], shading="auto",
                          cmap=sns.color_palette("mako", as_cmap=True), vmin=0, vmax=1.5, rasterized=True)
    plt.colorbar(mesh, label=r"$\phi(x,t)$", ticks=[0, 0.5, 1, 1.5])
    plt.xlabel(r"$x\ [m_\phi^{-1}]$")
    plt.ylabel(r"$t\ [m_\phi^{-1}]$")
    save("FieldEvolution")


def plot_2d_maps():
    """2D: maps of the kinetic energy density phi'^2/2 and of the gradient energy density at every snapshot.
    Both use one logarithmic colour scale, the same in all snapshots and for both quantities, so that they can
    be compared directly and put in a video: from the largest value down ENERGY_LOG_DECADES decades (the local
    energies grow strongly during the run). Smaller values get the lowest colour."""
    files = sorted(folder.glob("_phi_t=*.h5")) or sorted(folder.glob("_phi_t=*.dat"))

    def energies(f):
        """kinetic and gradient energy densities of snapshot f (when available) and its time"""
        fields, _, info = read_snapshot(f, p)
        maps = {}
        if "pi" in fields:
            maps["Ekin"] = 0.5 * fields["pi"] ** 2     # pi = phi' in 2D (no expansion)
        if "gradphi" in fields:
            maps["Egrad"] = fields["gradphi"]
        t = info.get("t")
        if t is None:  # older snapshots: time from the power spectrum file
            kfile = folder / f"_phi_k_t={f.stem.split('=')[-1]}.dat"
            t = read_spectrum(kfile)[0] if kfile.exists() else np.nan
        return maps, t

    # first pass: common colour scale
    vmax = 0.0
    for f in files:
        maps, _ = energies(f)
        for data in maps.values():
            vmax = max(vmax, data.max())
    if vmax <= 0:
        return
    norm = LogNorm(vmin=vmax * 10.0 ** (-ENERGY_LOG_DECADES), vmax=vmax)
    labels = {"Ekin": r"$\rho_\mathrm{kin}(x,y)$", "Egrad": r"$\rho_\mathrm{grad}(x,y)$"}

    # second pass: one map per quantity and snapshot
    for f in files:
        maps, t = energies(f)
        for name, data in maps.items():
            plt.figure(figsize=(6, 5))
            data = np.maximum(data.T, norm.vmin)   # zeros and tiny values: lowest colour
            im = plt.imshow(data, origin="lower", cmap="magma", norm=norm,
                            extent=(0, p["L"] * m, 0, p["L"] * m), aspect="equal")
            plt.colorbar(im, label=labels[name])
            plt.xlabel(r"$x\ [m_\phi^{-1}]$")
            plt.ylabel(r"$y\ [m_\phi^{-1}]$")
            plt.title(rf"$t = {t * m:.2f}\ m_\phi^{{-1}}$")
            plt.savefig(folder / f"{name}_{f.stem}.png")   # png only: one per snapshot
            plt.close()


#######################################################################
# Main
#######################################################################

parser = argparse.ArgumentParser(description="Summary plots of a CoolBubble run.")
parser.add_argument("--folder", required=True, type=Path, help="run folder (contains parameters.dat)")
folder = parser.parse_args().folder

p = read_parameters(folder)
E = np.loadtxt(folder / "EnergyDensities.dat", ndmin=2)
m = mass(p)                                   # mass in the true vacuum, unit of the plots
expanding = p.get("expansion", 0) == 1 and p["D"] == 3
tlabel = r"$t\ [m_\phi^{-1}]$" if p.get("ALPHAVAL", 0) == 0 else r"conformal time $\tau\ [m_\phi^{-1}]$"

plot_energies()
plot_eos()
plot_averages()
if p.get("FFT", 1) == 1:
    plot_power_spectra(comoving=False)
    if expanding:
        plot_power_spectra(comoving=True)
if expanding and E.shape[1] > CONSTRAINT:
    plot_expansion()
if p.get("WallTracking", 0) == 1 and (folder / "gamma.dat").exists():
    plot_gamma()
if p["D"] == 1 and (folder / "FullField.dat").exists():
    plot_1d_field_evolution()
if p["D"] == 2:
    plot_2d_maps()
