# CoolBubble -- Copyright (C) 2026 Henda Mansour, Yann Gouttenoire
# SPDX-License-Identifier: GPL-3.0-or-later (see the LICENSE file)
"""Reading the files of a CoolBubble run folder."""

import re
from pathlib import Path

import numpy as np


def read_parameters(folder):
    """Return the 'key = value' pairs of <folder>/parameters.dat as a dict of numbers."""
    params = {}
    pattern = re.compile(r"(\w+)\s*=\s*([-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?)")
    for line in (Path(folder) / "parameters.dat").read_text().splitlines():
        if line.strip().startswith("#"):
            continue
        for key, value in pattern.findall(line):
            if key == "b":  # runs from before the b -> c rename
                key = "c"
            params[key] = float(value) if any(ch in value for ch in ".eE") else int(value)
    return params


def mass(params):
    """Mass of the field in the true vacuum, m_phi = sqrt(2 (c + 6)). Plots use units of 1/m_phi."""
    return np.sqrt(2 * (params["c"] + 6))


def read_snapshot(path, params):
    """Read a field snapshot _phi_t=NNNN.h5 (or .dat from older runs with ASCII output).

    Returns (fields, axes, info):
      fields  dict of arrays with shape (n0, n1, ...) indexed like the lattice [i, j, k];
              always 'phi', plus 'pi' and 'gradphi' when they were saved
      axes    list of D coordinate arrays (simulation units), one per lattice axis
      info    dict with the time 't' and scale factor 'a' of the snapshot, if stored in the file
    """
    D, L, nGrids = params["D"], params["L"], params["nGrids"]
    path = Path(path)
    fields, info = {}, {}

    if path.suffix == ".h5":
        import h5py
        with h5py.File(path, "r") as f:
            for name in ("phi", "pi", "gradphi"):
                if name in f:
                    fields[name] = f[name][:]
            info = {key: float(f.attrs[key]) for key in ("t", "a", "reduc") if key in f.attrs}
        phi = fields["phi"]
        # snapshots are stored with shape (n,)*D; older files hold the same values as one flat array
        n = phi.shape[0] if phi.ndim == D else int(round(phi.size ** (1.0 / D)))
        stride = int(info.get("reduc", nGrids // n))       # saved every stride-th lattice point per axis
        axis = (np.arange(n) * stride + 0.5) * L / nGrids  # cell centres, as in xFunc()
        axes = [axis] * D
        fields = {k: v.reshape((n,) * D) for k, v in fields.items()}
    else:
        # ASCII snapshots of older runs: columns x [y [z]] phi [Ekin gradphi potEnergy]
        data = np.loadtxt(path, ndmin=2)
        axes = [np.unique(data[:, d]) for d in range(D)]
        shape = tuple(len(a) for a in axes)
        fields["phi"] = data[:, D].reshape(shape)
        if data.shape[1] > D + 2:
            fields["gradphi"] = data[:, D + 2].reshape(shape)

    return fields, axes, info


def read_spectrum(path):
    """Read a power spectrum file _phi_k_t=NNNN.dat (simulation units).

    Returns t, a, k_comoving, k_physical, Delta(k). Older files have only the
    columns k_physical and Delta; then k_comoving = a * k_physical.
    """
    t, a = np.loadtxt(path, skiprows=1, max_rows=1)
    data = np.loadtxt(path, skiprows=2, ndmin=2)
    if data.shape[1] >= 3:
        k_com, k_phys, delta = data[:, 0], data[:, 1], data[:, 2]
    else:
        k_phys, delta = data[:, 0], data[:, 1]
        k_com = a * k_phys
    return t, a, k_com, k_phys, delta
