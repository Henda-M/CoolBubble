/*==============================================================================
 *  CoolBubble -- lattice simulation of vacuum bubble nucleation, growth and
 *  collision during cosmological first-order phase transitions.
 *
 *  File: def.h
 *  Role: Compile-time simulation parameters. Values ending in _default can
 *        also be changed at run time with command-line options (see
 *        ./coolbubble --help); everything else is fixed when compiling, so
 *        rebuild with `make` after editing.
 *
 *  Units: all quantities are dimensionless, in units of v_phi (the vacuum
 *  expectation value in the true minimum), so phi = 1 in the true vacuum and
 *  the latent heat is DeltaV = v_phi^(D+1) = 1. The plots use units of the
 *  true-vacuum mass m_phi = sqrt(2(c+6)) v_phi.
 *
 *  Copyright (C) 2026 Henda Mansour, Yann Gouttenoire
 *  SPDX-License-Identifier: GPL-3.0-or-later
 *  This program is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License, version 3 or (at your
 *  option) any later version. See the LICENSE file for details.
 *==============================================================================*/

#ifndef DEF_H
#define DEF_H

/*---------------------------*/
// Lattice           (run-time options: --D, --nGrids, --L, --nt)
/*---------------------------*/
#define D_default 2          // number of spatial dimensions: 1, 2 or 3
#define nGrids_default 512  // grid points per dimension; lattice spacing dx = L/nGrids
#define L_default 32         // box size
#define nt_default 1000      // number of time steps
#define TIME_STEP 0.004      // time step dt; must satisfy dt < dx/sqrt(D)

/*---------------------------*/
// Potential and initial conditions   (run-time options: --c, --mode)
/*---------------------------*/
#define c_default 4.     // potential parameter c of V(phi) = c phi^2 - (2c+4) phi^3 + (c+3) phi^4 + 1 (see bubble.cpp)
#define readxns_default 1 // initial condition for 2D and 3D:
                          //   1: nBubbles bubbles at random positions
                          //   0: bubbles at the positions read from <run folder>/__xns.dat
                          //  -1: two bubbles on the box diagonal
                          //  -2: one bubble in the centre
                          //  -3: thermal fluctuations instead of bubbles
                          // In 1D the run always starts from two colliding bubbles.
#define nBubbles 5        // number of bubbles for mode 1
#define fixedseed 1       // 1: fixed random seed (reproducible runs), 0: seed from the clock
#define Trad 0.2          // mode -3: temperature of the fluctuations, T = Trad * m (m = mass in the true vacuum);
                          //   the fluctuations are around the false vacuum, with mass m_f = sqrt(2c)
#define DeltaR 0.         // mode -1: shift of both bubbles towards each other along the diagonal
#define TEST 0            // 1: start from a small sine wave around the true vacuum instead (for testing)

/*---------------------------*/
// Expansion
/*---------------------------*/
#define expansion 0  // 0: static background, 1: self-consistent expansion (3D only; 1D and 2D are always static)
#define ALPHAVAL 0   // time variable: 0 = cosmic time, 1 = conformal time
#define Mpl 200.0    // reduced Planck mass; smaller values mean stronger gravity

/*---------------------------*/
// Output             (run-time options: --out, --snapshots)
/*---------------------------*/
#define Path_default "./"   // folder in which the run folder is created
#define memorysave 0        // 0: also save the local kinetic, gradient and potential energy densities in the snapshots
#define OutEnergyEvery 40   // write the box-averaged energies every this many steps
#define nSnapshots 20       // number of field snapshots (run-time option --snapshots); the initial state is saved in addition
#define SnapshotSpacing 0   // 0: snapshots evenly spaced in time (best for videos)
                           // 1: denser until shortly after the first collisions: half of the snapshots evenly
                           //    up to t = L, the other half evenly from t = L to the end of the run
#define FFT 1               // 1: compute the power spectrum at every snapshot, 0: skip it
#define Full1DFieldEvery 100 // 1D only: save the full field in FullField.dat every this many steps
#define WallTracking 0      // 1: track the wall width and Lorentz factor of a bubble (reliable in 1D only)
#define WallTrackEvery 400  // wall tracking: update the wall position every this many steps
#define plotsnapshots 1     // 1: plot every snapshot (field and power spectrum) and make a video of them during the run,
                           // 0: only save the snapshot files; the summary plots at the end are made in any case
#define makevideos 0        // 1: make videos of the snapshots (kinetic and gradient energy) at the end of the run, 0: skip it

#endif // DEF_H
