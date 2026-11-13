/*==============================================================================
 *  CoolBubble -- lattice simulation of vacuum bubble nucleation, growth and
 *  collision during cosmological first-order phase transitions.
 *
 *  File: bubble.hpp
 *  Role: Declaration of the Bubble struct -- simulation state (field, grid,
 *        expansion, wall-tracking) and the methods implemented in bubble.cpp.
 *
 *  See README.md for build/usage instructions and def.h for run parameters.
 *
 *  Copyright (C) 2026 Henda Mansour, Yann Gouttenoire
 *  SPDX-License-Identifier: GPL-3.0-or-later
 *  This program is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License, version 3 or (at your
 *  option) any later version. See the LICENSE file for details.
 *==============================================================================*/

#ifndef BUBBLE_HPP
#define BUBBLE_HPP

#include <cmath>
#include <fstream>
#include <functional>
#include <string>
#include <vector>
#include <fftw3.h>

#include "def.h"
#include "gridfuncs.h"

// zero-padded snapshot number used in output file names, e.g. 7 -> "0007"
std::string dataName(int itPlot);

struct Bubble
{
    Bubble(std::string folder);

    std::string folder;             // output folder of this run

    //------------------------------
    // Lattice
    //------------------------------
    int N = std::pow(nGrids, D);                 // total number of grid points
    int n_cplx = nGrids / 2 + 1;                 // complex modes along the last axis of a real FFT
    double dx = (double)L / (double)nGrids;      // lattice spacing
    double dxD  = std::pow(dx, D);               // volume element dx^D
    double dxLD = std::pow(dx / L, D);           // volume element as a fraction of the box, (dx/L)^D
    double dxL = dx / L;
    int nkmax = (int)(std::sqrt(double(D)) * nGrids / 2); // number of |k| bins in the power spectrum

    //------------------------------
    // Momentum scales (written to parameters.dat)
    //------------------------------
    double kUV = std::sqrt(double(D)) * M_PI / dx; // largest momentum on the lattice
    double kIR = 2 * M_PI / double(L);            // smallest momentum in the box
    double kRStar;                                // momentum scale of the bubbles at collision (estimate)
    double kMeff;                                 // mass in the true vacuum

    //------------------------------
    // Potential and critical bubble
    //------------------------------
    double phi0, DeltaV, phiMax, mT; // true vacuum, latent heat, top of the barrier, mass in the true vacuum
    double phiZero;                  // exit point: V(phiZero) = V(0)
    double Rn, lw, phiC;             // radius, wall width and central field value of a nucleated bubble

    //------------------------------
    // Fields
    //------------------------------
    double *phi, *pi, *dpi;              // field, canonical momentum pi = a^(D-ALPHAVAL) phi', and the force on pi
    double *gradphi, *potEnergy, *Ekin;  // local energy densities (only allocated if memorysave == 0)

    //------------------------------
    // Background expansion
    //------------------------------
    double a = 1., ap = 0., H = 0.; // scale factor, its time derivative, Hubble rate
    double C = 0.;                  // Hubble constraint H^2 - rho/(3 Mpl^2), should stay small
    bool expanding() const { return expansion != 0 && D == 3; } // expansion is only implemented in 3D

    //------------------------------
    // Box averages, updated every step
    //------------------------------
    double t = 0., avphi = 0, avphiP = 0;
    double Etot = 0, Ekintot = 0, Epottot = 0, Egradtot = 0;
    double EkinHalfPrev = 0;        // kinetic energy at the previous half step (for the time-centred Ekintot)

    //------------------------------
    // Wall tracking
    //------------------------------
    double dini = 0, WallWidth = 0, Rmin = 0, Rmax = 0; // initial and current wall width, wall radii
    std::vector<int> wallCentre;    // lattice point at the centre of the tracked bubble

    bool CanTrackWall() const;      // wall tracking needs an isolated bubble: 1D, or mode -1 / -2 in 2D and 3D
    void TrackWall(int it);         // wall tracking at time step it
    int WallTrackStart() const;     // time step at which wall tracking starts
    void UpdateWall(int it);        // measure the wall now, without writing to gamma.dat (for the snapshot plots)
    void WallFinder();              // finds the wall of the tracked bubble: updates Rmin, Rmax and WallWidth
    double GammaFactor(double d);   // Lorentz factor of the wall, dini / d for wall width d

    //------------------------------
    // Initial conditions
    //------------------------------
    void Initialize1D();                // 1D: bubble on the box edge
    void InitializeWave();              // TEST = 1: sine wave around the true vacuum
    void NucleateOneBubble();           // mode -2
    void NucleateTwoDiagonalBubbles();  // mode -1
    void GenerateNucleationPoints();    // mode 1: random positions -> __xns.dat
    void NucleateBubblesAtGivenPoints();// modes 0 and 1: bubbles at the positions in __xns.dat
    void GenerateThermalFluc(double *field, double T, double m, std::function<double(double,double,double)> powerspectrum,
                             const std::string &specname, double shift = 0.0); // mode -3

    //------------------------------
    // Time evolution
    //------------------------------
    void oneStep();                 // one leapfrog step
    void computeForces();           // force on pi -> dpi
    double addot() const;           // a'' from the Friedmann equations
    void measureEnergies();         // box-averaged energies and field averages

    //------------------------------
    // Power spectrum
    //------------------------------
    void FourierTransform(int itPlot);   // power spectrum of snapshot itPlot -> _phi_k_t=NNNN.dat
    void ComputePowerSpectrum(fftw_complex *phi_k, std::fstream &dat);

    //------------------------------
    // Output
    //------------------------------
    void InitOutputFiles();         // headers of EnergyDensities.dat and gamma.dat
    void OutParameters();           // parameters.dat
    void outEP();                   // one line of EnergyDensities.dat
    void outGamma();                // one line of gamma.dat
    void outFullField();            // 1D: one line of FullField.dat
    void OutField(int itPlot);      // field snapshot _phi_t=NNNN.h5
    void CreateVideoIfFound(const std::string& type); // video of the snapshot plots (ffmpeg)
};

#endif // BUBBLE_HPP


