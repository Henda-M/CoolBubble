/*==============================================================================
 *  CoolBubble -- lattice simulation of vacuum bubble nucleation, growth and
 *  collision during cosmological first-order phase transitions.
 *
 *  File: bubble.cpp
 *  Role: Implementation of the Bubble class -- initial conditions, time
 *        evolution, energy measurement, power spectrum, wall tracking and
 *        file output.
 *
 *  See README.md for build/usage instructions and def.h for run parameters.
 *
 *  Copyright (C) 2026 Henda Mansour, Yann Gouttenoire
 *  SPDX-License-Identifier: GPL-3.0-or-later
 *  This program is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License, version 3 or (at your
 *  option) any later version. See the LICENSE file for details.
 *==============================================================================*/

#include <cstdio>
#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <omp.h>
#include <time.h>
#include <fftw3.h>
#include <filesystem>
#include <cstdlib>
#include <functional>

#include "H5Cpp.h"
#include "bubble.hpp"
#include <stdexcept>


namespace fs = std::filesystem;
using namespace H5;
using namespace std;

// spread of the box averages (square root of the sum of the squared cell contributions), written to EnergyDensities.dat
static double EtotVar, EkintotVar, EpottotVar, EgradtotVar, avphiPVar, avphiVar;

static int nthreads = omp_get_max_threads();

//------------------------------
// Helpers
//------------------------------

// zero-padded snapshot number used in output file names, e.g. 7 -> "0007"
string dataName(int itPlot)
{
    char buf[16];
    snprintf(buf, sizeof(buf), "%04d", itPlot);
    return string(buf);
}

//------------------------------
// Potential and critical bubble
//------------------------------

// Potential parametrisation of Jinno and Konstandin (2017). All quantities are in units of v_phi,
// the vacuum expectation value in the true minimum:
//   V(phi) = v^(d+1) [c phi^2 - (2c+4) phi^3 + (c+3) phi^4] + v^(d+1),   with v_phi = 1.
// True vacuum at phi = 1 with V = 0, false vacuum at phi = 0 with V = DeltaV = v^(d+1) = 1 (the latent heat),
// barrier top at phiMax. Higher cc gives a higher barrier --> thinner walls.
double Veps(double phi, double cc, double DeltaV)
{
    return cc * pow(phi, 2) - (2 * cc + 4) * pow(phi, 3) + (cc + 3) * pow(phi, 4) + DeltaV;
}
// dV/dphi
double Vepspr(double phi, double cc)
{
    return 2 * cc * pow(phi, 1) - 3 * (2 * cc + 4) * pow(phi, 2) + 4 * (cc + 3) * pow(phi, 3);
}

// Fits to the critical bubble as functions of c, used to set up the nucleated bubbles.
// field value at the centre
double phizero(double c)
{
    return 0.997 * (1 - exp(-c * 0.688));
}

// radius; in 3D the bubbles are nucleated 1.3 times larger
double CRadius(double c)
{
    if( D == 3 )
    {
        return 0.823 * pow(c, 0.203)*1.3; 
    }
    return 0.823 * pow(c, 0.203);
}
// wall width
double Clw(double c)
{
    return 0.496 * pow(c, -0.531);
}

// tanh profile of a nucleated bubble, at squared distance x2 from its centre
double ThickBubble(double x2, double c)
{
    double phiC = phizero(c);
    double Rn = CRadius(c);
    double lw = Clw(c);
    return (phiC / 2.) * (tanh((Rn - sqrt(x2)) / lw) + 1.);
}

//------------------------------
// Thermal fluctuations
//------------------------------
// Thermal power spectra of a free field with mass m at temperature T (Bose-Einstein occupation n_k,
// vacuum fluctuations not included): <|phi_k|^2> = n_k / omega_k and <|phi'_k|^2> = n_k omega_k,
// with omega_k = sqrt(k^2 + m^2), in the continuum normalisation used by GenerateThermalFluc.
double thermalspectrum(double T, double k2, double m)
{
    double omega_k = std::sqrt(k2 + m * m);
    double n_k = 1.0 / (std::exp(omega_k / T) - 1.0);
    double power = (n_k) / omega_k;     
    return power; 
}

double velocitythermalspectrum(double T, double k2, double m)
{
    double omega_k = std::sqrt(k2+ m* m);
    double n_k = 1.0 / (std::exp(omega_k / T) - 1.0);
    double power = (n_k) * omega_k;     
    return power; 
}

// Base seed for the thermal fluctuations. phi and its velocity get different
// streams (offset by the spectrum name), so they are statistically independent.
// Follows `fixedseed` in def.h: fixed (reproducible) or taken from the clock.
static uint64_t thermalSeed(const std::string &specname)
{
    static const uint64_t base = fixedseed ? 12345ULL : (uint64_t)time(NULL);
    uint64_t offset = 0;
    for (char ch : specname)
        offset = offset * 131 + (unsigned char)ch;
    return base + offset;
}

// test initial condition (TEST = 1): sine wave
double Wave(double x, double k, double amplitude)
{
    return amplitude * sin(k * x);
}


//------------------------------
// Constructor
//------------------------------

// Sets up the lattice, the initial condition, the initial energies and the expansion, and the output files.
Bubble::Bubble( std::string folder)
    : folder(folder)
{
    t = 0.;
    // potential and critical bubble
    DeltaV = 1.;                                                // latent heat v_phi^(D+1), = 1 in units of v_phi
    phi0 = 1.;                                                  // true vacuum
    phiC = phizero(c);                                          // field at the centre of a nucleated bubble
    Rn = CRadius(c);                                            // nucleation radius
    lw = Clw(c);                                                // wall width of critical bubble
    kRStar = sqrt(double(D)) * 4 * M_PI / (double(L - 2 * Rn)); // momentum scale of the bubbles at collision (estimate)
    phiMax = c / (6 + 2 * c);                  // top of the barrier
    phiZero = (2 + c - sqrt(4 + c)) / (3 + c); // exit point: V(phiZero) = V(0)
    mT = sqrt(2 * (c + 6));                    // mass in true vacuum
    kMeff = (double(mT));                      // same, as momentum scale for the plots

    // fields, all set to zero
    phi = (double *)calloc(N, sizeof(double)); // field
    pi = (double *)calloc(N, sizeof(double));  // canonical momentum pi = a^(D-ALPHAVAL) phi'
    dpi = (double *)calloc(N, sizeof(double)); // its time derivative, the force (see computeForces)
    // local energy densities, saved in the snapshots
    if (memorysave == 0) 
    {
        gradphi = (double *)calloc(N, sizeof(double));
        potEnergy = (double *)calloc(N, sizeof(double));
        Ekin = (double *)calloc(N, sizeof(double));
    }

    // initial condition: start from the false vacuum phi = 0, then add bubbles or fluctuations
#pragma omp parallel for num_threads(nthreads)
    for (int iGrid = 0; iGrid < N; iGrid++)
    {
        phi[iGrid] = 0.;
    }

    if (D == 1) // 1D: always two colliding bubbles
    {
        Initialize1D();
    }
    else if (D == 2 || D == 3)
    {
        if (TEST ==1) // test: sine wave around the true vacuum
        {
            InitializeWave();
        }
        else if (readxns == -2) // one bubble in the centre
        {
            NucleateOneBubble();
        }

        else if (readxns == -1) // two bubbles on the diagonal
        {
            NucleateTwoDiagonalBubbles();
        }
        else if (readxns == 1) // bubbles at random positions
        {
            GenerateNucleationPoints();
            NucleateBubblesAtGivenPoints();
        }
        else if (readxns == 0) // bubbles at the positions listed in __xns.dat
        {
            NucleateBubblesAtGivenPoints();
        }
        else if (readxns == -3) // thermal fluctuations
        {
            // fluctuations around the false vacuum phi = 0, where the mass is m_f = sqrt(V''(0)) = sqrt(2c);
            // the temperature is given in units of the true-vacuum mass (Trad in def.h)
            double Temp = Trad * mT;
            double mF = sqrt(2 * c);
            cout << "Generating thermal fluctuations around the false vacuum with T = " << Temp << ", m_f = " << mF << endl;
            GenerateThermalFluc(phi, Temp, mF, thermalspectrum, "thermalspec");
            GenerateThermalFluc(pi, Temp, mF, velocitythermalspectrum, "velocitythermalspec");
        } 
    }
    else 
    {
        cout << "Error: D must be 1, 2 or 3!" << endl;
        exit(1);
    }

// initial energies and averages (a = 1 at the start). The kinetic energy is zero for bubbles,
// which start at rest, but not for thermal fluctuations, where pi is drawn as well.
// The local energy densities (memorysave = 0) use the same definitions as measureEnergies().
#pragma omp parallel for reduction(+ : Epottot, Egradtot, Ekintot, avphi, avphiP)
    for (int iGrid = 0; iGrid < N; iGrid++)
    {
        double tmp = 0;
        if (D == 1)
        {
            tmp += pow((phi[iGridMFunc(iGrid,1)] - phi[iGrid]) / dx, 2);
        }
        else if (D == 2 || D == 3)
        {
            // forward differences: consistent with the compact Laplacian in computeForces
            for (int iShift = 0; iShift < D; iShift++)
                tmp += pow((phi[iGridPFunc(iGrid,iShift)] - phi[iGrid]) / dx, 2);
        }

        double egrad = 0.5 * tmp;
        double epot = Veps(phi[iGrid], c, DeltaV);
        double ekin = 0.5 * pi[iGrid] * pi[iGrid];
        if (memorysave == 0)
        {
            gradphi[iGrid] = egrad;
            potEnergy[iGrid] = epot;
            Ekin[iGrid] = ekin;
        }
        avphi += phi[iGrid] * dxLD;
        avphiP += pi[iGrid] * dxLD;
        Epottot += epot * dxLD;
        Egradtot += egrad * dxLD;
        Ekintot += ekin * dxLD;
    }
    Etot = Ekintot + Epottot + Egradtot;

    // expansion: a = 1 at the start, initial Hubble rate from the Friedmann equation
    // (H0 = 0 without expansion, so that a stays exactly 1)
    a = 1;
    double H0 = expanding() ? sqrt(pow(a, 2*ALPHAVAL) * Etot / (3.0 * Mpl * Mpl)) : 0.0;
    ap = a * H0;
    ap -= 0.5 * TIME_STEP * addot();   // a' lives at half steps: step it back to t = -dt/2

    // Leapfrog: phi lives at whole steps, pi at half steps. Step pi back from t = 0 to t = -dt/2,
    // pi_{-1/2} = pi_0 - dt/2 F(phi_0), with the force F(phi_0) at t = 0:
    computeForces();   // dpi = F(phi_0)
    double a_mhalf = a - 0.5 * TIME_STEP * ap; // a at t = -1/2
    double EkinHalf = 0.0;
#pragma omp parallel for reduction(+ : EkinHalf)
    for (int iGrid = 0; iGrid < N; iGrid++)
    {
        pi[iGrid] -= 0.5 * TIME_STEP * dpi[iGrid];
        EkinHalf += 0.5 * pow(a_mhalf, -2*D) * pi[iGrid] * pi[iGrid] * dxLD;
    }
    EkinHalfPrev = EkinHalf; // kinetic energy at t = -1/2

    H = H0; // Hubble rate at t = 0

    cout << "IC: a = " << a << " ap = " << ap << " ap/a = " << ap / a << endl;
    if (expanding())
    {
        cout << "Hubble Length 1/H = " << (1.0 / H) << endl; 
        if ((1.0 / H) < 10*L)
        {
            cout << "Warning: H^-1 < 10 L " << endl;
        }
    }
    else
    {
        cout << "No expansion: a = 1" << endl;
    }
    cout << "dxLD: " << dxLD << " dxL: " << dxL << endl;


    // centre of the bubble followed by the wall tracking
    int mid = (D == 3) ? nGrids / 2 : 0;                                         // third index only exists in 3D
    if (readxns == -2 && D >= 2)      wallCentre = {nGrids / 2, nGrids / 2, mid}; // bubble in the centre
    else if (readxns == -1 && D >= 2) wallCentre = {nGrids / 4, nGrids / 4, mid}; // first diagonal bubble (DeltaR = 0)
    else                              wallCentre = {0, 0, 0};                     // 1D: bubble at the box edge

    InitOutputFiles(); // headers of EnergyDensities.dat and gamma.dat
}

//------------------------------
// Initial conditions
//------------------------------

// 1D: a bubble centred on the box edge (profiles around x = 0 and x = L). In the periodic box its two
// walls meet again at x = L/2, as for two bubbles a distance L apart.
void Bubble::Initialize1D()
{
#pragma omp parallel for num_threads(nthreads)
    for (int iGrid = 0; iGrid < (int)nGrids; iGrid++)
    {
        double x = xFunc(iGrid, L)[0];  // position of the lattice point
        double xxn = x;
        phi[iGrid] = min(phi[iGrid] + (ThickBubble(xxn * xxn, c)), phiC); // add the profile, capped at phiC
        xxn = x - L;
        phi[iGrid] = min(phi[iGrid] + (ThickBubble(xxn * xxn, c)), phiC);
    }
}

// Test initial condition (TEST = 1, 2D and 3D): a small sine wave around the true vacuum
void Bubble::InitializeWave()
{
    double k = 2.0 * M_PI / (L / 40.0); // wave number: 40 wavelengths per box
    double amplitude = 0.01 * phi0;     // 1% of the true-vacuum value
    #pragma omp parallel for num_threads(nthreads)
    for (int iGrid = 0; iGrid < N; iGrid++)
    {
        vector<double> x = xFunc(iGrid, L); // position of the lattice point
        double wavevalue = 0.0;
        for (int d = 0; d < D; d++)
        {
            wavevalue += Wave(x[d], k, amplitude);
        }               
        phi[iGrid] = phi0 + wavevalue / double(D); // average over dimensions
    }
}

// Mode -2: one bubble in the centre of the box
void Bubble::NucleateOneBubble()
{
#pragma omp parallel for num_threads(nthreads)
    for (int iGrid = 0; iGrid < N; iGrid++)
    {
        vector<double> x = xFunc(iGrid, L);                                 // position of the lattice point
        vector<double> xxn = {x[0] - L / 2., x[1] - L / 2., x[2] - L / 2.}; // distance to the box centre
        double rn = (D == 2)
        ? (xxn[0]*xxn[0] + xxn[1]*xxn[1])
        : (xxn[0]*xxn[0] + xxn[1]*xxn[1] + xxn[2]*xxn[2]);

        phi[iGrid] = min(phi[iGrid] + (ThickBubble(rn, c)), phiC); // add the profile, capped at phiC
    }
}

// Mode -1: two bubbles on the diagonal, at (L/4, L/4) and (3L/4, 3L/4) (z = L/2 in 3D),
// both moved towards each other by DeltaR along x and y
void Bubble::NucleateTwoDiagonalBubbles()
{
#pragma omp parallel for num_threads(nthreads)
    for (int iGrid = 0; iGrid < N; iGrid++)
    {
        vector<double> x = xFunc(iGrid, L); // position of the lattice point
        vector<double> xxn1 = {x[0] - L / 4 - DeltaR, x[1] - L / 4. - DeltaR, x[2] - L / 2.};           // to the first bubble
        vector<double> xxn2 = {x[0] - 3. * L / 4. + DeltaR, x[1] - 3. * L / 4 + DeltaR, x[2] - L / 2.}; // to the second bubble

        double rn1 = (D == 2)
        ? (xxn1[0]*xxn1[0] + xxn1[1]*xxn1[1])
        : (xxn1[0]*xxn1[0] + xxn1[1]*xxn1[1] + xxn1[2]*xxn1[2]);

        double rn2 = (D == 2)
        ? (xxn2[0]*xxn2[0] + xxn2[1]*xxn2[1])
        : (xxn2[0]*xxn2[0] + xxn2[1]*xxn2[1] + xxn2[2]*xxn2[2]);

        phi[iGrid] = min(phi[iGrid] + (ThickBubble(rn1, c)), phiC); // add the profiles, capped at phiC
        phi[iGrid] = min(phi[iGrid] + (ThickBubble(rn2, c)), phiC);
    }
}

// Mode 1: draws nBubbles random positions and writes them to __xns.dat, together with the periodic
// images of bubbles that cross the box edge (fixed or clock-based seed, see fixedseed in def.h)
void Bubble::GenerateNucleationPoints()
{
    if (fixedseed)
    {
        srand(42);
    }
    else
    {
        srand((unsigned)time(NULL));
    }

    FILE *fp;
    fp = fopen((folder + "/__xns.dat").c_str(), "w");
    for (int iBubble = 0; iBubble < nBubbles; iBubble++)
    {
        vector<double> xn = {L * (double)rand() / (double)RAND_MAX, L * (double)rand() / (double)RAND_MAX, L * (double)rand() / (double)RAND_MAX};
        fprintf(fp, "%.5e %.5e %.5e\n", xn[0], xn[1], xn[2]);
        check_PBC(xn, L, Rn + 5 * lw, fp); // periodic images, also for bubbles whose wall tail crosses the box edge
    }
    fclose(fp);
}

// Modes 0 and 1: one bubble at each position listed in <run folder>/__xns.dat (one "x y z" per line)
void Bubble::NucleateBubblesAtGivenPoints()
{
    FILE *fp;
    fp = fopen((folder + "/__xns.dat").c_str(), "r");
    if (!fp)
    {
        cerr << "Error: mode 0 needs the bubble positions in " << folder << "/__xns.dat (one \"x y z\" per line)." << endl;
        exit(1);
    }
    vector<double> xn = {0, 0, 0};
    int nRead = 0;
    while (fscanf(fp, "%lf %lf %lf", &xn[0], &xn[1], &xn[2]) == 3) // stops at the end of the file or at a malformed line
    {
        nRead++;
        #pragma omp parallel for num_threads(nthreads)
        for (int iGrid = 0; iGrid < N; iGrid++)
        {
            vector<double> x = xFunc(iGrid, L);
            vector<double> xxn = {x[0] - xn[0], x[1] - xn[1], x[2] - xn[2]}; // distance to the bubble centre

            if (D == 2)
            {
                phi[iGrid] = min((phi[iGrid] + ThickBubble(xxn[0] * xxn[0] + xxn[1] * xxn[1], c)), phiC);
            }
            else if (D == 3)
            {
                phi[iGrid] = min((phi[iGrid] + ThickBubble(xxn[0] * xxn[0] + xxn[1] * xxn[1] + xxn[2] * xxn[2], c)), phiC);
            }
        }
    }
    fclose(fp);
    if (nRead == 0)
        cerr << "Warning: no bubble positions found in " << folder << "/__xns.dat" << endl;
}

// Gaussian random field for the thermal initial conditions (2D and 3D): either the field or its velocity.
// Each Fourier mode is drawn with <phi_k phi_k'^*> = (2 pi)^D delta(k - k') P(k), P(k) = powerspectrum(T, k^2, m)
// (continuum normalisation: the power per mode is the same for any lattice spacing; real-space quantities such
// as <phi^2> still depend on dx through the number of UV modes, k up to ~ pi/dx). `shift` is added in real space.
void Bubble::GenerateThermalFluc(double *field, double T, double m,
                                 std::function<double(double, double, double)> powerspectrum,
                                 const std::string &specname, double shift)
{
    // half-complex layout of FFTW: nRows rows (all values of the leading D-1 momenta) of n_cplx modes each
    int nRows = (D == 3) ? nGrids * nGrids : nGrids;
    fftw_complex *field_k = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * nRows * n_cplx);
    int dims[3] = {nGrids, nGrids, nGrids};
    fftw_plan plan = fftw_plan_dft_c2r(D, dims, field_k, field, FFTW_ESTIMATE);

    auto kphys = [](int i) { return 2.0 * M_PI * ((i <= nGrids / 2) ? i : i - nGrids) / L; };

    // each row gets its own random stream, seeded with (seed + row),
    // so the result does not depend on the number of threads
    #pragma omp parallel for num_threads(nthreads)
    for (int r = 0; r < nRows; ++r)
    {
        std::mt19937_64 gen(thermalSeed(specname) + r);
        std::normal_distribution<double> normal(0.0, 1.0);

        double k2rows = (D == 2) ? pow(kphys(r), 2)
                                 : pow(kphys(r / nGrids), 2) + pow(kphys(r % nGrids), 2);
        for (int j = 0; j < n_cplx; ++j)
        {
            double kj = 2.0 * M_PI * j / L;
            double sigma = std::sqrt(powerspectrum(T, k2rows + kj * kj, m) / 2.); // per real and imaginary part
            field_k[r * n_cplx + j][0] = sigma * normal(gen);
            field_k[r * n_cplx + j][1] = sigma * normal(gen);
        }
    }

    // A real field needs phi_{-k} = phi_k^*. For 0 < j < nGrids/2 FFTW only stores one of the two, but on the
    // planes j = 0 and j = nGrids/2 both partners are stored: copy the conjugate of the first one into the second.
    // Modes that are their own partner must be real, with the full variance P(k) in the real part.
    for (int j : {0, nGrids / 2})
    {
        for (int r = 0; r < nRows; ++r)
        {
            int partner = (D == 2) ? (nGrids - r) % nGrids
                                   : ((nGrids - r / nGrids) % nGrids) * nGrids + (nGrids - r % nGrids) % nGrids;
            double *mode = field_k[r * n_cplx + j];
            if (partner == r)
            {
                mode[0] *= std::sqrt(2.);
                mode[1] = 0.;
            }
            else if (partner < r)
            {
                mode[0] = field_k[partner * n_cplx + j][0];
                mode[1] = -field_k[partner * n_cplx + j][1];
            }
        }
    }

    // inverse FFT to real space. FFTW does not normalise; dividing by sqrt(V) = sqrt(L^D) gives the
    // continuum normalisation <|phi_k|^2> = V P(k) with phi_k = dx^D sum_x phi(x) e^{-ikx}
    fftw_execute(plan);
    double norm = 1.0 / std::sqrt(pow(L, D));
    #pragma omp parallel for num_threads(nthreads)
    for (int i = 0; i < N; ++i)
        field[i] = field[i] * norm + shift;

    fftw_destroy_plan(plan);
    fftw_free(field_k);
}

//------------------------------
// Wall tracking
//------------------------------

// Wall tracking (WallTracking = 1): the field is read along one line, from the bubble centre outwards along
// the first lattice axis. Rmin is the outermost radius where the field is still above its exit value phiZero,
// Rmax the first radius where it has dropped to the barrier top phiMax. The wall width is Rmax - Rmin, and the
// Lorentz factor of the wall is its initial width divided by the current one (see Fig. 9 in arXiv:1906.02588).
// Needs an isolated bubble: used in 1D, and in 2D/3D for a single bubble (mode -2) or the first of two
// diagonal bubbles (mode -1).
bool Bubble::CanTrackWall() const
{
    return D == 1 || readxns == -1 || readxns == -2;
}

// Tracking starts after half an oscillation period of the field in the true vacuum, when the bubble profile
// has settled; the wall width at that time is the reference for the Lorentz factor. Afterwards the wall is
// updated every WallTrackEvery steps until t = L/2.
void Bubble::TrackWall(int it)
{
    int itStart = WallTrackStart();
    if (it == itStart)
    {
        WallFinder();
        dini = WallWidth;
        cout << "initial wall width: " << dini << endl;
        outGamma();
    }
    else if (it > itStart && it % WallTrackEvery == 0 && t < L / 2)
    {
        WallFinder();
        outGamma();
        cout << "gamma factor: " << GammaFactor(WallWidth) << endl;
    }
}

// Time step at which wall tracking starts: half an oscillation period of the field in the true vacuum
int Bubble::WallTrackStart() const
{
    return int(M_PI / mT / TIME_STEP);
}

// Measures the wall at the current step without writing to gamma.dat, so that the markers in the
// snapshot plots show the current wall position (once tracking has started, until t = L/2)
void Bubble::UpdateWall(int it)
{
    if (it > WallTrackStart() && t < L / 2)
        WallFinder();
}

void Bubble::WallFinder()
{
    const vector<int> &centre = wallCentre;
    // radius at which the bubble meets another bubble or its own periodic image
    double MaximalBubbleRadius = (readxns == -1 && D >= 2) ? L / 4.0 - DeltaR : L / 2.0;
    // bubbles have already collided: no tracking needed anymore -> exit
    if (Rmax >= 0.98 * MaximalBubbleRadius || Rmin >= 0.98 * MaximalBubbleRadius)
    {
        cout << "bubble has almost left the box" << endl;
        return;
    }

    // field and distance from the centre at step i along the line
    auto phiAlong = [&](int i) {
        if (D == 1) return phi[i];
        int i0 = ((centre[0] + i) % nGrids + nGrids) % nGrids;
        return phi[ijkTOiGrids({i0, centre[1], centre[2]})];
    };
    // distance from the bubble centre: the centres (x = 0 in 1D, L/2 or L/4 in 2D/3D) lie on cell corners,
    // while the field values sit at the cell centres, half a lattice spacing away along every axis
    auto radius = [&](int i) { return dx * sqrt((i + 0.5) * (i + 0.5) + 0.25 * (D - 1)); };

    // start looking just inside the previous wall position (or inside the nucleated bubble)
    int istart = max(int((Rmin - lw) / dx), int(Rn / (2 * dx)));
    Rmin = 0;
    Rmax = 0;
    for (int i = istart; (Rmax == 0 || Rmin == 0) && radius(i) < MaximalBubbleRadius; i++)
    {
        if (phiAlong(i) > phiZero)                  // still inside: above the exit value of the field
            Rmin = radius(i);
        if (phiAlong(i) <= phiMax && Rmax == 0)     // outside: at or below the top of the barrier
            Rmax = radius(i);
    }

    if (Rmax < Rmin)
    {
        cout << "Rmax < Rmin ! Something went wrong ... " << endl;
        return;
    }
    if (Rmax * Rmin == 0 && Rmin + Rmax != 0)
        Rmax = Rmin;

    WallWidth = Rmax - Rmin;
    cout << "WallWidth is: " << WallWidth << endl;
}

// Lorentz factor of the wall: initial wall width / current width d (0 if not known yet)
double Bubble::GammaFactor(double d)
{
    if (WallTracking == 0 )
    {
        return 0;
    }

    if (d == 0)
    {
        return 0;
    }
    return dini / d;
}

//------------------------------
// Power spectrum
//------------------------------

// Power spectrum of the field at snapshot itPlot, written to _phi_k_t=NNNN.dat: a header, the time t and
// scale factor a, then one line per |k| bin with the comoving and physical momentum and Delta(k).
void Bubble::FourierTransform(int itPlot)
{
    if (FFT != 1)
        return;

    // real-to-complex FFT: N/nGrids rows of n_cplx complex modes (FFTW keeps phi unchanged)
    fftw_complex *phi_k = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * (N / nGrids) * n_cplx);
    int dims[3] = {nGrids, nGrids, nGrids};
    fftw_plan p = fftw_plan_dft_r2c(D, dims, phi, phi_k, FFTW_ESTIMATE);
    fftw_execute(p);

    fstream dat((folder + "/_phi_k_t=" + dataName(itPlot) + ".dat").c_str(), ios::out);
    dat << "# t    a        (then one line per |k| bin: k_comoving  k_physical  Delta(k))" << endl;
    dat << t << " " << a << endl;
    ComputePowerSpectrum(phi_k, dat);
    dat.close();

    fftw_destroy_plan(p);
    fftw_free(phi_k);
}

// Dimensionless power spectrum Delta(k), defined by <phi^2> = int Delta(k) dln k:
//   Delta(k) = k^D P(k) / (2^(D-1) pi^(D/2) Gamma(D/2))   (1D: k P / pi,  2D: k^2 P / 2pi,  3D: k^3 P / 2pi^2)
// with P(k) = dx^2D |phi_k|^2 / L^D. Both k and P are comoving here, so Delta is the same in comoving and
// physical coordinates. The modes are grouped in shells of width 2 pi / L in |k|; each output line gives the
// mean comoving momentum |k| of the modes in the shell, the physical momentum |k|/a, and their mean Delta.
void Bubble::ComputePowerSpectrum(fftw_complex *phi_k, std::fstream &dat)
{
    const int nRows = N / nGrids;     // rows of the half-complex array (1 in 1D)
    const int nBins = nkmax + 1;
    const double dk = 2.0 * M_PI / L; // comoving momentum spacing
    const double prefactor = (D == 1) ? 1.0 / M_PI : (D == 2) ? 1.0 / (2.0 * M_PI) : 1.0 / (2.0 * M_PI * M_PI);
    const double Pnorm = dxD * dxD / pow(L, D);
    auto signedIndex = [](int i) { return (i <= nGrids / 2) ? i : i - nGrids; };

    std::vector<double> power(nBins, 0.), weights(nBins, 0.), ksum(nBins, 0.);
    double *pw = power.data(), *ws = weights.data(), *ks = ksum.data();

    #pragma omp parallel for reduction(+ : pw[:nBins], ws[:nBins], ks[:nBins])
    for (int r = 0; r < nRows; ++r)
    {
        double n2rows = 0.; // squared mode numbers along the leading axes
        if (D == 2) n2rows = pow(signedIndex(r), 2);
        if (D == 3) n2rows = pow(signedIndex(r / nGrids), 2) + pow(signedIndex(r % nGrids), 2);

        for (int j = 0; j < n_cplx; ++j)
        {
            double n = sqrt(n2rows + double(j) * j); // |k| in units of dk
            int bin = int(n);
            if (bin >= nBins)
                continue;
            // modes with 0 < j < nGrids/2 also stand for their partner at -k, which FFTW does not store
            double w = (j == 0 || j == nGrids / 2) ? 1.0 : 2.0;
            const double *f = phi_k[(size_t)r * n_cplx + j];
            double kc = n * dk;
            pw[bin] += w * prefactor * pow(kc, D) * Pnorm * (f[0] * f[0] + f[1] * f[1]);
            ws[bin] += w;
            ks[bin] += w * kc;
        }
    }

    for (int bin = 0; bin < nBins; ++bin)
        if (weights[bin] > 0)
        {
            double kMean = ksum[bin] / weights[bin]; // comoving
            dat << kMean << "   " << kMean / a << "   " << power[bin] / weights[bin] << std::endl;
        }
}

//------------------------------
// Time evolution
//------------------------------

// Force on pi at the current whole step: dpi = a^(D-2+ALPHAVAL) Laplacian(phi) - a^(D+ALPHAVAL) V'(phi)
void Bubble::computeForces()
{
#pragma omp parallel for
    for (int iGrid = 0; iGrid < N; iGrid++)
    {
        double ddphi = 0.0;

        if (D == 1)
        {
            ddphi = phi[iGridMFunc(iGrid,1)]
                  - 2.0 * phi[iGrid]
                  + phi[iGridMFunc(iGrid,0)];
        }
        else
        {
            for (int iShift = 0; iShift < D; iShift++)
            {
                ddphi += phi[iGridPFunc(iGrid,iShift)]
                       - 2.0 * phi[iGrid]
                       + phi[iGridMFunc(iGrid,iShift)];
            }
        }

        dpi[iGrid] =
              pow(a, D-2+ALPHAVAL) * ddphi / (dx*dx)
            - pow(a, D+ALPHAVAL)   * Vepspr(phi[iGrid], c);
    }
}

// a'' from the Friedmann equations, with the energies at the current whole step (3D expansion only).
// With the pressure P = Ekin - Egrad/3 - Epot: (a'/a)^2 = Etot / (3 Mpl^2) and a''/a = -(Etot + 3P) / (6 Mpl^2).
double Bubble::addot() const
{
    if (!expanding()) return 0.0;
    if (ALPHAVAL == 0) {
        // cosmic time
        return - a * (2.0 * Ekintot - Epottot) / (3.0 * Mpl * Mpl);
    } else {
        // conformal time
        return   a*a*a * (2 * Epottot - Ekintot + Egradtot) / (3.0 * Mpl * Mpl);
    }
}

// Box-averaged energies and field averages at whole step n. Called in oneStep after the kick,
// when pi is already at n+1/2; the kinetic energy is averaged over the half steps n-1/2 and n+1/2.
void Bubble::measureEnergies()
{
    // pi lives at n+1/2 here (called after the kick, before the drift and the update of a),
    // so the kinetic energy needs a_{n+1/2}. ap is still ap_{n-1/2} -> estimate with O(dt^2) error
    double a_kin = a + 0.5 * TIME_STEP * ap;

    Etot = Ekintot = Epottot = Egradtot = 0.0;
    avphi = avphiP = 0.0;

    // reset the spreads
    EtotVar = EgradtotVar = EkintotVar = EpottotVar = 0.;
    avphiPVar = 0.;
    avphiVar = 0.;


    #pragma omp parallel for reduction(+ : Ekintot, Epottot, Egradtot, avphi, avphiP, EkintotVar, EpottotVar, EgradtotVar, EtotVar, avphiVar, avphiPVar)
    for (int iGrid = 0; iGrid < N; iGrid++)
    {
        double gradphix = 0.;
        double gradphiy = 0.;
        double gradphiz = 0.;
        double Egradgrid = 0 ; 
        double Ekingrid = 0 ; 
        double Epotgrid = 0 ;

        // forward differences: consistent with the compact Laplacian in computeForces,
        // so that the measured energy is the one conserved by the evolution
        if (D == 1)
        {
            gradphix = (phi[iGridMFunc(iGrid,1)] - phi[iGrid]) / (a*dx);
        }
        if (D == 2 || D == 3)
        {
            gradphix = (phi[iGridPFunc(iGrid,1)] - phi[iGrid]) / (a*dx);
            gradphiy = (phi[iGridPFunc(iGrid,0)] - phi[iGrid]) / (a*dx);
        }

        if (D == 3)
        {
            gradphiz = (phi[iGridPFunc(iGrid,2)] - phi[iGrid]) / (a*dx);
        }

        double gradphix2 = gradphix * gradphix ;
        double gradphiy2 = gradphiy * gradphiy ;
        double gradphiz2 = gradphiz * gradphiz ;

        Egradgrid = 0.5 * (gradphix2 + gradphiy2 + gradphiz2);
        Epotgrid = Veps(phi[iGrid], c, DeltaV); 
        // physical kinetic energy (dphi/dt)^2 / 2 in cosmic time t: pi = a^(D-ALPHAVAL) phi' with phi' = a^ALPHAVAL dphi/dt
        Ekingrid = 0.5 * pow(a_kin, -2*D) * pi[iGrid]*pi[iGrid];

        Ekintot += Ekingrid * dxLD;
        Epottot += Epotgrid * dxLD;
        Egradtot += Egradgrid * dxLD;

        avphi += phi[iGrid] * dxLD;
        avphiP += pi[iGrid]* pow(a_kin, -D)  * dxLD;   // <dphi/dt>

        // spreads (see EtotVar at the top of the file)
        EkintotVar += Ekingrid * dxLD*Ekingrid * dxLD;
        EpottotVar += Epotgrid * dxLD * Epotgrid * dxLD;
        EgradtotVar += Egradgrid * dxLD * Egradgrid * dxLD;
        EtotVar += (Ekingrid + Egradgrid + Epotgrid) * dxLD* (Ekingrid + Egradgrid + Epotgrid) * dxLD;
        avphiVar += phi[iGrid] * dxLD * phi[iGrid] * dxLD;
        avphiPVar += pi[iGrid]* pow(a_kin, -D) * dxLD * pi[iGrid]* pow(a_kin, -D) * dxLD;

        if (memorysave == 0)
        {
            gradphi[iGrid] = Egradgrid ; 
            potEnergy[iGrid] = Epotgrid ;
            Ekin[iGrid] = Ekingrid ;
        
        }
        

        }
        
    // time-centered kinetic energy at integer time n: average of the raw half-step values n-1/2 and n+1/2
    double EkinHalf = Ekintot;
    Ekintot = 0.5 * (EkinHalf + EkinHalfPrev);
    EkinHalfPrev = EkinHalf;

    Etot = Ekintot + Epottot + Egradtot; // at whole step n
}

// One leapfrog step from whole step n to n+1. phi and a live at whole steps, pi and a' at half steps:
// kick pi to n+1/2 with the force at n, measure the energies at n, update a' (n+1/2) and a (n+1),
// then drift phi to n+1.
void Bubble::oneStep()
{
    computeForces();    // force at step n, from phi_n and a_n

    // kick: pi from n-1/2 to n+1/2
    #pragma omp parallel for
    for (int iGrid = 0; iGrid < N; iGrid++)
        pi[iGrid] += TIME_STEP * dpi[iGrid];

    measureEnergies();  // energies at step n

    double ap_prev = ap;                    // a' at n-1/2
    ap += TIME_STEP * addot();              // a' at n+1/2
    double ap_atn = (ap_prev + ap)/2.0 ;    // a' at n
    double a_old = a;                       // a at n
    a += TIME_STEP * ap;                    // a at n+1
    double a_half = 0.5 * (a + a_old);      // a at n+1/2

    // drift: phi from n to n+1, with phi' = pi / a^(D-ALPHAVAL)
    #pragma omp parallel for
    for (int iGrid = 0; iGrid < N; iGrid++)
        phi[iGrid] += TIME_STEP * pi[iGrid] * pow(a_half, -D + ALPHAVAL);
    t += TIME_STEP;

    H = ap_atn / a_old ;  // Hubble rate at n
    // Hubble constraint at n: should stay small, of order dt^2 (zero without expansion)
    C = expanding() ? H*H - pow(a_old, 2*ALPHAVAL) * Etot / (3.0 * Mpl *Mpl) : 0.0;
}


//------------------------------
// Output
//------------------------------

// Saves field snapshot itPlot as _phi_t=NNNN.h5: datasets phi and pi (and Ekin, gradphi, potEnergy if
// memorysave = 0) with shape (n), (n, n) or (n, n, n) in lattice order [i][j][k], as float. Large lattices are
// thinned to every reduc-th point per axis (2D above 2048^2, 3D above 512^3) to keep files manageable.
// The file attributes give D, L, nGrids, reduc, c, t and a, so a snapshot can be used on its own.
void Bubble::OutField(int itPlot)
{
    std::string filepath = folder + "/_phi_t=" + dataName(itPlot) + ".h5";

    int reduc = 1;
    if (D == 2 && nGrids > 2048) reduc = nGrids / 2048;
    if (D == 3 && nGrids > 512)  reduc = nGrids / 512;
    int n = (nGrids + reduc - 1) / reduc;   // saved points per axis (rounded up, as in the loops below)

    // copy the saved points into float arrays, in the order of the lattice index iGrid
    size_t total = (size_t)pow((double)n, D);
    std::vector<const double *> sources = {phi, pi};
    std::vector<std::string> names = {"phi", "pi"};
    if (memorysave != 1)
    {
        sources.insert(sources.end(), {Ekin, gradphi, potEnergy});
        names.insert(names.end(), {"Ekin", "gradphi", "potEnergy"});
    }
    std::vector<std::vector<float>> data(sources.size(), std::vector<float>(total));

    size_t count = 0;
    int n0 = (D >= 3) ? nGrids : 1, n1 = (D >= 2) ? nGrids : 1;
    for (int i0 = 0; i0 < n0; i0 += (D >= 3 ? reduc : 1))
        for (int i1 = 0; i1 < n1; i1 += (D >= 2 ? reduc : 1))
            for (int i2 = 0; i2 < nGrids; i2 += reduc)
            {
                size_t iGrid = ((size_t)i0 * n1 + i1) * nGrids + i2;
                for (size_t f = 0; f < sources.size(); ++f)
                    data[f][count] = static_cast<float>(sources[f][iGrid]);
                ++count;
            }

    try
    {
        H5File file(filepath.c_str(), H5F_ACC_TRUNC);

        hsize_t dims[3] = {(hsize_t)n, (hsize_t)n, (hsize_t)n};
        DataSpace space(D, dims);
        // compressed in chunks of whole rows (1D: whole field, 2D: groups of rows, 3D: planes)
        hsize_t chunk[3] = {(hsize_t)n, (hsize_t)n, (hsize_t)n};
        if (D == 2) chunk[0] = std::min<hsize_t>(n, std::max<hsize_t>(1, 65536 / n));
        if (D == 3) chunk[0] = 1;
        DSetCreatPropList plist;
        plist.setChunk(D, chunk);
        plist.setShuffle();   // improves compression
        plist.setDeflate(9);  // gzip, maximal compression
        for (size_t f = 0; f < sources.size(); ++f)
            file.createDataSet(names[f], PredType::NATIVE_FLOAT, space, plist)
                .write(data[f].data(), PredType::NATIVE_FLOAT);

        auto attribute = [&](const std::string &name, double value) {
            file.createAttribute(name, PredType::NATIVE_DOUBLE, DataSpace(H5S_SCALAR))
                .write(PredType::NATIVE_DOUBLE, &value);
        };
        attribute("D", D);
        attribute("L", L);
        attribute("nGrids", nGrids);
        attribute("reduc", reduc);
        attribute("c", c);
        attribute("t", t);
        attribute("a", a);
    }
    catch (H5::Exception &error)
    {
        error.printErrorStack();
        std::cerr << "Error: could not write the snapshot " << filepath << std::endl;
    }
}

// Creates the time-series output files with their headers: EnergyDensities.dat (written by outEP)
// and, with wall tracking, gamma.dat (written by outGamma)
void Bubble::InitOutputFiles()
{
    FILE* fp = fopen((folder + "/EnergyDensities.dat").c_str(), "w");
    if (!fp) throw std::runtime_error("Could not open EnergyDensities.dat for writing");
    fprintf(fp, "# c and DeltaV are: %.2f %.4f \n", c, DeltaV);
    fprintf(fp, "# columns: 1 t, 2 a, 3 Etot, 4 Ekin, 5 Egrad, 6 Epot, 7-10 spread of Etot/Ekin/Egrad/Epot,\n");
    fprintf(fp, "#          11 <phi>, 12 spread of <phi>, 13 <phi'>, 14 spread of <phi'>, 15 H, 16 Hubble constraint C = H^2 - rho/(3 Mpl^2)\n");
    fclose(fp);

    if (WallTracking == 1)
    {
        fp = fopen((folder + "/gamma.dat").c_str(), "w");
        if (!fp) throw std::runtime_error("Could not open gamma.dat for writing");
        fprintf(fp, "#t, gamma factor, wall width, rmin, rmax\n");
        fclose(fp);
    }
}

// Appends one line with the box-averaged energies, field averages and Hubble rate to EnergyDensities.dat
void Bubble::outEP()
{
    FILE *fp;
    fp = fopen((folder + "/EnergyDensities.dat").c_str(), "a");
    fprintf(fp, "%.6e %.10e %.12e %.12e %.12e %.12e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.8e %.8e\n",
    t, a, Etot, Ekintot, Egradtot, Epottot,
    sqrt(EtotVar), sqrt(EkintotVar),
    sqrt(EgradtotVar), sqrt(EpottotVar),
    avphi, sqrt(avphiVar),
    avphiP, sqrt(avphiPVar),
    H, C);
    fclose(fp);
}


// 1D: appends the full field as one line to FullField.dat (every Full1DFieldEvery steps)
void Bubble::outFullField()
{
    FILE *fp;
    fp = fopen((folder + "/FullField.dat").c_str(), "a");

    for (int iGrid = 0; iGrid < nGrids; ++iGrid)
    {
        fprintf(fp, "%.5f ", phi[iGrid]);
    }
    fprintf(fp, "\n");

    fflush(fp);
    fclose(fp);
}

// appends t, gamma, wall width, Rmin and Rmax to gamma.dat
void Bubble::outGamma()
{
    FILE *fp;
    fp = fopen((folder + "/gamma.dat").c_str(), "a");
    fprintf(fp, "%.3f %.4f %.3f %.3f %.3f\n ", t, GammaFactor(WallWidth), WallWidth, Rmin, Rmax);
    fclose(fp);
}

// writes all run parameters to parameters.dat (read by the plotting scripts)
void Bubble::OutParameters()
{
    FILE* fp = fopen((folder + "/parameters.dat").c_str(), "w");
    if (!fp) throw std::runtime_error("Could not open parameters.dat for writing");

    fprintf(fp, "# The data in this folder correspond to the following:\n");
    // %.10g: full precision, so that the plotting scripts get the exact values
    fprintf(fp, " D = %d , L = %.10g,  nGrids = %d , TIME_STEP = %.10g  \n", D, L, nGrids, (double)TIME_STEP);
    fprintf(fp, " c = %.10g , nt = %d,  readxns = %d , expansion = %d  \n", c, nt, readxns, expansion);
    fprintf(fp, " nb = %d , FFT = %d,  ALPHAVAL = %d , WallTracking = %d  \n", nBubbles, FFT, ALPHAVAL, WallTracking);
    fprintf(fp, "Mpl = %.10g, memorysave = %d, TEST = %d  \n", (double)Mpl, memorysave, TEST);
    fprintf(fp, " Full1DFieldEvery = %d , Trad = %.10g , fixedseed = %d \n", Full1DFieldEvery, (double)Trad, fixedseed);
    fprintf(fp, " kUV = %.10g , kIR = %.10g , kRStar = %.10g , kMeff = %.10g \n", kUV, kIR, kRStar, kMeff);
    fclose(fp);
}

// Makes a video <type>.mp4 from the snapshot plots <type>*.png in the run folder (needs ffmpeg)
void Bubble::CreateVideoIfFound(const std::string &type)
{
    bool found = false;
    for (const auto &entry : fs::directory_iterator(folder))
    {
        std::string filename = entry.path().filename().string();
        if (filename.find(type) == 0 && entry.path().extension() == ".png")
        {
            found = true;
            break;
        }
    }

    if (found && std::system("command -v ffmpeg > /dev/null 2>&1") != 0)
    {
        std::cout << "ffmpeg not found, skipping video creation.\n";
        return;
    }
    if (found)
    {
        // H.264 needs an even width and height: pad the frames by one pixel where needed
        std::string video = folder + "/" + type + ".mp4";
        std::string command = "ffmpeg -y -framerate 5 -pattern_type glob -i '" +
                              folder + "/" + type + "*.png' -vf 'pad=ceil(iw/2)*2:ceil(ih/2)*2' "
                              "-c:v libx264 -pix_fmt yuv420p '" + video + "' > /dev/null 2>&1";
        if (std::system(command.c_str()) != 0)
        {
            std::cout << "Video creation with ffmpeg failed.\n";
            fs::remove(video); // no empty video file
        }
    }
}
