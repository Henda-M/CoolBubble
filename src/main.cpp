/*==============================================================================
 *  CoolBubble -- lattice simulation of vacuum bubble nucleation, growth and
 *  collision during cosmological first-order phase transitions.
 *
 *  File: main.cpp
 *  Role: Program entry point -- command-line options, output folder setup,
 *        and the main time-evolution loop.
 *
 *  Supports 1+1D, 2+1D and 3+1D runs, in a static or expanding background.
 *  See README.md for build/usage instructions and def.h for run parameters.
 *
 *  Copyright (C) 2026 Henda Mansour, Yann Gouttenoire
 *  SPDX-License-Identifier: GPL-3.0-or-later
 *  This program is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License, version 3 or (at your
 *  option) any later version. See the LICENSE file for details.
 *==============================================================================*/

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "def.h"       // default values of the parameters
#include "gridfuncs.h" // lattice indexing helpers
#include "bubble.hpp"  // the Bubble class

namespace fs = std::filesystem;
using namespace std;

// Folder holding the plotting scripts. Set by the makefile to <repo>/python.
#ifndef PYTHON_DIR
#define PYTHON_DIR "python"
#endif

//-----------------------------
// Run parameters shared with the Bubble class and the grid helpers (see gridfuncs.h).
// They start from the defaults in def.h and can be overridden on the command line.
//-----------------------------
int D = D_default;           // number of spatial dimensions
int nGrids = nGrids_default; // number of grid points per dimension
double L = L_default;        // box size
int nt = nt_default;         // number of time steps
int readxns = readxns_default; // initial condition / nucleation mode
double c = c_default;        // potential parameter

// Options that only matter in this file
static string outDir = Path_default;  // where the run folder is created
static string pythonDir = PYTHON_DIR; // where the plotting scripts live
static bool makePlots = true;         // call the plotting scripts during the run
static int nSnap = nSnapshots;        // number of field snapshots after the initial one

//------------------------------
// Helpers
//------------------------------

// number as a string with a fixed number of decimals, used in the folder name
static string to_precistr(double num, int precisionVal = 0)
{
    return to_string(num).substr(0, to_string(num).find(".") + precisionVal + 1);
}

static void printUsage(const char *prog)
{
    cout <<
        "Usage: " << prog << " [options]\n"
        "\n"
        "Every option is optional; missing ones take their default from src/def.h.\n"
        "\n"
        "  --D <1|2|3>       number of spatial dimensions            (default " << D_default << ")\n"
        "  --L <value>       box size, in units of 1/v_phi           (default " << L_default << ")\n"
        "  --nGrids <n>      grid points per dimension               (default " << nGrids_default << ")\n"
        "  --nt <n>          number of time steps                    (default " << nt_default << ")\n"
        "  --c <value>       potential parameter                     (default " << c_default << ")\n"
        "  --mode <m>        initial condition (2D/3D):              (default " << readxns_default << ")\n"
        "                      1  nBubbles bubbles at random positions\n"
        "                      0  bubbles at the positions listed in <run folder>/__xns.dat\n"
        "                     -1  two bubbles on the box diagonal\n"
        "                     -2  one bubble in the centre\n"
        "                     -3  thermal fluctuations instead of bubbles\n"
        "                    (in 1D the run always starts from two colliding bubbles)\n"
        "  --snapshots <n>   number of field snapshots, besides the initial one (default " << nSnapshots << ")\n"
        "  --out <dir>       folder in which the run folder is created (default " << Path_default << ")\n"
        "  --no-plots        do not call the Python plotting scripts\n"
        "  --python-dir <d>  folder with the plotting scripts (default " << PYTHON_DIR << ")\n"
        "  -h, --help        show this help\n"
        "\n"
        "Example: " << prog << " --D 2 --L 64 --nGrids 1024 --nt 1000 --mode 1 --out runs/\n";
}

// Reads the command-line options into the globals above. Exits on --help or on an unknown option.
static void parseArguments(int argc, const char *argv[])
{
    for (int i = 1; i < argc; i++)
    {
        string opt = argv[i];
        if (opt == "-h" || opt == "--help") { printUsage(argv[0]); exit(0); }
        if (opt == "--no-plots") { makePlots = false; continue; }

        if (i + 1 >= argc)
        {
            cerr << "Error: option " << opt << " needs a value (see --help)\n";
            exit(1);
        }
        string val = argv[++i];
        try
        {
            if      (opt == "--D")          D = stoi(val);
            else if (opt == "--L")          L = stod(val);
            else if (opt == "--nGrids")     nGrids = stoi(val);
            else if (opt == "--nt")         nt = stoi(val);
            else if (opt == "--c")          c = stod(val);
            else if (opt == "--mode")       readxns = stoi(val);
            else if (opt == "--snapshots")  nSnap = stoi(val);
            else if (opt == "--out")        outDir = val;
            else if (opt == "--python-dir") pythonDir = val;
            else
            {
                cerr << "Error: unknown option " << opt << " (see --help)\n";
                exit(1);
            }
        }
        catch (const std::exception &)
        {
            cerr << "Error: cannot read the value '" << val << "' of option " << opt << "\n";
            exit(1);
        }
    }

    if (D < 1 || D > 3)                    { cerr << "Error: --D must be 1, 2 or 3\n"; exit(1); }
    if (nGrids < 2 || L <= 0 || nt < 1)    { cerr << "Error: --nGrids, --L and --nt must be positive\n"; exit(1); }
    if (readxns < -3 || readxns > 1)       { cerr << "Error: --mode must be between -3 and 1\n"; exit(1); }
    if (nSnap < 1)                         { cerr << "Error: --snapshots must be at least 1\n"; exit(1); }
}

// Calls the Python script that plots one field snapshot and its power spectrum
static void PlotFieldAndPS(const string &folder, const string &filename, const Bubble &s, double gamma = 1.0)
{
    if (!makePlots || plotsnapshots == 0) return;

    string command = "python3 \"" + pythonDir + "/plot_FieldAndPowerSpectrum.py\"";
    command += " --folder \"" + folder + "\"";
    command += " --filename " + filename;

    if (WallTracking == 1) {
        command += " --gamma " + to_string(gamma);
        command += " --rmin " + to_string(s.Rmin);
        command += " --rmax " + to_string(s.Rmax);
    }

    if (system(command.c_str()) != 0)
        cerr << "Warning: plotting failed for " << filename << "\n";
}

// Calls the Python script that makes the summary plots at the end of the run
static void SummaryPlots(const string &folder)
{
    if (!makePlots) return;

    string command = "python3 \"" + pythonDir + "/plot_SummaryPlots.py\" --folder \"" + folder + "\"";
    if (system(command.c_str()) != 0)
        cerr << "Warning: summary plots failed\n";
}

// The leapfrog scheme is stable only if dt < dx/sqrt(D)
static bool CheckStability(double L, int D, int nGrids)
{
    double dx = double(L) / nGrids;
    double dtmax = dx / sqrt(D);
    if (TIME_STEP < dtmax)
    {
        cout << "Stability condition dt < dx/sqrt(D) is satisfied." << endl;
        return true;
    }
    cout << "Stability condition dt < dx/sqrt(D) is not satisfied (dt = " << TIME_STEP
         << ", dx/sqrt(D) = " << dtmax << "). Stopping." << endl;
    return false;
}

// File name of field snapshot n
static string SnapshotFile(int n)
{
    return "_phi_t=" + dataName(n) + ".h5";
}

// Saves field snapshot n with its power spectrum, and plots them
static void SaveSnapshot(Bubble &s, int n, const string &folder)
{
    cout << "snapshot " << n << " at t = " << s.t << endl;
    s.OutField(n);
    s.FourierTransform(n);
    PlotFieldAndPS(folder, SnapshotFile(n), s, WallTracking == 1 ? s.GammaFactor(s.WallWidth) : 0.0);
}

// Time steps of snapshots 1..n (see SnapshotSpacing in def.h). The last one is always the final step.
static vector<int> SnapshotSchedule(int n)
{
    double T = nt * TIME_STEP; // duration of the run
    vector<double> times;
    if (SnapshotSpacing == 1 && T > L && n >= 2)
    {
        int nEarly = (n + 1) / 2, nLate = n - nEarly;
        for (int k = 1; k <= nEarly; k++) times.push_back(L * k / nEarly);            // up to t = L
        for (int k = 1; k <= nLate; k++)  times.push_back(L + (T - L) * k / nLate);   // from t = L to the end
    }
    else
    {
        for (int k = 1; k <= n; k++) times.push_back(T * k / n);
    }

    vector<int> steps;
    for (double time : times)
    {
        int step = min(nt, max(1, (int)lround(time / TIME_STEP)));
        if (steps.empty() || step > steps.back()) steps.push_back(step); // skip repeats in very short runs
    }
    return steps;
}

/*-------------------------------------------------------------*/
/*---               Here begins the main code               ---*/
/*-------------------------------------------------------------*/

int main(int argc, const char *argv[])
{
    auto start = std::chrono::high_resolution_clock::now();

    parseArguments(argc, argv);

    if (makePlots && !fs::exists(fs::path(pythonDir) / "plot_FieldAndPowerSpectrum.py"))
    {
        cerr << "Warning: plotting scripts not found in " << pythonDir
             << ", running without plots (use --python-dir to point to them).\n";
        makePlots = false;
    }

    /*     Generating a folder where files and plots will be saved      */
    string Sname = "Sim"+to_string(D)+"_D_"+"nt_"+to_string(nt)+
    "_dt_"+to_precistr(TIME_STEP,5)+"_nGrids_"+to_string(nGrids)+
    "_L_"+to_precistr(L,1)+"_Mode_"+to_string(readxns)+"_nB_"
    +to_string(nBubbles)+"_c_"+to_precistr(c,2)+"_FFT_"+to_string(FFT)
    +(D == 1 ? "_Full1DFieldEvery_"+to_string(Full1DFieldEvery) : "")+"_Mpl_"+to_string(int(Mpl));

    string folder = (fs::path(outDir) / Sname).string();

    /* Check that choice of time step and grid spacing is alright */
    if (!CheckStability(L, D, nGrids))
        return 1;

    /* Create folder */
    fs::create_directories(folder);
    cout << "Output folder: " << folder << endl;

    /*         initialize a Bubble           */
    Bubble s(folder);
    s.OutParameters();  // save the parameters in parameters.dat
    s.outEP();          // initial energy densities
    SaveSnapshot(s, 0, folder);
    if (D == 1) s.outFullField();

    vector<int> snapshotSteps = SnapshotSchedule(nSnap);
    size_t nextSnapshot = 0;
    bool trackWall = WallTracking == 1 && s.CanTrackWall();
    if (WallTracking == 1 && !trackWall)
        cout << "Note: wall tracking only works in 1D or with --mode -1 or -2; skipping it." << endl;
    int progressEvery = max(1, nt / 20);

    /*----------    Time Evolution Loop     -----------*/
    for (int it = 1; it <= nt; it++)
    {
        s.oneStep();

        if (trackWall)                          s.TrackWall(it);
        if (D == 1 && it % Full1DFieldEvery == 0) s.outFullField();
        if (it % OutEnergyEvery == 0)           s.outEP();
        if (nextSnapshot < snapshotSteps.size() && it == snapshotSteps[nextSnapshot])
        {
            if (trackWall) s.UpdateWall(it);   // current wall position for the plot markers
            SaveSnapshot(s, ++nextSnapshot, folder);
        }
        if (it % progressEvery == 0)
            cout << "step " << it << " / " << nt << " (" << 100 * it / nt << "%)" << endl;
    }

    // 1D: a Fourier transform in time and space can be made afterwards from FullField.dat (phi(x, t) every
    // Full1DFieldEvery steps). Since the evolution is not periodic in time, it needs a window function.

    SummaryPlots(folder);
    // videos of the snapshot plots (field) and of the 2D maps from the summary plots (kinetic and gradient energy)
    if (makePlots && plotsnapshots == 1 && makevideos == 1) s.CreateVideoIfFound("_phi_t");
    if (makePlots) { s.CreateVideoIfFound("Ekin__phi_t"); s.CreateVideoIfFound("Egrad__phi_t"); }

    //Compute how long the run took
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    cout << "Elapsed time: " << elapsed.count() << " seconds\n";

    FILE* fp = fopen((folder + "/parameters.dat").c_str(), "a");
    if (fp)
    {
        fprintf(fp, "Elapsed time %f seconds\n", elapsed.count());
        fclose (fp);
    }

    return 0;
}
