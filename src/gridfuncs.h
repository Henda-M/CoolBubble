/*==============================================================================
 *  CoolBubble -- lattice simulation of vacuum bubble nucleation, growth and
 *  collision during cosmological first-order phase transitions.
 *
 *  File: gridfuncs.h
 *  Role: Declarations for periodic-lattice helpers -- neighbor indexing and
 *        conversions between flat grid index, (i,j,k) indices and physical
 *        coordinates.
 *
 *  See README.md for build/usage instructions and def.h for run parameters.
 *
 *  Copyright (C) 2026 Henda Mansour, Yann Gouttenoire
 *  SPDX-License-Identifier: GPL-3.0-or-later
 *  This program is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License, version 3 or (at your
 *  option) any later version. See the LICENSE file for details.
 *==============================================================================*/

#ifndef GRID_FUNCS_H
#define GRID_FUNCS_H


#include <cstdio>
#include <vector>

// Run parameters, defined in main.cpp (defaults from def.h, overridable on the command line)
extern int D;
extern int nGrids;
extern double L;
extern int nt; 
extern int readxns;
extern double c; 


// Function declarations
int iGridMFunc(int iGrid, int iShift);
int iGridPFunc(int iGrid, int iShift);
std::vector<int> iGridsFunc(int iGrid); //gives a vector {i,j,k}: lattice indices corresponding to lattice iGrid
int ijkTOiGrids(std::vector<int> x); //gives the corresponding iGrids to a vector {i,j,k}
std::vector<double> xFunc(int iGrid, double lL);
void check_PBC(std::vector<double> xn, double shift, double Rn, std::FILE *fp); //makes sure that randomely nucleated bubbles respect the periodic boundary conditions



#endif // GRID_FUNCS_H
