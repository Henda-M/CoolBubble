/*==============================================================================
 *  CoolBubble -- lattice simulation of vacuum bubble nucleation, growth and
 *  collision during cosmological first-order phase transitions.
 *
 *  File: gridfuncs.cpp
 *  Role: Implementation of periodic-lattice helpers -- neighbor indexing and
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

#include <cstdio>
#include <vector>

#include "gridfuncs.h"
#include "def.h"

// Index conventions: in 2D the flat index is iGrid = nGrids*i + j, in 3D
// iGrid = nGrids^2*i + nGrids*j + k. Neighbours wrap around periodically.
// In 1D only iGridMFunc is used: iShift = 1 gives the next point, iShift = 0 the previous one.

// Index of the neighbor of iGrid one step "below/left" along axis iShift
// (iShift: 0=y, 1=x, 2=z), wrapping around periodically at the box edge.
int iGridMFunc(int iGrid, int iShift)
{
    if(D==1)
    {
        if(iShift == 1) return (iGrid + nGrids + 1)%nGrids; //moves forward
        if(iShift == 0) return (iGrid + nGrids - 1)%nGrids; //moves back
    } 
    if(D == 2) return nGrids*((iGrid/nGrids + nGrids - (iShift == 0 ? 1 : 0))%nGrids) + (iGrid + nGrids - (iShift == 1 ? 1 : 0))%nGrids;
    if(D == 3) return nGrids*(nGrids*((iGrid/(nGrids*nGrids) + nGrids - (iShift == 0 ? 1 : 0))%nGrids) + ((iGrid/nGrids + nGrids - (iShift == 1 ? 1 : 0))%nGrids)) + (iGrid + nGrids - (iShift == 2 ? 1 : 0))%nGrids;
    return 0;
}

// Index of the neighbor of iGrid one step "above/right" along axis iShift
// (iShift: 0=y, 1=x, 2=z), wrapping around periodically at the box edge.
// Not implemented for D=1 (unused there, see comment above).
int iGridPFunc(int iGrid, int iShift)
{
    if(D == 2) return nGrids*((iGrid/nGrids + (iShift == 0 ? 1 : 0))%nGrids) + (iGrid + (iShift == 1 ? 1 : 0))%nGrids;
    if(D == 3) return nGrids*(nGrids*((iGrid/(nGrids*nGrids) + (iShift == 0 ? 1 : 0))%nGrids) + (iGrid/nGrids + (iShift == 1 ? 1 : 0))%nGrids) + (iGrid + (iShift == 2 ? 1 : 0))%nGrids;
    return 0;
}


// Converts a flat array index iGrid into its {i,j,k} lattice indices
// (unused trailing components are 0 for D<3).
std::vector<int> iGridsFunc(int iGrid) //gives a vector {i,j,k}: lattice indices corresponding to lattice iGrid
{
    std::vector<int> result = {0,0,0}; 
    if(D == 1) result= {iGrid%nGrids, 0, 0}; //gives the indices of the grid/2D array elements of size nGrid^2
    if(D == 2) result= {(iGrid/nGrids)%nGrids, iGrid%nGrids, 0}; //gives the indices of the grid/2D array elements of size nGrid^2
    if(D == 3) result= {(iGrid/(nGrids*nGrids))%nGrids, (iGrid/nGrids)%nGrids, iGrid%nGrids}; //gives the indices of the grid/3D array elements
    return result; 
}

// Inverse of iGridsFunc: converts {i,j,k} lattice indices back into the
// corresponding flat array index. Not implemented for D=1.
int ijkTOiGrids(std::vector<int> x) //gives the corresponding iGrids to a vector {i,j,k}
{
    if(D == 2) return x[0]*nGrids+x[1]; //gives the indices of the grid/2D array elements of size nGrid^2
    if(D == 3) return  x[0]*nGrids*nGrids+x[1]*nGrids+x[2]; //gives the indices of the grid/3D array elements
    return 0;
}


// Converts a flat array index iGrid into cell-centered physical coordinates
// (x,y,z) in a box of size lL, i.e. x_i = lL*(i+0.5)/nGrids per dimension.
std::vector<double> xFunc(int iGrid, double lL)
{
    std::vector<int> iGrids = iGridsFunc(iGrid); //gives a 3D vector with the indices of the array elements and zero in the last component: (n,m,0)
    std::vector<double> result = {0,0,0}; 
    if (D == 1) result= {lL*((double)iGrids[0] + 0.5)/(double)nGrids, 0 ,0 };
    if(D == 2) result= {lL*((double)iGrids[0] + 0.5)/(double)nGrids, lL*((double)iGrids[1] + 0.5)/(double)nGrids, 0.};
    if(D == 3) result= {lL*((double)iGrids[0] + 0.5)/(double)nGrids, lL*((double)iGrids[1] + 0.5)/(double)nGrids, lL*((double)iGrids[2] + 0.5)/(double)nGrids};
    return result;
}

// For a randomly nucleated bubble centred at xn with radius Rn in a periodic box of
// size `shift`, writes to fp the periodic images of xn needed so that a bubble
// crossing the box boundary is nucleated on the opposite side(s) as well.
// Loops over every combination of shifts (-1, 0, +1) along the D axes: an image
// is written if the bubble crosses the boundary on every shifted axis. In 2D this
// gives up to 3 images (edges and corner), in 3D up to 7 (faces, edges and corner).
void check_PBC(std::vector<double> xn, double shift, double Rn, std::FILE *fp)
{
    // allowed shift per axis: -1 if the bubble sticks out at the top, +1 at the bottom, 0 otherwise
    int allowed[3] = {0, 0, 0};
    for (int d = 0; d < D; d++)
    {
        if (xn[d] + Rn > shift) allowed[d] = -1;
        else if (xn[d] - Rn < 0) allowed[d] = +1;
    }

    // each axis is either left alone (bit = 0) or shifted by allowed[d] (bit = 1)
    for (int mask = 1; mask < (1 << D); mask++)
    {
        std::vector<double> image = xn;
        bool valid = true;
        for (int d = 0; d < D; d++)
        {
            if (mask & (1 << d))
            {
                if (allowed[d] == 0) { valid = false; break; }
                image[d] += allowed[d] * shift;
            }
        }
        if (valid)
            fprintf(fp, "%.6e %.6e %.6e\n", image[0], image[1], image[2]);
    }
}
