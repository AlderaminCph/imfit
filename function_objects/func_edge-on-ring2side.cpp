/* FILE: func_edge-on-ring2side.cpp ------------------------------------ */
/* 
 *   Highly experimental class (derived from FunctionObject; function_object.h)
 * which produces a crude edge-on ring model; this version allows for different
 * Gaussian sigma for r < r_ring vs r > r_ring.
 *
 *   
 *   BASIC IDEA:
 *      Setup() is called as the first part of invoking the function;
 *      it pre-computes various things that don't depend on x and y.
 *      GetValue() then completes the calculation, using the actual value
 *      of x and y, and returns the result.
 *      So for an image, we expect the user to call Setup() once at
 *      the start, then loop through the pixels of the image, calling
 *      GetValue() to compute the function results for each pixel coordinate
 *      (x,y).
 *
 *   NOTE: Currently, we assume input PA is in *degrees* [and then we
 * convert it to radians] relative to +x axis.
 *
 *   MODIFICATION HISTORY:
 *     [v0.1]  28 Sept 2010: Created as modification of func_edge-on-ring.cpp.
 */

// Copyright 2010--2016 by Peter Erwin.
// 
// This file is part of Imfit.
// 
// Imfit is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your
// option) any later version.
// 
// Imfit is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
// 
// You should have received a copy of the GNU General Public License along
// with Imfit.  If not, see <http://www.gnu.org/licenses/>.


/* ------------------------ Include Files (Header Files )--------------- */
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <string>

#include "func_edge-on-ring2side.h"

using namespace std;


/* ---------------- Definitions ---------------------------------------- */
const int   N_PARAMS = 6;
const char  PARAM_LABELS[][20] = {"PA", "I_0", "r", "sigma_r_in", "sigma_r_out", "sigma_z"};
const char  FUNCTION_NAME[] = "Edge-on Ring function";
const double  DEG2RAD = 0.017453292519943295;
const int  SUBSAMPLE_R = 10;

const char EdgeOnRing2Side::className[] = "EdgeOnRing2Side";


/* ---------------- CONSTRUCTOR ---------------------------------------- */

EdgeOnRing2Side::EdgeOnRing2Side( )
{
  string  paramName;
  
  nParams = N_PARAMS;
  functionName = FUNCTION_NAME;
  shortFunctionName = className;

  // Set up the vector of parameter labels
  for (int i = 0; i < nParams; i++) {
    paramName = PARAM_LABELS[i];
    parameterLabels.push_back(paramName);
  }
  
  doSubsampling = true;
}


/* ---------------- PUBLIC METHOD: Setup ------------------------------- */

void EdgeOnRing2Side::Setup( double params[], int offsetIndex, double xc, double yc )
{
  x0 = xc;
  y0 = yc;
  PA = params[0 + offsetIndex];
  A = params[1 + offsetIndex ];
  r1 = -params[2 + offsetIndex ];   // left-side radius of ring (negative value)
  r2 = -r1;                        // right-side radius of ring (positive value)
  sigma_r_inner = params[3 + offsetIndex ];
  sigma_r_outer = params[4 + offsetIndex ];
  sigma_z = params[5 + offsetIndex ];

  // pre-compute useful things for this round of invoking the function
  // convert PA to +x-axis reference
  PA_rad = (PA + 90.0) * DEG2RAD;
  cosPA = cos(PA_rad);
  sinPA = sin(PA_rad);

  twosigma_rin_squared = 2.0 * sigma_r_inner*sigma_r_inner;
  twosigma_rout_squared = 2.0 * sigma_r_outer*sigma_r_outer;
  twosigma_z_squared = 2.0 * sigma_z*sigma_z;
}


/* ---------------- PRIVATE METHOD: CalculateIntensity ----------------- */
// This function calculates the intensity for an edge-on ring 2D function, at
// cylindrical radius r (along the major-axis of the component) and height z 
// (perpendicular to the major-axis of the component).
// NOTE: This function requires that z be *non-negative*!
double EdgeOnRing2Side::CalculateIntensity( double r, double z )
{
  double  I1, I2;
  double  r1_diff = r - r1;   // distance from left-side ring center
  double  r2_diff = r - r2;   // distance from right-side ring center
  double  z_diff = z;
  double  z_diff_ratio = z_diff*z_diff/twosigma_z_squared;
  double  twosigma_r1_squared, twosigma_r2_squared;
  
  // decide which sigma to use for each Gaussian component
  if (r < r1) {  // outside ring to the left
    twosigma_r1_squared = twosigma_rout_squared;
    twosigma_r2_squared = twosigma_rin_squared;
  }
  else {  // r > r1, so we're to the right of left-hand Gaussian
    if (r < r2) {  // interior of the ring
      twosigma_r1_squared = twosigma_rin_squared;
      twosigma_r2_squared = twosigma_rin_squared;
    }
    else {  // outside ring to the right
      twosigma_r1_squared = twosigma_rin_squared;
      twosigma_r2_squared = twosigma_rout_squared;
    }
  }
  
  I1 = A * exp(-(r1_diff*r1_diff)/twosigma_r1_squared - z_diff_ratio);
  I2 = A * exp(-(r2_diff*r2_diff)/twosigma_r2_squared - z_diff_ratio);
  return I1 + I2;
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */
// NOTE: x and y are orthognal image coordinates; R and z are orthogonal
// coordinates *in the component reference frame* (corresponding to xp and
// yp in other, non-edge-on function objects).
// Note that both CalculateIntensity() and CalculateSubsamples() assume that
// R and z are *non-negative*!
double EdgeOnRing2Side::GetValue( double x, double y )
{
  double  x_diff = x - x0;
  double  y_diff = y - y0;
  double  R, z, totalIntensity;
  int  nSubsamples;
  
  // Calculate R,z (= x,y in component reference frame)
  R = x_diff*cosPA + y_diff*sinPA;    // "R" is x in the component reference frame
  z = fabs(-x_diff*sinPA + y_diff*cosPA);   // "z" is y in the component reference frame

  nSubsamples = CalculateSubsamples(R, z);
  if (nSubsamples > 1) {
    // Do subsampling
    // start in center of leftmost/bottommost sub-picel
    double deltaSubpix = 1.0 / nSubsamples;
    double x_sub_start = x - 0.5 + 0.5*deltaSubpix;
    double y_sub_start = y - 0.5 + 0.5*deltaSubpix;
    double theSum = 0.0;
    for (int ii = 0; ii < nSubsamples; ii++) {
      double x_ii = x_sub_start + ii*deltaSubpix;
      for (int jj = 0; jj < nSubsamples; jj++) {
        double y_ii = y_sub_start + jj*deltaSubpix;
        x_diff = x_ii - x0;
        y_diff = y_ii - y0;
        R = x_diff*cosPA + y_diff*sinPA;
        z = fabs(-x_diff*sinPA + y_diff*cosPA);
        theSum += CalculateIntensity(R, z);
      }
    }
    totalIntensity = theSum / (nSubsamples*nSubsamples);
  }
  else
    totalIntensity = CalculateIntensity(R, z);

  return totalIntensity;
}


/* ---------------- PROTECTED METHOD: CalculateSubsamples ------------------------- */
// Function which determines the number of pixel subdivisions for sub-pixel integration,
// given that the current pixel is a distance of r away from the center of the
// r=0 line and/or z away from the disk plane.
// This function returns the number of x and y subdivisions; the total number of subpixels
// will then be the return value *squared*.
int EdgeOnRing2Side::CalculateSubsamples( double r, double z )
{
  int  nSamples = 1;
  return nSamples;
}



/* END OF FILE: func_edge-on-ring2side.cpp ----------------------------- */
