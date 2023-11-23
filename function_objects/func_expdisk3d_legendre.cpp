/* FILE: func_expdisk3d_legendre.cpp -------------------------------------------- */
/* 
 *   Function object class for a 3D exponential disk (luminosity
 * density = radial exponential with scale length h and vertical sech^(2/n) profile
 * with scale heigh h_z), seen at position angle PA and inclination inc.
 * Integration performed within integration bounds using gauss legendre quadrature.
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
 */

// Copyright 2012--2022 by Peter Erwin.
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
#include <functional>
#include <iostream>
#include <ostream>

#include "func_expdisk3d_legendre.h"
#include "legendre_roots_weights.h"


using namespace std;

/* ---------------- Definitions ---------------------------------------- */
const int   N_PARAMS = 6;
const char  PARAM_LABELS[][20] = {"PA", "inc", "J_0", "h", "n", "z_0"};
const char  PARAM_UNITS[][30] = {"deg (CCW from +y axis)", "deg", "counts/voxel", 
                                "pixels", "", "pixels"};
const char  FUNCTION_NAME[] = "ExponentialDisk3DLegendre function";
const double  DEG2RAD = 0.017453292519943295;
const int  SUBSAMPLE_R = 10;

const double  COSH_LIMIT = 100.0;

const double  INTEGRATION_MULTIPLIER = 20;

const char ExponentialDisk3DLegendre::className[] = "ExponentialDisk3DLegendre";




/* ---------------- Local Functions ------------------------------------ */

double LuminosityDensityLegendre( double s, void *params );

double ExponentialDisk3DLegendre::DistanceAlongMajorAxis(double x_diff, double y_diff){
    double h_los, projected_distance, y_p;
    h_los = sqrt(pow((h * sinInc),2.0) + pow((z_0 * cosInc), 2.0));
    y_p = -x_diff*cosPA +y_diff*sinPA;
    projected_distance = fabs(y_p);
  return projected_distance / h_los;
  };

typedef struct
{ 
  // Parameters of integration.
  double interval; // the size of interval
  double inner_scale, outer_scale;
  int number_of_intervals, roots_number;
  std::function<double(double, void*)> integrand_function;// density function
} IntegrationParameters;

struct Point{
  // 2 dimensional point.
  double x;
  double y;
};

typedef struct Bounds{
  //Limits of integration.
  double lower, upper;
  Bounds(double lower, double upper){
    this->lower = lower;
    this->upper = upper;
  }
} Bounds;

double GaussLegendre(
            std::function<double(double, void*)> luminosity_density,
            Bounds limits,
            double xyParameters[11],
            std::vector<double> roots, 
            std::vector<double> weights
    ){
        /*Integrate luminosity_density.
        *
        * Integration performed within integration bounds using
        * gauss legendre quadrature.
        *
        * Args:
        *    parameters: the list of parameters for integrand function
        *    luminosity_density: integrand function
        *    limits: a Bounds object with lower and upper integration limit
        * Return:
        *    the integral value
        */
        double new_root;
        double half_difference = (limits.upper - limits.lower)/2.0;
        double half_sum = (limits.upper + limits.lower)/2.0;
        double WeightedLuminosityDensity = 0.0;

        for (int i = 0 ; i< roots.size(); i++){
          new_root = half_difference*roots[i] + half_sum;
          WeightedLuminosityDensity +=  weights[i]*LuminosityDensityLegendre(new_root, xyParameters);
        }
        
        return half_difference*WeightedLuminosityDensity;
      };


IntegrationParameters CreateQuadratureIntegration(double distance_to_center,
                                                  double base_interval,
                                                  std::function<double(double, void*)> integrand_function){
  using namespace std::placeholders;
  IntegrationParameters quadrature_integration;
  quadrature_integration.interval = base_interval;
  quadrature_integration.integrand_function = std::bind(integrand_function, std::placeholders::_1, std::placeholders::_2);
  if (distance_to_center < 2.0){
    quadrature_integration.inner_scale = 1.0;
    quadrature_integration.outer_scale = 36.0;
    quadrature_integration.number_of_intervals = 4;
    quadrature_integration.roots_number = 50;
  } else if ((distance_to_center >=2.0) && (distance_to_center < 2.8)) {
    quadrature_integration.inner_scale = 1.0;
    quadrature_integration.number_of_intervals = 2;
    quadrature_integration.outer_scale = 50.0;
    quadrature_integration.roots_number = 100;
  } else if ((distance_to_center >=2.8) && (distance_to_center < 5.2)){
    quadrature_integration.inner_scale = 1.0;
    quadrature_integration.number_of_intervals = 2;
    quadrature_integration.outer_scale = 50.0;
    quadrature_integration.roots_number = 50;
  } else if ((distance_to_center>=5.2) && (distance_to_center < 12.8)){
    quadrature_integration.inner_scale = 1.0;
    quadrature_integration.number_of_intervals = 2;
    quadrature_integration.outer_scale = 50.0;
    quadrature_integration.roots_number = 40;
  } else {
    quadrature_integration.inner_scale = 1.0;
    quadrature_integration.number_of_intervals = 2;
    quadrature_integration.outer_scale = 50.0;
    quadrature_integration.roots_number = 30;
  }            
  return quadrature_integration;
};


double IntegrateLuminosityDensity(double s_plane, 
                                  IntegrationParameters integration_parameters,
                                  double xyParameters[11]
                                  ){

  double inner = integration_parameters.interval;
  double outer = integration_parameters.interval;
  std::vector<double> roots, weights;
  std::vector<Bounds> integrate_limits;
  inner *= integration_parameters.inner_scale;
  outer *= integration_parameters.outer_scale;
  int interval_number = integration_parameters.number_of_intervals;
  if (interval_number == 4){
    integrate_limits = {
                Bounds(s_plane - outer, s_plane - inner),
                Bounds(s_plane - inner, s_plane),
                Bounds(s_plane, s_plane + inner),
                Bounds(s_plane + inner, s_plane + outer),
    };
  } else if (interval_number == 2){
    integrate_limits = {
                Bounds(s_plane - outer, s_plane),
                Bounds(s_plane, s_plane + outer),
    };
  } else{
    integrate_limits = {Bounds(-100, 100)};
  };
  double integral_value = 0.0;
  for (Bounds limit : integrate_limits) {
    switch(integration_parameters.roots_number) {
      case 30:
        roots = ROOTS_30;
        weights = WEIGHTS_30;
        break;
      case 40:
        roots = ROOTS_40;
        weights = WEIGHTS_40;
        break;
      case 50:
        roots = ROOTS_50;
        weights = WEIGHTS_50;
        break;
      case 100:
        roots = ROOTS_100;
        weights = WEIGHTS_100;
        break;
      default:
        roots = ROOTS_50;
        weights = WEIGHTS_50;
        break;
    }
    integral_value += GaussLegendre(
      integration_parameters.integrand_function,
      limit,
      xyParameters,
      roots, weights);

}   
  return integral_value;
};


/* ---------------- CONSTRUCTOR ---------------------------------------- */

ExponentialDisk3DLegendre::ExponentialDisk3DLegendre( )
{
  
  nParams = N_PARAMS;
  functionName = FUNCTION_NAME;
  shortFunctionName = className;

  // Set up vectors of parameter labels and units
  for (int i = 0; i < nParams; i++) {
    parameterLabels.push_back(PARAM_LABELS[i]);
    parameterUnits.push_back(PARAM_UNITS[i]);
  }
  parameterUnitsExist = true;
  
  doSubsampling = false;
}


/* ---------------- PUBLIC METHOD: Setup ------------------------------- */

void ExponentialDisk3DLegendre::Setup( double params[], int offsetIndex, double xc, double yc )
{
  x0 = xc;
  y0 = yc;
  PA = params[0 + offsetIndex];
  inclination = params[1 + offsetIndex];
  J_0 = params[2 + offsetIndex ];
  h = params[3 + offsetIndex ];
  n = params[4 + offsetIndex ];
  z_0 = params[5 + offsetIndex ];

  // pre-compute useful things for this round of invoking the function
  // convert PA to +x-axis reference
  PA_rad = (PA + 90.0) * DEG2RAD;
  cosPA = cos(PA_rad);
  sinPA = sin(PA_rad);
  inc_rad = inclination * DEG2RAD;
  cosInc = cos(inc_rad);
  sinInc = sin(inc_rad);
  
  alpha = 2.0/n;
  scaledZ0 = alpha*z_0;
  two_to_alpha = pow(2.0, alpha);

  base_interval = max(h * sinInc, z_0 * cosInc);
}


/* ---------------- PUBLIC METHOD: GetValue ---------------------------- */

double ExponentialDisk3DLegendre::GetValue( double x, double y )
{
  double  x_diff = x - x0;
  double  y_diff = y - y0;
  double  xp, yp, x_d0, y_d0, z_d0, totalIntensity;
  double  xyParameters[11];
  double distance;
  double s_plane;
  
  // Calculate x,y in component (projected sky) reference frame
  xp = x_diff*cosPA + y_diff*sinPA;
  yp = -x_diff*sinPA + y_diff*cosPA;

  // Calculate (x,y,z)_start in component's native xyz reference frame, corresponding to
  // intersection of line-of-sight ray with projected sky frame
  x_d0 = xp;
  y_d0 = yp * cosInc;
  z_d0 = yp * sinInc;

  // Set up parameter vector for the integration (everything that stays unchanged
  // for this particular xp,yp location)
  xyParameters[0] = x_d0;
  xyParameters[1] = y_d0;
  xyParameters[2] = z_d0;
  xyParameters[3] = cosInc;
  xyParameters[4] = sinInc;
  xyParameters[5] = J_0;
  xyParameters[6] = h;
  xyParameters[7] = z_0;
  xyParameters[8] = scaledZ0;
  xyParameters[9] = two_to_alpha;
  xyParameters[10] = alpha;

  

  distance = DistanceAlongMajorAxis(x_diff, y_diff);

  IntegrationParameters integration_parameters = CreateQuadratureIntegration(distance,
                                                        base_interval,
                                                        LuminosityDensityLegendre);

  if (inclination < 90.0){
    s_plane = yp * tan(inc_rad);
  }
  else{
    s_plane = 0.0;
  };

  totalIntensity = IntegrateLuminosityDensity(s_plane, integration_parameters, xyParameters);

  return totalIntensity;
}





/* ----------------------------- OTHER FUNCTIONS -------------------------------- */








/* Compute luminosity density for a location (x_d,y_d,z_d) which is at line-of-sight 
 * distance s from start point (x_d0, y_d0, z_d0), where midplane of component (e.g.,
 * disk of galaxy) is oriented at angle (90 - inclination) to the line of sight vector. 
 */ 
double LuminosityDensityLegendre( double s, void *params )
{
  double  y_d, z_d, z, R, lumDensity;
  double  verticalScaling, sech;
  double  *paramsVect = (double *)params;
  double x_d0 = paramsVect[0];
  double y_d0 = paramsVect[1];
  double z_d0 = paramsVect[2];
  double cosInc = paramsVect[3];
  double sinInc = paramsVect[4];
  double J_0 = paramsVect[5];
  double h = paramsVect[6];
  double z_0 = paramsVect[7];
  double scaledZ0 = paramsVect[8];
  double two_to_alpha = paramsVect[9];
  double alpha = paramsVect[10];
  
  // Given s and the pre-defined parameters, determine our 3D location (x_d,y_d,z_d)
  // [by construction, x_d = x_d0]
  y_d = y_d0 + s*sinInc;
  z_d = z_d0 - s*cosInc;
  
  // Convert 3D Cartesian coordinate to R,z coordinate
  R = sqrt(x_d0*x_d0 + y_d*y_d);
  z = fabs(z_d);

  // if combination of n*z/z_0 is large enough, switch to simple exponential
  if ((z/scaledZ0) > COSH_LIMIT)
    verticalScaling = two_to_alpha * exp(-z/z_0);
  else {
    sech = 1.0 / cosh(z/scaledZ0);
    verticalScaling = pow(sech, alpha);
  }

  lumDensity = J_0 * exp(-R/h) * verticalScaling;
  return lumDensity;
}



/* END OF FILE: func_expdisk3d.cpp ------------------------------------- */
