/* FILE: model_object.cpp ---------------------------------------------- */
/* VERSION 0.3
 *
 * This is intended to be an abstract base class for the various
 * "model" objects (e.g., image data + fitting functions).
 *   
 *     [v0.3]:  4 Dec 2009: Added handling of mask images.
 *     [v0.2]: 27 Nov 2009: Modified to include AddDataVectors function, which
 * will be used by derived class ModelObject1D
 *     [v0.1]: 13--15 Nov 2009: Created; initial development.
 *
 */


/* ------------------------ Include Files (Header Files )--------------- */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <iostream>

#include "definitions.h"
#include "model_object.h"
#include "mp_enorm.h"


/* ---------------- Definitions ---------------------------------------- */
static string  UNDEFINED = "<undefined>";


/* ---------------- CONSTRUCTOR ---------------------------------------- */

ModelObject::ModelObject( )
{
  dataValsSet = weightValsSet = false;
  parameterBoundsSet = false;
  parameterBounds = NULL;
  modelVector = NULL;
  modelVectorAllocated = false;
  weightVectorAllocated = false;
  modelImageComputed = false;
  maskExists = false;
  nFunctions = 0;
  nFunctionSets = 0;
  nFunctionParams = 0;
  nParamsTot = 0;
}


/* ---------------- PUBLIC METHOD: AddFunction ------------------------- */

void ModelObject::AddFunction( FunctionObject *newFunctionObj_ptr )
{
  int  nNewParams;
  
  functionObjects.push_back(newFunctionObj_ptr);
  nFunctions += 1;
  nNewParams = newFunctionObj_ptr->GetNParams();
  paramSizes.push_back(nNewParams);
  nFunctionParams += nNewParams;
  newFunctionObj_ptr->GetParameterNames(parameterLabels);
}



/* ---------------- PUBLIC METHOD: DefineFunctionSets ------------------ */

void ModelObject::DefineFunctionSets( vector<int>& functionStartIndices )
{
  int  nn, i;
  
  nFunctionSets = functionStartIndices.size();
    // define array of [false, false, false, ...]
  setStartFlag = (bool *)calloc(nFunctions, sizeof(bool));
  for (i = 0; i < nFunctionSets; i++) {
    nn = functionStartIndices[i];
    // function number n is start of new function set; 
    // change setStartFlag[n] to true
    setStartFlag[nn] = true;
  }
  
  // total number of parameters = number of parameters for individual functions
  // plus x0/y0 pair for each function set
  nParamsTot = nFunctionParams + 2*nFunctionSets;
}



/* ---------------- PUBLIC METHOD: AddDataVectors --------------------- */

void ModelObject::AddDataVectors( int nDataValues, double *xValVector, double *yValVector,
																	bool magnitudeData )
{
  // Just a placeholder for now (needs to be modified & overridden in derived
  // class ModelObject1D
  nDataVals = nDataValues;
}


/* ---------------- PUBLIC METHOD: AddImageDataVector ------------------ */

void ModelObject::AddImageDataVector( int nDataValues, int nImageColumns,
                                      int nImageRows, double *pixelVector )
{
  nDataVals = nValidDataVals = nDataValues;
  nColumns = nImageColumns;
  nRows = nImageRows;
  dataVector = pixelVector;
  dataValsSet = true;
  
  SetupModelImage(nDataVals, nColumns, nRows);
}


/* ---------------- PUBLIC METHOD: SetupModelImage -------------------- */
// Called by AddImageDataVector(); can also be used by itself in make-image
// mode
void ModelObject::SetupModelImage( int nDataValues, int nImageColumns,
                                      int nImageRows )
{
  nDataVals = nDataValues;
  nColumns = nImageColumns;
  nRows = nImageRows;
  
  // Allocate modelimage vector
  modelVector = (double *) calloc((size_t)nDataVals, sizeof(double));
  modelVectorAllocated = true;
}


/* ---------------- PUBLIC METHOD: AddErrorVector ---------------------- */

void ModelObject::AddErrorVector( int nDataValues, int nImageColumns,
                                      int nImageRows, double *pixelVector,
                                      int inputType )
{
  assert( (nDataValues == nDataVals) && (nImageColumns == nColumns) && 
          (nImageRows == nRows) );
  weightVector = pixelVector;
  
  // Convert noise values into weights, if needed.  Our standard approach is
  // to compute weights as 1/sigma; this assumes that whatever function calls
  // ComputeDeviates() will then square the individual (weighted) deviate values
  // in order to get the proper chi^2 result.
  // Currently, we assume three possibilities for weight-map pixel values:
  //    sigma (std.dev.); variance (sigma^2); and plain weights
  //    Note that correct interpretation of chi^2 values depends on weights
  //    being based on sigmas or variances!
  switch (inputType) {
    case WEIGHTS_ARE_SIGMAS:
      for (int z = 0; z < nDataVals; z++) {
        weightVector[z] = 1.0 / weightVector[z];
      }
      break;
    case WEIGHTS_ARE_VARIANCES:
      for (int z = 0; z < nDataVals; z++) {
        weightVector[z] = 1.0 / sqrt(weightVector[z]);
      }
      break;
    default:
      // do nothing, since input values *are* weights (e.g., the input image is
      // a weight map with pixel values = 1/sigma already)
      ;
  }
      
  weightValsSet = true;
}


/* ---------------- PUBLIC METHOD: AddErrorVector1D -------------------- */
// This is a stub function; it is meant to be properly defined in the derived
// class ModelObject1d
void ModelObject::AddErrorVector1D( int nDataValues, double *inputVector,
                                      int inputType )
{
  ;
}



/* ---------------- PUBLIC METHOD: GenerateErrorVector ----------------- */
// Generate an error vector based on Poisson statistics.
void ModelObject::GenerateErrorVector( double gain, double readNoise, double skyValue )
{
  double  sky_plus_readNoise, noise_squared;
  
  // Allocate storage for weight image:
  weightVector = (double *) calloc((size_t)nDataVals, sizeof(double));
  weightVectorAllocated = true;
  
  sky_plus_readNoise = skyValue/gain + readNoise*readNoise/(gain*gain);
  for (int z = 0; z < nDataVals; z++) {
    noise_squared = dataVector[z]/gain + sky_plus_readNoise;
    weightVector[z] = 1.0 / sqrt(noise_squared);
  }

  weightValsSet = true;
}



/* ---------------- PUBLIC METHOD: AddMaskVector ----------------------- */

void ModelObject::AddMaskVector( int nDataValues, int nImageColumns,
                                      int nImageRows, double *pixelVector,
                                      int inputType )
{
  assert( (nDataValues == nDataVals) && (nImageColumns == nColumns) && 
          (nImageRows == nRows) );

  maskVector = pixelVector;
  nValidDataVals = 0;   // Since there's a mask, not all pixels from the original
                        // image will be valid
    
  // We need to convert the mask values so that good pixels = 1 and bad
  // pixels = 0.
  switch (inputType) {
    case MASK_ZERO_IS_GOOD:
      // This is our "standard" input mask: good pixels are zero, bad pixels
      // are positive integers
      printf("ModelObject::AddMaskVector -- zero is good ...\n");
      for (int z = 0; z < nDataVals; z++) {
        if (maskVector[z] > 0.0) {
          maskVector[z] = 0.0;
        } else {
          maskVector[z] = 1.0;
          nValidDataVals++;
        }
      }
      break;
    case MASK_ZERO_IS_BAD:
      // Alternate form for input masks: good pixels are 1, bad pixels are 0
      printf("ModelObject::AddMaskVector -- zero is bad ...\n");
      for (int z = 0; z < nDataVals; z++) {
        if (maskVector[z] < 1.0)
          maskVector[z] = 0.0;
        else {
          maskVector[z] = 1.0;
          nValidDataVals++;
        }
      }
      break;
    default:
      printf("ModelObject::AddMaskVector -- WARNING: unknown inputType detected!\n\n");
      exit(-1);
  }
      
  maskExists = true;
}


/* ---------------- PUBLIC METHOD: ApplyMask --------------------------- */

void ModelObject::ApplyMask( )
{
  if ( (weightValsSet) && (maskExists) ) {
    for (int z = 0; z < nDataVals; z++) {
      weightVector[z] = maskVector[z] * weightVector[z];
    }
    printf("ModelObject: mask vector applied to weight vector. ");
    printf("(%d valid pixels remain)\n", nValidDataVals);
  }
  else {
    printf(" ** ALERT: ModelObject::ApplyMask() called, but we are missing either\n");
    printf("    error image or mask image, or both!  ApplyMask() ignored ...\n");
  }
}


/* ---------------- PUBLIC METHOD: AddParameterLimits ------------------ */

void ModelObject::AddParameterLimits( double *paramLimits )
{
  parameterBounds = paramLimits;
  parameterBoundsSet = true;
}


/* ---------------- PUBLIC METHOD: CreateModelImage -------------------- */

void ModelObject::CreateModelImage( double params[] )
{
  double  x0, y0, x, y, newVal;
  int  offset = 0;
  
  // Separate out the individual-component parameters and tell the
  // associated function objects to do setup work.
  // The first component's parameters start at params[0]; the second's
  // start at params[paramSizes[0]], the third at 
  // params[paramSizes[0] + paramSizes[1]], and so forth...
  for (int n = 0; n < nFunctions; n++) {
    if (setStartFlag[n] == true) {
      // start of new function set: extract x0,y0 and then skip over them
      x0 = params[offset];
      y0 = params[offset + 1];
      offset += 2;
      printf("  Function %d = start of new set; x0 = %g, y0 = %g\n",
              n, x0, y0);
    }
    functionObjects[n]->Setup(params, offset, x0, y0);
    offset += paramSizes[n];
  }
  
  // populate modelVector with the model image
  for (int i = 0; i < nRows; i++) {   // step by row number = y
    y = (double)(i + 1);
    for (int j = 0; j < nColumns; j++) {   // step by column number = x
      x = (double)(j + 1);
      newVal = 0.0;
      for (int n = 0; n < nFunctions; n++)
        newVal += functionObjects[n]->GetValue(x, y);
      modelVector[i*nColumns + j] = newVal;
      //printf("   x = %g, y = %g, value = %g\n", x, y, modelVector[i*nColumns + j]);
    }
  }
  modelImageComputed = true;
}


/* ---------------- PUBLIC METHOD: ComputeDeviates --------------------- */
/* This function computes the vector of weighted deviates (differences between
 * model and data pixel values).  Note that a proper chi^2 sum requires *squaring*
 * each deviate value before summing them; we assume this is done elsewhere, by 
 * whatever function calls ComputeDeviates().
 */
void ModelObject::ComputeDeviates( double yResults[], double params[] )
{

  CreateModelImage(params);
  
  for (int z = 0; z < nDataVals; z++) {
    //printf("weightVector[%d] = %f, ", z, weightVector[z]);
    yResults[z] = weightVector[z] * (dataVector[z] - modelVector[z]);
    //printf("yResults[%d] = %f  ", z, yResults[z]);
  }
}


/* ---------------- PUBLIC METHOD: ChiSquared -------------------------- */
/* Function for calculating chi^2 value for a model.  Current version is meant
 * to be used once (e.g., after fitting is done); for general, repetitive use
 * (e.g., with Diff'l Evoln.), we should allocate deviates[] separately, rather
 * than allocating and then freeing it as part of this function.
 */
double ModelObject::ChiSquared( double params[] )
{
  double  *deviates;
  double  chi;
  
  CreateModelImage(params);
  
  deviates = (double *) malloc(nDataVals * sizeof(double));

  for (int z = 0; z < nDataVals; z++)
    deviates[z] = weightVector[z] * (dataVector[z] - modelVector[z]);
  chi = mp_enorm(nDataVals, deviates);
  free(deviates);
  
  return chi*chi;
}


/* ---------------- PUBLIC METHOD: GetDescription ---------------------- */

void ModelObject::GetDescription( )
{
  printf("Model Object: %d data values (pixels)\n", nDataVals);
}


/* ---------------- PUBLIC METHOD: PrintImage ------------------------- */

void ModelObject::PrintImage( )
{

  if (! dataValsSet) {
    printf("* ModelObject: No image data supplied!\n\n");
    return;
  }
  
  printf("The whole image, row by row:\n");
  // The following fetches pixels row-by-row, starting with the bottom
  // row (i.e., what we would normally like to think of as the first row)
  for (int i = 0; i < nRows; i++) {   // step by row number = y
    for (int j = 0; j < nColumns; j++)   // step by column number = x
      printf(" %f", dataVector[i*nColumns + j]);
    printf("\n");
  }
  printf("\n");
}


/* ---------------- PUBLIC METHOD: PrintModelImage -------------------- */

void ModelObject::PrintModelImage( )
{

  if (! modelImageComputed) {
    printf("* ModelObject: Model image has not yet been computed!\n\n");
    return;
  }
  
  printf("The model image, row by row:\n");
  // The following fetches pixels row-by-row, starting with the bottom
  // row (i.e., what we would normally like to think of as the first row)
  for (int i = 0; i < nRows; i++) {   // step by row number = y
    for (int j = 0; j < nColumns; j++)   // step by column number = x
      printf(" %f", modelVector[i*nColumns + j]);
    printf("\n");
  }
  printf("\n");
}


/* ---------------- PUBLIC METHOD: GetParameterName -------------------- */

string& ModelObject::GetParameterName( int i )
{
  if (i < nParamsTot)
    return parameterLabels[i];
  else
	return UNDEFINED;
}


/* ---------------- PUBLIC METHOD: GetNParams -------------------------- */

int ModelObject::GetNParams( )
{
  return nParamsTot;
}


/* ---------------- PUBLIC METHOD: GetNValidPixels --------------------- */

int ModelObject::GetNValidPixels( )
{
  return nValidDataVals;
}


/* ---------------- PUBLIC METHOD: GetModelImageVector ----------------- */

double * ModelObject::GetModelImageVector( )
{
  if (! modelImageComputed) {
    printf("* ModelObject: Model image has not yet been computed!\n\n");
    return NULL;
  }
  
  return modelVector;
}


/* ---------------- DESTRUCTOR ----------------------------------------- */

ModelObject::~ModelObject()
{
  if (modelVectorAllocated)
    free(modelVector);
  if (weightVectorAllocated)
    free(weightVector);
  
  if (nFunctions > 0)
    for (int i = 0; i < nFunctions; i++)
      delete functionObjects[i];
  free(setStartFlag);
}



/* END OF FILE: model_object.cpp --------------------------------------- */