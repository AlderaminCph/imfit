/* FILE: add_functions.cpp ----------------------------------------------- */
/*
 * Function which takes a vector of strings listing function names and generates
 * the corresponding FunctionObjects, passing them to the input ModelObject
 *
 */

#include <string>
#include <vector>

#include "model_object.h"

#include "function_object.h"
#include "func_gaussian.h"
#include "func_sersic.h"
#include "func_exp.h"

using namespace std;


int AddFunctions( ModelObject *theModel, vector<string> &functionNameList,
                  vector<int> &functionSetIndices )
{
  int  nFunctions = functionNameList.size();
  string  currentName;
  FunctionObject  *thisFunctionObj;
  
  for (int i = 0; i < nFunctions; i++) {
    currentName = functionNameList[i];
    printf("Function: %s\n", currentName.c_str());
    if (currentName == "Exponential") {
      thisFunctionObj = new Exponential();
      theModel->AddFunction(thisFunctionObj);
      continue;
    }
    if (currentName == "Sersic") {
      thisFunctionObj = new Sersic();
      theModel->AddFunction(thisFunctionObj);
      continue;
    }
    if (currentName == "Gaussian") {
      thisFunctionObj = new Gaussian();
      theModel->AddFunction(thisFunctionObj);
      continue;
    }
    // If we reach here, then something went wrong
    printf("*** AddFunctions: unidentified function name (\"%s\")\n", currentName.c_str());
    return - 1;
  }
  
  // Tell model object about arrangement of functions into common-center sets
  theModel->DefineFunctionSets(functionSetIndices);
  return 0;
}


/* END OF FILE: add_functions.cpp ---------------------------------------- */