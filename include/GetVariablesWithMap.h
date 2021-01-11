#ifndef GETVARIABLESWithMap_H
#define GETVARIABLESWithMap_H
#include <vector>
#include "CVUniverse.h"
#include "PlotUtils/PhysicsVariables.h"

namespace CCProtonPi0
{
  class VariableWithMap;
}

std::vector<CCProtonPi0::VariableWithMap*> GetVariablesWithMap(const std::vector<std::string> variables, const std::vector<std::string> tags );
#endif
