#ifndef GETVARIABLES_H
#define GETVARIABLES_H
#include <vector>
#include "CVUniverse.h"
#include "PlotUtils/PhysicsVariables.h"

namespace CCProtonPi0
{
  class Variable;
}

std::vector<CCProtonPi0::Variable*> GetVariables(std::string tag = "" );
#endif
