#include "GetVariablesWithMap.h"
#include "VariableWithMap.h"
#include "CVUniverse.h"
#include "CCQEBinning.h"

//==============================================================================
// Define variables
// Hists are defined inside of the CCProtonPi0::Variable class definition (Variable.h).
//==============================================================================
//std::vector<CCProtonPi0::VariableWithMap*> GetVariables(variables, tags) {
//  // variables allows you to code a lot of functions but only chose the ones with names you pick
//  // tags allows you to associate variables with groups of cuts or other selections
//  //=========================================
//  // Constructor wants: name, x-axis label,
//  //                    binning,
//  //                    CVUniverse reco and
//  //                    truth functions
//  //=========================================
//
//
//
//  example of an exclusive variable not used here...
//  // Pion kinetic energy
//  // An exclusive variable is one that corresponds to a track, cluster, michel,
//  // etc.
//  // It is one that accepts a track-identifer integer index.
//  //  std::vector<double> tpi_bins = {35., 68., 100., 133., 166., 200., 350.};
//  //  typedef ExclusiveVariable1Arg<CVUniverse, CCProtonPi0::Variable> EVar1;
//  //  EVar1* tpi = new EVar1("tpi", "T#pi (MeV)", tpi_bins, &CVUniverse::GetTpi,
//  //                         &CVUniverse::GetTpi);
//


// version with tags you can use

std::vector< CCProtonPi0::VariableWithMap*> GetVariablesWithMap(const std::vector<std::string> vars, const std::vector<std::string> tags) {
  //=========================================
  // Constructor wants: name, x-axis label,
  //                    binning,
  //                    CVUniverse reco and
  //                    truth functions
  //=========================================
  std::vector< CCProtonPi0::VariableWithMap*> allvariables;  // set this internally to span the set of variables
  std::vector< CCProtonPi0::VariableWithMap*> variables; // this is the set that actually gets returned
  
  // ok, so configure a bunch of histograms here.
  std::vector<double> enu_bins = {0., 1., 2. , 3. , 4. , 5. , 6. , 7. , 8. , 9. , 10. , 15. };
  std::vector<double> pt_bins = {0,0.075,0.15,0.25,0.325,0.4,0.475,0.55,0.7,0.85,1,1.25,1.5,2.5};
  std::vector<double> pi0Mass_bins = {0,0.02,0.04,0.06,0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.16,0.18,0.20}; 
 
    std::vector<double> logq2_bins = GetLogQ2BinsGeV();
    std::vector<double> q2_bins = GetQ2BinsGeV();

    CCProtonPi0::VariableWithMap* enuQE = new CCProtonPi0::VariableWithMap("enuQE", "E_{#nu}(QE) (GeV)", enu_bins,&CVUniverse::GetEnuCCQEGeV, &CVUniverse::GetTrueEnuGeV);
    CCProtonPi0::VariableWithMap* pT = new CCProtonPi0::VariableWithMap("pT", "muon p_{T} (GeV/c)", pt_bins,&CVUniverse::GetPperpMuGeV, &CVUniverse::GetTruePperpMuGeV);
    CCProtonPi0::VariableWithMap* total = new CCProtonPi0::VariableWithMap("total", "muon p_{T} (GeV/c)",1,1.5,15.,&CVUniverse::GetPparMuGeV, &CVUniverse::GetTruePparMuGeV);
    CCProtonPi0::VariableWithMap* pZ = new CCProtonPi0::VariableWithMap("pZ", "muon p_{||} (GeV/c)", enu_bins,&CVUniverse::GetPparMuGeV, &CVUniverse::GetTruePparMuGeV);
    CCProtonPi0::VariableWithMap* Log10Q2QE = new CCProtonPi0::VariableWithMap("Log10Q2QE", "Log_{10} Q^{2} (QE) (GeV^{2})", logq2_bins,&CVUniverse::GetLog10Q2QEGeV, &CVUniverse::GetTrueLog10Q2QEGeV);
    CCProtonPi0::VariableWithMap* Q2QE = new CCProtonPi0::VariableWithMap("Q2QE", "Q^{2} (QE) (GeV^{2})", q2_bins,&CVUniverse::GetQ2QEGeV, &CVUniverse::GetTrueQ2QEGeV);
    CCProtonPi0::VariableWithMap* thetaX = new CCProtonPi0::VariableWithMap("thetaX", "#theta_{x} radians",40,-.1,.1,&CVUniverse::GetThetaXmu,&CVUniverse::GetTrueThetaXmu);
    CCProtonPi0::VariableWithMap* thetaY = new CCProtonPi0::VariableWithMap("thetaY", "#theta_{y} radians",40,-.1,.1,&CVUniverse::GetThetaYmu,&CVUniverse::GetTrueThetaYmu);
    CCProtonPi0::VariableWithMap* theta = new CCProtonPi0::VariableWithMap("theta", "#theta radians",20,0.,.1,&CVUniverse::GetThetamu,&CVUniverse::GetTrueThetamu);
    CCProtonPi0::VariableWithMap* recoil = new CCProtonPi0::VariableWithMap("recoil", "recoil (GeV)",25,0.,0.5,&CVUniverse::GetRecoilEnergyGeV,&CVUniverse::GetTrueQ0GeV);
    CCProtonPi0::VariableWithMap* pi0Mass = new CCProtonPi0::VariableWithMap("pi0Mass", "pi0Mass (GeV)", pi0Mass_bins, &CVUniverse::GetRecoPi0Mass, &CVUniverse::GetTruePi0Mass);
  // make an internal list of all the points which you will then prune based on variables you actually want for this run
  
  allvariables = {total,enuQE,pT,pZ,Log10Q2QE,Q2QE,thetaX,thetaY,theta,recoil,pi0Mass};//,enuQE,Q2QE,recoil,thetaX,thetaY};
 
  for (auto var:vars){
    bool found = false;
    for (auto variable:allvariables){
      if(var == variable->GetName()){
        std::cout << " study variable " << var << std::endl;
        // this is the point where you add the tags.  Saves space this way.
        variable->AddTags(tags);
        variables.push_back(variable);
        found = true;
        break;
      }
    }
    if (!found){
      std::cout << " Warning - have requested an unimplemented variable in GetVariableWithMap " << var << std::endl;
    }
  }
  return variables;
}
