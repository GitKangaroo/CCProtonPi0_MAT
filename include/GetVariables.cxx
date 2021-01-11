#include "GetVariables.h"
#include "Variable.h"
#include "CVUniverse.h"
#include "CCQEBinning.h"

//==============================================================================
// Define variables
// Hists are defined inside of the CCProtonPi0::Variable class definition (Variable.h).
//==============================================================================
//std::vector<CCProtonPi0::Variable*> GetVariables() {
//  //=========================================
//  // Constructor wants: name, x-axis label,
//  //                    binning,
//  //                    CVUniverse reco and
//  //                    truth functions
//  //=========================================
//
//
//  // Neutrino energy
//  std::vector<double> enu_bins = {0., 1., 2. , 3. , 4. , 5. , 6. , 7. , 8. , 9. , 10. , 15. };
//  std::vector<double> logq2_bins = GetLogQ2BinsGeV();
//  CCProtonPi0::Variable* enuQE = new CCProtonPi0::Variable("enuQE", "E_{#nu}(QE) (GeV)", enu_bins,&CVUniverse::GetEnuCCQEGeV, &CVUniverse::GetTrueEnuGeV);
//  CCProtonPi0::Variable* Q2QE = new CCProtonPi0::Variable("Log10Q2QE", "Log_{10} Q^{2} (QE) (GeV^{2})", logq2_bins,&CVUniverse::GetLog10Q2QEGeV, &CVUniverse::GetTrueLog10Q2QEGeV);
//  CCProtonPi0::Variable* thetaX = new CCProtonPi0::Variable("thetaX", "#theta_{x} radians",40,-.1,.1,&CVUniverse::GetThetaXmu,&CVUniverse::GetTrueThetaXmu);
//  CCProtonPi0::Variable* thetaY = new CCProtonPi0::Variable("thetaY", "#theta_{y} radians",40,-.1,.1,&CVUniverse::GetThetaYmu,&CVUniverse::GetTrueThetaYmu);
//  CCProtonPi0::Variable* theta = new CCProtonPi0::Variable("theta", "#theta radians",20,0.,.1,&CVUniverse::GetThetamu,&CVUniverse::GetTrueThetamu);
//  CCProtonPi0::Variable* recoil = new CCProtonPi0::Variable("recoil", "recoil (GeV)",25,0.,0.5,&CVUniverse::GetRecoilEnergyGeV,&CVUniverse::GetTrueQ0GeV);
//
//
//  // Pion kinetic energy
//  // An exclusive variable is one that corresponds to a track, cluster, michel,
//  // etc.
//  // It is one that accepts a track-identifer integer index.
//  //  std::vector<double> tpi_bins = {35., 68., 100., 133., 166., 200., 350.};
//  //  typedef ExclusiveVariable1Arg<CVUniverse, CCProtonPi0::Variable> EVar1;
//  //  EVar1* tpi = new EVar1("tpi", "T#pi (MeV)", tpi_bins, &CVUniverse::GetTpi,
//  //                         &CVUniverse::GetTpi);
//
//  std::vector<CCProtonPi0::Variable*> variables = {enuQE,thetaX,thetaY,theta,Q2QE,recoil};
//  return variables;
//}

// version with tags you can use

std::vector< CCProtonPi0::Variable*> GetVariables( std::string tag) {
  //=========================================
  // Constructor wants: name, x-axis label,
  //                    binning,
  //                    CVUniverse reco and
  //                    truth functions
  //=========================================
  std::vector< CCProtonPi0::Variable*> variables;
  std::vector<double> enu_bins = {0., 1., 2. , 3. , 4. , 5. , 6. , 7. , 8. , 9. , 10. , 15. };
  std::vector<double> pt_bins = {0,0.075,0.15,0.25,0.325,0.4,0.475,0.55,0.7,0.85,1,1.25,1.5,2.5};
     
  std::string tagb;
    if(tag != "")  tagb = tag + "_";
    std::vector<double> logq2_bins = GetLogQ2BinsGeV();
    CCProtonPi0::Variable* enuQE = new CCProtonPi0::Variable(tagb+"enuQE", "E_{#nu}(QE) (GeV)", enu_bins,&CVUniverse::GetEnuCCQEGeV, &CVUniverse::GetTrueEnuGeV);
    CCProtonPi0::Variable* pT = new CCProtonPi0::Variable(tagb+"pT", "muon pT (GeV/c)", pt_bins,&CVUniverse::GetPperpMuGeV, &CVUniverse::GetTruePperpMuGeV);
    CCProtonPi0::Variable* Q2QE = new CCProtonPi0::Variable(tagb+"Log10Q2QE", "Log_{10} Q^{2} (QE) (GeV^{2})", logq2_bins,&CVUniverse::GetLog10Q2QEGeV, &CVUniverse::GetTrueLog10Q2QEGeV);
    CCProtonPi0::Variable* thetaX = new CCProtonPi0::Variable(tagb+"thetaX", "#theta_{x} radians",40,-.1,.1,&CVUniverse::GetThetaXmu,&CVUniverse::GetTrueThetaXmu);
    CCProtonPi0::Variable* thetaY = new CCProtonPi0::Variable(tagb+"thetaY", "#theta_{y} radians",40,-.1,.1,&CVUniverse::GetThetaYmu,&CVUniverse::GetTrueThetaYmu);
    CCProtonPi0::Variable* theta = new CCProtonPi0::Variable(tagb+"theta", "#theta radians",20,0.,.1,&CVUniverse::GetThetamu,&CVUniverse::GetTrueThetamu);
    CCProtonPi0::Variable* recoil = new CCProtonPi0::Variable(tagb+"recoil", "recoil (GeV)",25,0.,0.5,&CVUniverse::GetRecoilEnergyGeV,&CVUniverse::GetTrueQ0GeV);
    variables = {pT,enuQE,Q2QE,recoil,thetaX,thetaY};
  return variables;
}
