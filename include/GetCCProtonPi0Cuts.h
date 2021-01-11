#include "PlotUtils/Cut.h"
#include "CVUniverse.h"
#include "CCProtonPi0Cuts.h"
#include "CCProtonPi0Signal.h"
#ifndef GETCCPROTONPI0CUTS_H
#define GETCCPROTONPI0CUTS_H

// GetCCProtonPi0Cuts.h was created by Kang Yang kang.yang@physics.ox.ac.uk

// ---------------------- reco namespace start ------------------------
namespace reco {
//============================================================================
// This function instantiates each of the above cuts and adds them to a
// container, over which we'll loop during our event selection to apply the
// cuts.
//
// The return type for this function is a `cuts_t<UNIVERSE, EVENT>`, which is
// shorthand for std::vector<std::unique_ptr<PlotUtils::Cut<UNIVERSE, EVENT>>>;
//============================================================================
template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
PlotUtils::cuts_t<UNIVERSE, EVENT> GetCCProtonPi0Cuts()
{
  PlotUtils::cuts_t<UNIVERSE, EVENT> ccprotonpi0_cuts;
  std::cout << " Set up the cuts for CCProtonPi0" << std::endl;
  // event cuts
  ccprotonpi0_cuts.emplace_back(new NoDeadtime<UNIVERSE, EVENT>(1, "Deadtime"));
  ccprotonpi0_cuts.emplace_back(new Apothem<UNIVERSE>(850.));
  ccprotonpi0_cuts.emplace_back(new ZRange<UNIVERSE>("Tracker", 5980, 8422));
  // muon cuts
  ccprotonpi0_cuts.emplace_back(new HasMINOSMatch<UNIVERSE, EVENT>());
  ccprotonpi0_cuts.emplace_back(new IsNeutrino<UNIVERSE, EVENT>(true));
  ccprotonpi0_cuts.emplace_back(new MaxMuonAngle<UNIVERSE, EVENT>(20.));
  ccprotonpi0_cuts.emplace_back(new MuonMomentumMin<UNIVERSE, EVENT>(1.5, "Minimum Pmu"));
  ccprotonpi0_cuts.emplace_back(new MuonMomentumMax<UNIVERSE, EVENT>(20., "Maximum Pmu"));
  //pi0 mass
  //ccprotonpi0_cuts.emplace_back(new Pi0MassMin<UNIVERSE, EVENT>(0.02, "Minimum pi0Mass"));
  // recoil cuts
  //ccprotonpi0_cuts.emplace_back(new MaxMultiplicity<UNIVERSE, EVENT>(1));
  //ccprotonpi0_cuts.emplace_back(new GoodRecoil<UNIVERSE, EVENT>());
  return ccprotonpi0_cuts;
}

} //----- End of reco namespace



// ---------------------- truth namespace start ------------------------
namespace truth{
//-----------signal
template <class UNIVERSE>
PlotUtils::constraints_t<UNIVERSE> GetCCProtonPi0PhaseSpace()
{
  PlotUtils::constraints_t<UNIVERSE> signalDef;
  
  signalDef.emplace_back(new MuonAngle<UNIVERSE>(20.));
  signalDef.emplace_back(new PmuMin<UNIVERSE>(1.5));
  signalDef.emplace_back(new PmuMax<UNIVERSE>(20.));
  signalDef.emplace_back(new Apothem<UNIVERSE>(850.));
  signalDef.emplace_back(new ZRange<UNIVERSE>("Tracker", 5980, 8422));
  
  return signalDef;
}

// CCProtonPi0 ------------------------
template <class UNIVERSE>
PlotUtils::constraints_t<UNIVERSE> GetCCProtonPi0Signal(double kecut=NSFDefaults::TrueProtonKECutCentral)
{
  PlotUtils::constraints_t<UNIVERSE> signalDef;
  signalDef.emplace_back(new IsNeutrino<UNIVERSE>(14));
  signalDef.emplace_back(new IsCC<UNIVERSE>());
  //signalDef.emplace_back(new IsCCProtonPi0<UNIVERSE>(kecut , true));
  signalDef.emplace_back(new IsFidVol<UNIVERSE>()); 
  signalDef.emplace_back(new IsSignal<UNIVERSE>(true)); 
  return signalDef;
}

template <class UNIVERSE>
PlotUtils::constraints_t<UNIVERSE> GetCCProtonPi0Bkg(double kecut=NSFDefaults::TrueProtonKECutCentral)
{
  PlotUtils::constraints_t<UNIVERSE> signalDef;
  signalDef.emplace_back(new IsNeutrino<UNIVERSE>(14));
  signalDef.emplace_back(new IsCC<UNIVERSE>());
  //signalDef.emplace_back(new IsCCProtonPi0<UNIVERSE>(kecut , false));
  signalDef.emplace_back(new IsFidVol<UNIVERSE>());
  signalDef.emplace_back(new IsSignal<UNIVERSE>(false));
  return signalDef;
}


} //----- End of truth namespace
#endif //GETCCPROTONPI0CUTS_H
