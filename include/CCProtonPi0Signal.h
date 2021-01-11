//File: CCProtonPi0Signal.h
//Brief: Signal definition constraints for a CCProtonPi0 analysis.
//Author: Andrew Olivier aolivier@ur.rochester.edu

// CCProtonPi0Signal.h was created by Kang Yang kang.yang@physics.ox.ac.uk

// Public-use truth::Cuts Classes here defined:
// -------------------------------------------
// IsCCProtonPi0
// IsNeutrino
// IsCC
// MuonAngle
// PmuMin
// PmuMax
// Apothem
// ZRange

#ifndef CCPROTONPI0SIGNAL_H
#define CCPROTONPI0SIGNAL_H

namespace truth {
  // This function is needed when calling IsCCProtonPi0
  bool passTrueCCProtonPi0(bool neutrinoMode, std::vector<int> mc_FSPartPDG, std::vector<double> mc_FSPartE, int mc_nFSPart, double proton_KECut){
    int genie_n_muons = 0;
    int genie_n_mesons        = 0;
    int genie_n_heavy_baryons = 0;
    int genie_n_photons       = 0;
    int genie_n_protons       = 0; 
    int genie_n_pi0s          = 0;
  
    for(int i = 0; i < mc_nFSPart; ++i) {
      int pdg =  mc_FSPartPDG[i];
      double energy = mc_FSPartE[i];
      double KEp = energy-MinervaUnits::M_p;
      // implement 120 MeV KE cut, if needed
    
    
      //The photon energy cut is hard-coded at 10 MeV at present. We're happy to make it general, if the need arises !
      if( abs(pdg) == 13 ) genie_n_muons++;
      else if( pdg == 22 && energy >10 ) genie_n_photons++;
      else if( abs(pdg) == 211 || abs(pdg) == 321 || abs(pdg) == 323 || pdg == 130 || pdg == 310 || pdg == 311 || pdg == 313 ) genie_n_mesons++;
      else if( pdg == 3112 || pdg == 3122 || pdg == 3212 || pdg == 3222 || pdg == 4112 || pdg == 4122 || pdg == 4212 || pdg == 4222 || pdg == 411 || pdg == 421) genie_n_heavy_baryons++;
      else if( pdg == 2212 && KEp > proton_KECut) genie_n_protons++;
      else if( pdg == 111) genie_n_pi0s++;
    }
    //Definition of CCProtonPi0: 1 muon (from neutrino) and at least one pi0 in final state
    //Any number of final state nucleons (protons or neutrons) allowed - except no protons above 120 GeV (check for CCProtonPi0)
    //Photons from nuclear de-excitation are kept. These tend to be < 10 MeV. 
    //GENIE simulates nuclear de-excitations only for Oxygen atoms at present.
    
    if(neutrinoMode){ // Note: no limit on final state nucleons
      if( genie_n_muons         == 1 &&
          genie_n_mesons        == 0 &&
          genie_n_heavy_baryons == 0 &&
          genie_n_photons       == 0 &&
          genie_n_pi0s          >= 1  ) return true;
    }
    else{
      if( genie_n_muons         == 1 &&
          genie_n_mesons        == 0 &&
          genie_n_heavy_baryons == 0 &&
          genie_n_photons       == 0 &&
          genie_n_pi0s          >= 1  ) return true; 
    }
    return false;
  }
 
// IsCCProtonPi0 -- check if this is a CCProtonPi0 signal (see the definition above) 
template <class UNIVERSE>
class IsCCProtonPi0: public PlotUtils::SignalConstraint<UNIVERSE>
{
public:
  IsCCProtonPi0(double protonKEcut,bool isSignal): PlotUtils::SignalConstraint<UNIVERSE>("IsCCProtonPi0"),fKE(protonKEcut),fsig(isSignal)
  {
  }
  
private:
  bool checkConstraint(const UNIVERSE& univ) const override
  {
    std::vector<int>mc_FSPartPDG = univ.GetVecInt("mc_FSPartPDG");
    std::vector<double>mc_FSPartE = univ.GetVecDouble("mc_FSPartE");
    bool neutrinoMode = univ.GetTruthNuPDG() > 0;
    int mc_nFSPart = univ.GetInt("mc_nFSPart");
    int mc_incoming = univ.GetInt("mc_incoming");
    int mc_current = univ.GetInt("mc_current");
 
    bool passes = ( passTrueCCProtonPi0(neutrinoMode,  mc_FSPartPDG, mc_FSPartE, mc_nFSPart, fKE));
  
    return passes == fsig;
  }
  const bool fsig;
  const int fKE;
};


// IsNeutrino
template <class UNIVERSE>
class IsNeutrino: public PlotUtils::SignalConstraint<UNIVERSE>
{
public:
  IsNeutrino(int pdg): PlotUtils::SignalConstraint<UNIVERSE>("IsNeutrino"),fpdg(pdg)
  {
  }
  
private:
  bool checkConstraint(const UNIVERSE& univ) const override
  {
    return univ.GetTruthNuPDG() == fpdg;
  }
  const int fpdg;
};


// IsCC (charge current interaction)
template <class UNIVERSE>
class IsCC: public PlotUtils::SignalConstraint<UNIVERSE>
{
public:
  IsCC(): PlotUtils::SignalConstraint<UNIVERSE>("IsCC")
  {
  }
  
private:
  bool checkConstraint(const UNIVERSE& univ) const override
  {
    return univ.GetCurrent() == 1;
  }
};


// MuonAngle 
template <class UNIVERSE>
class MuonAngle: public PlotUtils::SignalConstraint<UNIVERSE>
{
public:
  MuonAngle(const double angleMax): PlotUtils::SignalConstraint<UNIVERSE>(std::string("Muon Angle ") + std::to_string(angleMax)), fMax(angleMax*M_PI/180.)
  {
  }
  
private:
  bool checkConstraint(const UNIVERSE& univ) const override
  {
  return univ.GetThetalepTrue() <= fMax;
  }
  
  const double fMax;
};

// Muon Minimum Momentum
template <class UNIVERSE>
class PmuMin: public PlotUtils::SignalConstraint<UNIVERSE>
{
public:
  PmuMin(const double min): PlotUtils::SignalConstraint<UNIVERSE>(std::string("Pmu > ") + std::to_string(min)), fMin(min)
  {
  }
  
private:
  bool checkConstraint(const UNIVERSE& univ) const override
  {
  return univ.GetTruePmuGeV() >= fMin;
  }
  
  const double fMin;
};

// Muon Maximum Momentum
template <class UNIVERSE>
class PmuMax: public PlotUtils::SignalConstraint<UNIVERSE>
{
public:
  PmuMax(const double max): PlotUtils::SignalConstraint<UNIVERSE>(std::string("Pmu <= ") + std::to_string(max)), fmax(max)
  {
  }
  
private:
  bool checkConstraint(const UNIVERSE& univ) const override
  {
  return univ.GetTruePmuGeV() <= fmax;
  }
  
  const double fmax;
};

// Apothem
template <class UNIVERSE>
class Apothem: public PlotUtils::SignalConstraint<UNIVERSE>
{
public:
  Apothem(const double apothem): PlotUtils::SignalConstraint<UNIVERSE>(std::string("Apothem ") + std::to_string(apothem)), fApothem(apothem)
  {
  }
  
private:
  bool checkConstraint(const UNIVERSE& univ) const override
  {
  const auto vertex = univ.GetTrueVertex();
  return (fabs(vertex.y()) < fSlope*fabs(vertex.x()) + 2.*fApothem/sqrt(3.))
  && (fabs(vertex.x()) < fApothem);
  }
  
  const double fApothem;
  const double fSlope = -1./sqrt(3.); //A regular hexagon has angles of 2*M_PI/3, so I can find this is 1/tan(M_PI/3.)
};

// ZRange
template <class UNIVERSE>
class ZRange: public PlotUtils::SignalConstraint<UNIVERSE>
{
public:
  ZRange(const std::string& name, const double zMin, const double zMax): PlotUtils::SignalConstraint<UNIVERSE>(name), fMin(zMin), fMax(zMax)
  {
  }
  
private:
  bool checkConstraint(const UNIVERSE& univ) const override
  {
  return univ.GetTrueVertex().z() >= fMin && univ.GetTrueVertex().z() <= fMax;
  }
  
  const double fMin;
  const double fMax;
};

// --------------------- Quantities migrated from CCProtonPi0 NTuple TruthAnalysis ---------------------

// IsFidVol
template <class UNIVERSE>
class IsFidVol: public PlotUtils::SignalConstraint<UNIVERSE>
{
public:
  IsFidVol(): PlotUtils::SignalConstraint<UNIVERSE>("IsFidVol")
  {
  }
private:
  bool checkConstraint(const UNIVERSE& univ) const override
  {
    return univ.IsFidVol() == 1;
  }
};

// IsSignal
template <class UNIVERSE>
class IsSignal: public PlotUtils::SignalConstraint<UNIVERSE>
{
public:
  IsSignal(bool isSignal): PlotUtils::SignalConstraint<UNIVERSE>("IsSignal"), fsig(isSignal), nSignal(0)
  { 
  }
private:
  bool checkConstraint(const UNIVERSE& univ) const override
  {
    return univ.IsSignal() == fsig;
  }
  void countSignal(const UNIVERSE& univ)
  {
    if (univ.IsSignal() == 1) nSignal++;
  }
  const bool fsig;
  int nSignal;
};

}
#endif //CCPROTONPI0SIGNAL_H
