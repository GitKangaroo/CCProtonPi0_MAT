//==============================================================================
//In this file several CCProtonPi0 cuts are defined.
//
//Each cut is a class that inherits from PlotUtils::Cut.
//
//At minimum, cuts must have a name and they must override the checkCut
//function.
//
//You cut also takes two template parameters:
//* UNIVERSE template parameter: Cuts are built on Universe objects. Plug
//in your CVUniverse as a template parameter so the cut can access your
//branches and do so correctly within a systematic universe.
//
// * EVENT template parameter: sometimes cuts need to do more than just return
//a bool when you call the checkCut function. The event object is a
//user-specifiable object that can hold onto information that is learned within
//checkCut. E.g. A HasMichel cut determines which tracks are potential pion
//tracks. See CCPionCuts.h for a fleshed out example making use of EVENT.
//==============================================================================

// Original Authors: Andrew Olivier aolivier@ur.rochester.edu and Heidi Schellman heidi.schellman@oregonstate.edu

// CCProtonPi0Cuts.h was created by Kang Yang kang.yang@physics.ox.ac.uk

//Public-use reco::Cuts Classes here defined:
// -------------------------------------------
// HasMINOSMatch
// MuonMomentumMin
// MuonMomentumMax
// RecoilMax
// Pi0MassMin (useless?)
// NoDeadtime
// IsNeutrino
// MaxMuonAngle
// MaxMultiplicity
// ZRange
// Apothem

#include "PlotUtils/Cut.h"
#include "include/CVUniverse.h"

#ifndef CCPROTONPI0CUTS_H
#define CCPROTONPI0CUTS_H

namespace reco {
//============================================================================
//Example 1: The simplest cut example. Just derive from Cut<> base class.
//============================================================================
template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
class HasMINOSMatch: public PlotUtils::Cut<UNIVERSE, EVENT>
{
public:
  // Constructor
  HasMINOSMatch(): PlotUtils::Cut<UNIVERSE, EVENT>("Has MINOS Match") {}
  
private:
  // THE cut function
  bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
  {
  // Call a CVUniverse member function to make the cut (see funtion definitions in include/CVUniverse.h)
  //std::cout << " muon match " <<univ.IsMinosMatchMuon() << std::endl;
  return univ.IsMinosMatchMuon() == 1;
  }
};

//============================================================================
//Example 2: Many cuts require variables to be above, below, or at some value.
//To simply help avoid typos and reduce a little typing, use Minimum,
//Maximum, and IsSame helper templates.
//In this case the muon energy is required to be at least X, at minimum.
//The specific value of X gets set when you instantiate this cut (in
//getCCqeCuts).
//============================================================================

// Muon Minimum Momentum Cut 
template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
using MuonMomentumMin = PlotUtils::Minimum<UNIVERSE, double, &UNIVERSE::GetPmuGeV, EVENT>;

// Muon Maximum Momentum Cut
template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
using MuonMomentumMax = PlotUtils::Maximum<UNIVERSE, double, &UNIVERSE::GetPmuGeV, EVENT>;

// Maximum Recoil Energy Cut
//==========================================================================
template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
using RecoilMax = PlotUtils::Maximum<UNIVERSE, double, &UNIVERSE::GetRecoilEnergyGeV, EVENT>;

//Minimum Pi0 Invariant Mass
template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
using Pi0MassMin = PlotUtils::Minimum<UNIVERSE, double, &UNIVERSE::GetRecoPi0Mass, EVENT>;

//============================================================================
//Example 3: The first, commented-out dead time cut inherits directly from Ben's
//CVUniverse class. (In form, it resembles case 1, by the way.) While this is
//valid, it means that no one else can use this cut (because they don't also
//use my CVUniverse).
//============================================================================
//============================================================================
//Example 3 contd: It's better to make the cut a template so that any analyzer's
//CVUniverse will work with it.
//Furthermore, look how much space we save by using the Maximum helper
//template.
//============================================================================

template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
using NoDeadtime = PlotUtils::Maximum<UNIVERSE, int, &UNIVERSE::GetDeadTime>;

//============================================================================
//Example 4: Another case like Example 1.
//Physics-wise, we're making sure that the event is from a neutrino as
//opposed to an antineutrino.
//============================================================================

// I used muon_charge defined in CCProtonPi0 NTuples, instead of using NuHelicity (Need to confirm if this is okay)
// Actually all my events have muon charge -1 (no need for this cuts ?)
template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
class IsNeutrino: public PlotUtils::Cut<UNIVERSE, EVENT>
{
public:
  IsNeutrino(const bool isneutrino): PlotUtils::Cut<UNIVERSE, EVENT>("Muon Charge Sign"),fisneutrino(isneutrino) {}
  
private:
  bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
  {
  //if (fisneutrino) return univ.GetNuHelicity()==1;
  //return univ.GetNuHelicity()==2;
  if (fisneutrino) return univ.GetMuonCharge()==-1;
  return univ.GetMuonCharge()==1;
  }
  const bool fisneutrino;
};

//==========================================================================

/*!!!!!Not Clear what is the recoil cut in CCProtonPi0, need to check the old codes.

double f_recoil_cut(double Q2) {
  
  double result = 0;
  double offset=0.0;
  
  
  if( Q2 < 0.175 )
    return result = 0.08;
  else if( Q2 < 1.4 )
    return result = 0.03 + 0.3*Q2;
  else
    return result = 0.45; //antinu
  
};

template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
class GoodRecoil: public PlotUtils::Cut<UNIVERSE, EVENT>
{
public:
  GoodRecoil(): PlotUtils::Cut<UNIVERSE, EVENT>("RecoilCut") {}
  
private:
  bool checkCut(const UNIVERSE& univ, EVENT& /*evt*//*) const override
  {
  return univ.GetRecoilEnergyGeV()<f_recoil_cut(univ.GetQ2QEGeV());
  }
};
*/
//==========================================================================

// Muon Maximum Angle Cut 
template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
class MaxMuonAngle: public PlotUtils::Cut<UNIVERSE, EVENT>
{
public:
  MaxMuonAngle(const double max): PlotUtils::Cut<UNIVERSE, EVENT>(std::string("Muon Theta < ") + std::to_string(max) + std::string("degrees")), fMax(max*M_PI/180.) {}
  
private:
  bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
  {
  return univ.GetThetamu() < fMax;
  }
  
  const double fMax; //radians
};

//==========================================================================

// Multiplicity Cut
template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
class MaxMultiplicity: public PlotUtils::Cut<UNIVERSE, EVENT>
{
public:
  MaxMultiplicity( const int max): PlotUtils::Cut<UNIVERSE, EVENT>(std::string("Multiplicity <= ") + std::to_string(max)), fMax(max) {}
  
private:
  bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
  {
  return univ.GetMultiplicity() <= fMax;
  }
  const double fMax;
};


//==========================================================================

// Z Range Cut (need to check for CCProtonPi0)
template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
class ZRange: public PlotUtils::Cut<UNIVERSE, EVENT>
{
public:
  ZRange(const std::string& name, const double zMin, const double zMax): PlotUtils::Cut<UNIVERSE, EVENT>(name), fMin(zMin), fMax(zMax)
  {
  }
  
private:
  bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
  {
  return univ.GetVertexZ() >= fMin && univ.GetVertexZ() <= fMax;
  }
  
  const double fMin;
  const double fMax;
};

//==========================================================================

// Apothem Cut (need to check for CCProtonPi0)
template <class UNIVERSE, class EVENT = PlotUtils::detail::empty>
class Apothem: public PlotUtils::Cut<UNIVERSE, EVENT>
{
public:
  Apothem(const double apothem): PlotUtils::Cut<UNIVERSE, EVENT>(std::string("Apothem ") + std::to_string(apothem)), fApothem(apothem)
  {
  }
  
private:
  bool checkCut(const UNIVERSE& univ, EVENT& /*evt*/) const override
  {
  const auto vertex = univ.GetVertex();
  return (fabs(vertex.y()) < fSlope*fabs(vertex.x()) + 2.*fApothem/sqrt(3.))
  && (fabs(vertex.x()) < fApothem);
  }
  
  const double fApothem;
  const double fSlope = -1./sqrt(3.); //A regular hexagon has angles of 2*M_PI/3, so I can find this is 1/tan(M_PI/3.)
};
  }
#endif //CCPROTONPI0CUTS_H
