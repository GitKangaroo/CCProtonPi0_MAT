// ========================================================================
// Base class for an un-systematically shifted (i.e. CV) universe.
// Implement "Get" functions for all the quantities that you need for your
// analysis.
//
// This class inherits from PU::sUniverse. PU::DCVU may already define
// your "Get" functions the way you want them. In that case, you  don't need to
// re-write them here.
// ========================================================================
#ifndef CVUNIVERSE_H
#define CVUNIVERSE_H

#include <iostream>

#include "PlotUtils/ChainWrapper.h"
#include "PlotUtils/MinervaUniverse.h"
#include "PlotUtils/PhysicsVariables.h"
#include "TVector3.h"
#include "PlotUtils/GeantHadronSystematics.h"


class CVUniverse : public PlotUtils::MinervaUniverse {
public:
  #include "PlotUtils/SystCalcs/WeightFunctions.h"
  #include "PlotUtils/SystCalcs/MuonFunctions.h"
  #include "PlotUtils/SystCalcs/TruthFunctions.h"
  #include "PlotUtils/SystCalcs/RecoilEnergyFunctions.h"
  //#include "PlotUtils/SystCalcs/ProtonKECutFunctions.h"
  
  double MeVGeV = 0.001; // lazy conversion from MeV to GeV before filling histos
  bool m_useNeutronCVReweight = true;

  // ========================================================================
  // Constructor/Destructor
  // ========================================================================
  
  // fake a default constructor
  
  CVUniverse(){PlotUtils::MinervaUniverse();};
  
  CVUniverse(PlotUtils::ChainWrapper* chw, double nsigma = 0)
  : PlotUtils::MinervaUniverse(chw, nsigma) {}

 
  virtual ~CVUniverse() {}
  
  
  // Geant weight stuff overriden here as it doesn't have its own header yet.
  
  
  void InitGeantWeighter( std::string project_name,
                             double min_z, double max_z, double apothem,
                             bool useNeutronCVReweight, bool useElastics)
     {
       std::string projectname = PlotUtils::weight_hadron<PlotUtils::ChainWrapper*>( dynamic_cast<PlotUtils::ChainWrapper*>(m_chw), project_name, min_z, max_z, apothem, useNeutronCVReweight, useElastics ).getProjectName();
       std::cout<<"Filename "<<PlotUtils::weight_hadron<PlotUtils::ChainWrapper*>().makefullfilename("",projectname)<<std::endl;
       m_useNeutronCVReweight = useNeutronCVReweight;
     }
     virtual double GetGeantHadronWeight() const
     {
       return m_useNeutronCVReweight ? PlotUtils::weight_hadron<PlotUtils::ChainWrapper*>().reweightNeutronCV( *this ) : 1;
     }

  // ========================================================================
  // Get Weight
  // ========================================================================
  //
  virtual double GetWeight() const {
    double wgt_flux_and_cv = 1., wgt_genie = 1., wgt_2p2h = 1.;
    double wgt_rpa = 1., wgt_mueff = 1.;
    double wgt_geant = 1.0;
    // flux + cv
    wgt_flux_and_cv = GetFluxAndCVWeight();

    // genie
    wgt_genie = GetGenieWeight();

    // 2p2h
    wgt_2p2h = GetLowRecoil2p2hWeight();

    // rpa
    wgt_rpa = GetRPAWeight();
    
//    if (!IsTruth()){ wgt_geant = PlotUtils::GetGeantHadronWeight();
//     // if (wgt_geant != 1.0) std::cout << ShortName() << wgt_geant << std::endl;
//    }

    // MINOS muon tracking efficiency
    if (!IsTruth() && IsMinosMatchMuon()) wgt_mueff = GetMinosEfficiencyWeight();

    return wgt_flux_and_cv * wgt_genie * wgt_2p2h * wgt_rpa * wgt_mueff * wgt_geant;
  }

  // ========================================================================
  // Write a "Get" function for all quantities access by your analysis.
  // For composite quantities (e.g. Enu) use a calculator function.
  //
  // Currently PlotUtils::MinervaUniverse may already define some functions
  // that you want. The future of these functions is uncertain. You might want
  // to write all your own functions for now.
  // ========================================================================
 
  virtual int GetEventID() const{
    return GetInt("eventID");
  } 
  virtual int GetMultiplicity() const{
    return GetInt("nTracks_Close");
  }
  
  virtual int GetDeadTime() const{
    return GetInt("phys_n_dead_discr_pair_upstream_prim_track_proj") ;
  }
  // ----------------------- Analysis-related Variables ------------------------

  virtual double GetEnuHadGeV() const {
  // virtual double GetEnuGeV() const {
    return GetEmuGeV()+GetHadronEGeV();
  }

  virtual double GetTrueEnuGeV() const {
    return GetDouble("mc_incomingE")*MeVGeV;
  }

  virtual double GetEnuCCQEGeV() const {
    int charge = GetAnalysisNuPDG() > 0? 1:-1;
    //std::cout << " before try" << GetEmu() <<  std::endl;
    double val = PlotUtils::nuEnergyCCQE( GetEmu(), GetPmu(), GetThetamu(), charge ); // un-PlotUtil's this?
    //std::cout << " try to getEnuCCQE" << val << std::endl;
    return val*MeVGeV;
  } // both neutrino and antinu

  virtual double GetTrueEnuCCQEGeV() const {
    int charge = GetAnalysisNuPDG() > 0? 1:-1;
    //std::cout << " before try" << GetEmu() <<  std::endl;
    double val = PlotUtils::nuEnergyCCQE( GetElepTrue(), GetPlepTrue(), GetThetalepTrue(), charge ); // un-PlotUtil's this?
    //std::cout << " try to getEnuCCQE" << val << std::endl;
    return val*MeVGeV;
  } // may be a better way to implement this

  virtual double GetQ2QEGeV() const {
    const double q2min = 0.001;
    int charge = GetAnalysisNuPDG() > 0? 1:-1;
    if (GetEnuCCQEGeV()<=0.0)return 0.0;
    //std::cout<<"CVUniverse::GetQ2QE Cannot calculate neutrino energy "<<std::endl;
    else{
    // Q2 = 2.0*GetEnuCCQE() * ( GetEmu() - GetPmu() * cos( GetThetamu() ) ) - pow( MinervaUnits::M_mu, 2 );
      double q2 = PlotUtils::qSquaredCCQE( GetEmu(), GetPmu(), GetThetamu(), charge )*MeVGeV*MeVGeV;
    //  if (q2 < q2min)q2 =  q2min;
      return q2;
    }
    // return Q2;
  }

  virtual double GetTrueQ2QEGeV() const {
      const double q2min = 0.001;
      int charge = GetAnalysisNuPDG() > 0? 1:-1;
      double q2 = PlotUtils::qSquaredCCQE( GetElepTrue(), GetPlepTrue(), GetThetalepTrue(), charge )*MeVGeV*MeVGeV;
     // if (q2 < q2min)q2 =  q2min;
      return q2;
    }

  virtual double GetLog10Q2QEGeV() const {
    return std::log10( GetQ2QEGeV() );
  }

  virtual double GetTrueLog10Q2QEGeV() const {
    return std::log10( GetTrueQ2QEGeV() );
  }


  // ------------------------------ Muon Variables -----------------------------

  virtual double GetEmuGeV() const {
    return std::sqrt(GetPmu()*GetPmu() + pow( MinervaUnits::M_mu, 2 ))*MeVGeV;
  }

  virtual double GetTrueEmuGeV() const {
    return GetElepTrue(); // not sure if this is right
  }

  virtual double GetPmuGeV() const {
    return GetPmu()*MeVGeV;
  }

  virtual double GetTruePmuGeV() const {
    return GetPlepTrue()*MeVGeV;
  }

  virtual double GetPparMuGeV() const {
  // virtual double GetPparMuGeV(bool istruth=false) const {
    return GetPmuGeV()*std::cos( GetThetamu() );
  }

  virtual double GetTruePparMuGeV() const {
    return GetTruePmuGeV()*std::cos( GetTrueThetamu() );
  }

  virtual double GetPperpMuGeV() const {
  // virtual double GetPperpMu(bool istruth=false) const {
    return GetPmuGeV()*std::sin( GetThetamu() );
  }

  virtual double GetTruePperpMuGeV() const {
    return GetTruePmuGeV()*std::sin( GetTrueThetamu() );
  }

  virtual double GetTrueThetaXmu() const {
    TVector3 p3lep( GetVecElem("mc_primFSLepton",0), GetVecElem("mc_primFSLepton",1), GetVecElem("mc_primFSLepton",2) );
      p3lep.RotateX(MinervaUnits::numi_beam_angle_rad);
      double px =  p3lep[0];
      double pz =  p3lep[2];
      //std::cout << "tx " << px << " " << pz << std::endl;
      return std::atan2(px,pz);
  }
  virtual double GetTrueThetaYmu() const {
    TVector3 p3lep( GetVecElem("mc_primFSLepton",0), GetVecElem("mc_primFSLepton",1), GetVecElem("mc_primFSLepton",2) );
      p3lep.RotateX(MinervaUnits::numi_beam_angle_rad);
      double py =  p3lep[1];
      double pz =  p3lep[2];
      //std::cout << "tx " << px << " " << pz << std::endl;
      return std::atan2(py,pz);
  }

  virtual double GetTrueThetamu() const {
     return GetThetalepTrue();
   }


  // ----------------------------- Proton Variables ----------------------------

  virtual double GetHadronEGeV() const {
    return (GetCalRecoilEnergyGeV() + GetNonCalRecoilEnergyGeV());
  }


  // ----------------------------- Recoil Variables ----------------------------

  
  virtual double GetCalRecoilEnergy() const {
    bool neutrinoMode = GetAnalysisNuPDG() > 0;
    if(neutrinoMode) return (GetDouble("CCProtonPi0_hadron_recoil_CCInc")/*"nonvtx_iso_blobs_energy")+GetDouble("dis_id_energy")*/); // several definitions of this, be careful
    else {
      if (GetVecDouble("recoil_summed_energy").size()==0) return -999.; // protect against bad input,
      return (GetVecDouble("recoil_summed_energy")[0] - GetDouble("recoil_energy_nonmuon_vtx100mm"));
    }
  }

  virtual double GetCalRecoilEnergyGeV() const {
    return GetCalRecoilEnergy()*MeVGeV;
  }
  
  
  virtual double GetNonCalRecoilEnergy() const {
    // not certain why I want to implement this but there ya go.
    return 0;
  }
  
  virtual double GetNonCalRecoilEnergyGeV() const {
    return GetNonCalRecoilEnergy()*MeVGeV;
  }

  virtual double GetRecoilEnergyGeV() const {
    // return GetCalRecoilEnergy();
    return GetRecoilEnergy()*MeVGeV;
  }

  // virtual double GetTrueRecoilEnergyGeV() const {
  //   // don't know if this needs to exist
  //   // return 0;
  //   return GetTrueQ0GeV();
  // }

  virtual double GetTrueQ0GeV() const {
    static std::vector<double> mc_incomingPartVec;
    static std::vector<double> mc_primFSLepton;
    mc_incomingPartVec = GetVecDouble("mc_incomingPartVec");
    mc_primFSLepton = GetVecDouble("mc_primFSLepton");
    double q0 = mc_incomingPartVec[3] - mc_primFSLepton[3];
    return q0*MeVGeV;
  }
  virtual double GetTrueQ3GeV() const {
    static std::vector<double> mc_incomingPartVec;
    static std::vector<double> mc_primFSLepton;
    mc_incomingPartVec = GetVecDouble("mc_incomingPartVec");
    mc_primFSLepton = GetVecDouble("mc_primFSLepton");
    double px = mc_primFSLepton[0] - mc_incomingPartVec[0];
    double py = mc_primFSLepton[1] - mc_incomingPartVec[1];
    double pz = mc_primFSLepton[2] - mc_incomingPartVec[2];
    double q3 = std::sqrt( px*px + py*py + pz*pz );
    return q3*MeVGeV;
  }
  virtual double GetTrueQ2GeV() const {
    double q3 = GetTrueQ3GeV();
    double q0 = GetTrueQ0GeV();
    return q3*q3 - q0*q0;
  }

  // ----------------------------- Pi0 Variables -----------------------------

  virtual double GetRecoPi0Mass() const {
    return GetDouble("pi0_invMass")*MeVGeV;
  }
  // Hard coded pi0 mass,Need to change it in future
  virtual double GetTruePi0Mass() const {
    return 0.13497;
  }
  // ----------------------------- Other Variables -----------------------------

  //  virtual double GetWgenie() const { return GetDouble("mc_w"); }
  virtual int GetMCIntType() const {
    return GetInt("mc_intType");
  }

  virtual int GetTruthNuPDG() const {
    return GetInt("mc_incoming");
  }

  virtual int GetCurrent() const {
    return GetInt("mc_current");
  }

  virtual double GetTpiGeV( const int hadron ) const {
    double TLA = GetVecElem("hadron_track_length_area", hadron);
    return (2.3112 * TLA + 37.03)*MeVGeV; // what are these numbers
  }


  // --------------------- Quantities only needed for cuts ---------------------
  // Although unlikely, in principle these quanties could be shifted by a
  // systematic. And when they are, they'll only be shifted correctly if we
  // write these accessor functions.
  virtual bool IsMinosMatchMuon() const {
    return GetInt("muon_hasMinosMatchTrack") == 1;  // does this flip for neutrino?
  } // This isn't even used anymore, there's something else. This is left over from Amit's analysis

  // No information for this variable in CCProtonPi0 Ntuples (declared but not filled)
  virtual int GetNuHelicity() const {
    return GetInt("CCProtonPi0_nuHelicity");
  }
  
  // Muon Charge information to check if it's a neutrino interaction
  virtual int GetMuonCharge() const{
    return GetInt("muon_charge");
  }

  virtual int GetTDead() const {
    return GetInt("tdead");
  } // Dead electronics, a rock muon removal technique. Amit's analysis doesn't
  // have that cut most likely. Not even in that tuple, only reason it survives is bc it's not called anymore.
  // Can't find it? Just a hard coded constant in CCQENuCutsNSF.h

  // These cuts are already made in the CCProtonPi0 AnaTuple, may be unnecessary
  ROOT::Math::XYZTVector GetVertex() const
  {
  ROOT::Math::XYZTVector result;
  result.SetCoordinates(GetVec<double>("vtx").data());
  return result;
  }

  ROOT::Math::XYZTVector GetTrueVertex() const
  {
  ROOT::Math::XYZTVector result;
  result.SetCoordinates(GetVec<double>("mc_vtx").data());
  return result;
  }


  // Some stuff Heidi added to test out some issues with the NTuples

  virtual int GetRun() const {
    return GetInt("ev_run");
  }

  virtual int GetSubRun() const {
    return GetInt("ev_subrun");
  }

  virtual int GetGate() const {
    return GetInt("ev_gate");
  }

  virtual int GetTrueRun() const {
     return GetInt("mc_run");
   }

  virtual int GetTrueSubRun() const {
     return GetInt("mc_subrun");
   }

   int GetTrueGate() const {  // mot certain if this is stored
     return GetInt("mc_nthEvtInFile")+1;
   }

   // --------------------- Quantities migrated from CCProtonPi0 NTuple TruthAnalysis ---------------------
   
   // IsFidVol
   virtual int IsFidVol() const {
     return GetInt("truth_isFidVol");
   }
   // IsSignal
   virtual int IsSignal() const {
     return GetInt("truth_isSignal");
   }

   // getTrueBkgWithPi0
   virtual std::vector<bool> GetTrueBkgWithPi0() const {
     std::vector<bool> TrueBkgWithPi0;
     TrueBkgWithPi0.push_back(GetInt("truth_isBckg_NoPi0"));
     TrueBkgWithPi0.push_back(GetInt("truth_isBckg_SinglePi0"));
     TrueBkgWithPi0.push_back(GetInt("truth_isBckg_MultiPi0"));
     return TrueBkgWithPi0;
   }


 
};
#endif
