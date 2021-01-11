
//==============================================================================
// Adapted from the MINERvA Analysis Toolkit Example
// Compared to previous examples, this script is a more faithful attempt to
// perform a realistic event selection for the cross section pipeline.
//
// It attempts to use all/most of the new MAT tools as of 2020-09.
// These include:
// * CVUniverse + Histwrapper
// * PlotUtils::MacroUtil
// * Variable class
// * PlotUtils::Cutter
//
// This script follows the canonical event-looping structure:
// Setup (I/O, variables, histograms, systematics)
// Loop events
//   loop universes
//     make cuts
//     loop variables
//       fill histograms
// Plot and Save
//==============================================================================
#include <iostream>
#include <vector>
#include "utils/NuConfig.h"

#include "PlotUtils/MacroUtil.h"
#include "include/CVUniverse.h"
#include "include/Systematics.h"  // GetStandardSystematics
#ifndef __CINT__
#include "PlotUtils/Cut.h"
#include "PlotUtils/Cutter.h"
//#include "include/GetCCQEAntiNuCuts.h"
//#include "include/CCQEAntiNuSignal.h"
#include "include/GetCCProtonPi0Cuts.h"
#include "include/CCProtonPi0Signal.h"
#include "include/VariableWithMap.h"
#include "include/GetVariablesWithMap.h"
//TODO: Variable2D
#include "include/plotting_functions.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TError.h"
#endif

// needs to be global
double prescale;

// Forward declare my variables because we're hiding the header.
namespace CCProtonPi0 {
class VariableWithMap;
}

enum EDataMCTruth {kData, kMC, kTruth, kNDataMCTruthTypes};

//==============================================================================
// Fill Data or MC histograms
//==============================================================================
void FillData(std::string tag, CVUniverse* univ,
              std::vector<CCProtonPi0::VariableWithMap*> variables) {
  for (auto v : variables) {
    if (v->hasData[tag]){
    double fill_val = v->GetRecoValue(*univ, 0);
    v->m_selected_data.Fill(tag, univ, fill_val);
    }
  }
}

void FillMC(std::string tag, CVUniverse* univ, double weight,
            std::vector<CCProtonPi0::VariableWithMap*> variables) {
  for (auto v : variables) {
    if (v->hasMC[tag]){
    double fill_val = v->GetRecoValue(*univ, 0);
    v->m_selected_mc_reco.Fill(tag, univ, fill_val, weight);
    double true_val = v->GetTrueValue(*univ, 0);
    v->m_selected_mc_truth.Fill(tag, univ, true_val, weight);
    }
  }
}

void FillResponse(std::string tag,CVUniverse* univ, double weight,
            std::vector<CCProtonPi0::VariableWithMap*> variables) {
  for (auto v : variables) {
    if (v->hasResponse[tag]){
    double reco_val = v->GetRecoValue(*univ, 0);
    double true_val = v->GetTrueValue(*univ, 0);
    v->FillResponse(tag,univ,reco_val,true_val, weight);
    }
  }
}

void FillSignalTruth(std::string tag, CVUniverse* univ, double weight,
               std::vector<CCProtonPi0::VariableWithMap*> variables) {
  
  int run = univ->GetTrueRun();
  int subrun = univ->GetTrueSubRun();
  int gate = univ->GetTrueGate();
  //if (run == 123000 && subrun == 288) return;
   for (auto v : variables) {
     if (v->hasTruth[tag]){
      double true_val = v->GetTrueValue(*univ, 0);
      //if (univ->ShortName() == "cv" ) std::cout << v->GetName() << " " << univ->GetEventID() << " " << run << " " << subrun << " "  <<  " "  << true_val << " " << weight << std::endl;
      v->m_signal_mc_truth.Fill(tag,univ, true_val, weight);
    }
   }
}

//==============================================================================
// Loop and fill
//==============================================================================
#ifndef __CINT__ //for PlotUtils::cuts_t<>
void LoopAndFillEventSelection(std::string tag,
    const PlotUtils::MacroUtil& util,
    std::map<std::string, std::vector<CVUniverse*> > error_bands,
    std::vector<CCProtonPi0::VariableWithMap*>& variables,
    EDataMCTruth data_mc_truth,
    PlotUtils::Cutter<CVUniverse>& selection) {
  // Prepare loop
  MinervaUniverse::SetTruth(false);
  int nentries = -1;
  if (data_mc_truth == kData ){
    nentries = util.GetDataEntries();
  }
  else if (data_mc_truth == kMC ){
    nentries = util.GetMCEntries();
  }
  else{
    
      nentries = util.GetTruthEntries() ;
      MinervaUniverse::SetTruth(true);

  }

  std::cout << " starting loop " << data_mc_truth << " " << nentries << std::endl;
  // Begin entries loop
  for (int i = 0; i < nentries; i++) {
    if(data_mc_truth != kData) i+= prescale-1;
    if (i+1 % 1000 == 0) std::cout << (i / 1000) << "k " << std::endl;
    // Loop bands and universes
    for (auto band : error_bands) {
      std::vector<CVUniverse*> error_band_universes = band.second;
      for (auto universe : error_band_universes) {
        
        universe->SetEntry(i);

        // Process this event/universe
        double weight = data_mc_truth == kData ? 1. : universe->GetWeight();

        //TODO: Put the IsVertical() short circuit into PlotUtils::Cutter directly
        PlotUtils::detail::empty event;

        //=========================================
        // Fill
        //=========================================
        
        if(data_mc_truth == kMC){
          if(selection.isMCSelected(*universe, event, weight).all()
             && selection.isSignal(*universe)) {
            FillMC(tag, universe, weight, variables);
            FillResponse(tag, universe,weight,variables);
          }
        }
        else if (data_mc_truth == kTruth){
          
          if(selection.isEfficiencyDenom(*universe, weight)){
            FillSignalTruth(tag, universe, weight, variables);
          }
        }
        else{ //kData
          
          if(selection.isDataSelected(*universe, event).all()) { 
            FillData(tag, universe, variables);
          }
        }
      
      } // End universes
    } // End error bands
  } // End entries loop
}
#endif //__CINT__

//==============================================================================
// Main
//==============================================================================

int main(const int argc, const char *argv[] ) {

gROOT->ProcessLine("gErrorIgnoreLevel = kWarning;");
//++++++++++++++++++=  Initialization +++++++++++++++++++++++++
  
  std::string pl = "1A";
  if (argc > 1){
    pl = std::string(argv[1]);
  }
  else{
    std::cout << " arguments are:\n loop <config> [<prescale from config file>] " << std::endl;
    exit(0);
  }
  std::string configfilename(pl+".json");
  NuConfig config;
  config.Read(configfilename);
  
  if (argc > 2){
    prescale = std::stoi(argv[2]);
  }
  else{
    prescale = config.GetDouble("prescale");
  }
  bool LOCAL = config.GetBool("local");
  double kecut = config.GetDouble("protonKECutMeV");

//=========================================
  // MacroUtil (makes your anatuple chains)
  //=========================================
 
  const std::string mc_file_list(config.GetString("mcIn"));
  const std::string data_file_list(config.GetString("dataIn"));

  
  const std::string plist_string(config.GetString("playlist"));
  const std::string reco_tree_name("CCProtonPi0");
  const std::string truth_tree_name("Truth");
  const bool do_truth = true;
  const bool is_grid = false;

  PlotUtils::MacroUtil util(reco_tree_name, mc_file_list, data_file_list,
                    plist_string, do_truth, is_grid);
  
  //Data, MC reco, and Truth trees
  

  util.PrintMacroConfiguration("runEventLoop");

  //=========================================
  // Systematics
  // GetStandardSystematics in Systematics.h
  //=========================================
  std::map<std::string, std::vector<CVUniverse*> > mc_error_bands =
      systematics::GetStandardSystematics(util.m_mc,config);

  std::map<std::string, std::vector<CVUniverse*> > truth_error_bands =
      systematics::GetStandardSystematics(util.m_truth,config);
  
  std::cout << " mc error bands is " << mc_error_bands.size() << std::endl;
  std::cout << " truth error bands is " << truth_error_bands.size() << std::endl;
  // Setup systematics-related constants
  // TODO These should be handled by PU::MacroUtil
  PlotUtils::MinervaUniverse::SetNonResPiReweight(true);
  
  PlotUtils::MinervaUniverse::SetNuEConstraint(false); //Needs to be on to match Dan's ME 2D inclusive analysis
  PlotUtils::MinervaUniverse::SetPlaylist(plist_string);
  
  PlotUtils::MinervaUniverse::SetAnalysisNuPDG(14);

  // If we're doing data we need a single central value universe.
  // TODO Take care of this in PU::MacroUtil.
  CVUniverse* data_universe = new CVUniverse(util.m_data);
  std::vector<CVUniverse*> data_band = {data_universe};
  std::map<std::string, std::vector<CVUniverse*> > data_error_bands;
 
  data_error_bands["cv"] = data_band;

  //Selection Criteria
   
  
  PlotUtils::cuts_t<CVUniverse> noSidebands;
  PlotUtils::constraints_t<CVUniverse> signal, phaseSpace;
  std::map<std::string,PlotUtils::Cutter<CVUniverse> *> selectionCriteria;
  selectionCriteria["CCProtonPi0"] = new PlotUtils::Cutter<CVUniverse> (reco::GetCCProtonPi0Cuts<CVUniverse>(),std::move(noSidebands),truth::GetCCProtonPi0Signal<CVUniverse>(kecut),truth::GetCCProtonPi0PhaseSpace<CVUniverse>());
  selectionCriteria["CCProtonPi0Bkg"] = new PlotUtils::Cutter<CVUniverse> (reco::GetCCProtonPi0Cuts<CVUniverse>(),std::move(noSidebands),truth::GetCCProtonPi0Bkg<CVUniverse>(kecut),truth::GetCCProtonPi0PhaseSpace<CVUniverse>());
  
  selectionCriteria["Data"] = new PlotUtils::Cutter<CVUniverse> (reco::GetCCProtonPi0Cuts<CVUniverse>(),
                                                  std::move(noSidebands),
                                                  truth::GetCCProtonPi0Signal<CVUniverse>(kecut),
                                                  truth::GetCCProtonPi0PhaseSpace<CVUniverse>());

  //=========================================
  // Get variables and initialize their hists
  //=========================================
 
  
  std::vector<std::string> tags;
  
  for (auto cuts:selectionCriteria){
    tags.push_back(cuts.first);
  }
  std::vector<std::string> vars = config.GetStringVector("AnalyzeVariables");
  
  std::vector<CCProtonPi0::VariableWithMap*> variables = GetVariablesWithMap(vars,tags);
  
  // associate different datasets with tags so you don't have to loop over everything for everything
  
  std::vector<std::string> mctags = {"CCProtonPi0","CCProtonPi0Bkg"};  //bkg only needed for recontructed MC
  std::vector<std::string> datatags = {"Data"};
  std::vector<std::string> truthtags = {"CCProtonPi0"};
  std::vector<std::string> responsetags = {"CCProtonPi0"};
  for (auto v : variables) {
    v->InitializeMCHistograms(mc_error_bands,mctags);
    v->InitializeDataHistograms(data_error_bands,datatags);
    v->AddMCResponse(responsetags);
    v->InitializeTruthHistograms(truth_error_bands,truthtags);
  }

  for (auto tag:datatags){
  //=========================================
  // Entry loop and fill
  //=========================================
    std::cout << "Loop and Fill Data for " << tag << "\n" ;
    LoopAndFillEventSelection(tag, util, data_error_bands, variables, kData, *selectionCriteria[tag]);
    
    std::cout << "\nCut summary for Data:" <<  tag << "\n" << *selectionCriteria[tag] << "\n";
    selectionCriteria[tag]->resetStats();
  }/*
  for (auto tag:mctags){
    std::cout << "Loop and Fill MC Reco  for " <<  tag << "\n";
    LoopAndFillEventSelection(tag, util, mc_error_bands, variables, kMC, *selectionCriteria[tag]);
  
    std::cout << "\nCut summary for MC Reco:" <<  tag << "\n" << *selectionCriteria[tag] << "\n";
    selectionCriteria[tag]->resetStats();
  }*/
  for (auto tag:truthtags){
    std::cout << "Loop and Fill MC Truth  for " <<  tag << "\n";
    LoopAndFillEventSelection(tag, util, truth_error_bands, variables, kTruth, *selectionCriteria[tag]);
    std::cout << "\nCut summary for MC Truth:" <<  tag << "\n";
    std::cout < *selectionCriteria[tag];
    std::cout << "\n";
    selectionCriteria[tag]->resetStats();
  }

  for (auto tag:mctags){
    std::cout << "Loop and Fill MC Reco  for " <<  tag << "\n";
    LoopAndFillEventSelection(tag, util, mc_error_bands, variables, kMC, *selectionCriteria[tag]);

    std::cout << "\nCut summary for MC Reco:" <<  tag << "\n" << *selectionCriteria[tag] << "\n";
    selectionCriteria[tag]->resetStats();
  }

    //std::cout << "\nCut summary for MC:\n" << selectionCriteria << "\n";

  // Sync CV hists
  
  for (auto v : variables) v->SyncAllHists();

  //=========================================
  // Plot
  //=========================================
  std::cout << "Done filling. Begin plotting.\n";
  
  for (auto tag:tags){
    for (auto v : variables) {
    if (v->hasMC[tag] && v->hasData[tag]){
      PlotDataMCAndError(v->m_selected_data.GetHist(tag), v->m_selected_mc_reco.GetHist(tag),
                        (util.m_data_pot/util.m_mc_pot*prescale), v->GetName()+"_eff");
      PlotErrorSummary(v->m_selected_mc_reco.GetHist(tag), v->GetName()+"_eff");
    }
    if (v->hasMC[tag] && v->hasTruth[tag]){PlotDataMCAndError(v->m_selected_mc_truth.GetHist(tag), v->m_signal_mc_truth.GetHist(tag), 1, v->GetName()+"_ratio");
    }
    
   // PlotStacked(v->m_selected_data.hist,v->m_selected_mc_by_channel.GetHistArray(), (util.m_data_pot/util.m_mc_pot)*prescale, v->m_selected_mc_by_channel.GetName(), v->m_selected_mc_by_channel.GetName());
  }
  }

 
  
  std::string outname = config.GetString("outRoot")+"_"+pl+".root";
  TFile* out=TFile::Open(outname.c_str(),"RECREATE");
  for (auto v : variables){
    v->WriteAllHistogramsToFile(*out);
  }
  
  PlotUtils::MnvH1D* h_pot = new PlotUtils::MnvH1D("POT_summary","data, mc, Prescaled mc", 3, 0., 3.);
  h_pot->Fill(0.5, util.m_data_pot);
  h_pot->Fill(1.5, util.m_mc_pot);
  h_pot->Fill(2.5, util.m_mc_pot/prescale);  // compensate for prescale
  h_pot->Print();
  h_pot->Write();
  out->Close();
  std::cout << "Success" << std::endl;
}
