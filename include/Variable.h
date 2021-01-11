#ifndef VARIABLE_H
#define VARIABLE_H

#include "CVUniverse.h"
#include "PlotUtils/HistFolio.h"
#include "PlotUtils/HistWrapper.h"
#include "MinervaUnfold/MnvResponse.h"


#ifndef __CINT__  // CINT doesn't know about std::function
#include "PlotUtils/VariableBase.h"
#endif  // __CINT__

namespace CCProtonPi0 {

class Variable : public PlotUtils::VariableBase<CVUniverse> {
private:
  typedef PlotUtils::HistWrapper<CVUniverse> HW;
  typedef MinervaUnfold::MnvResponse RESPONSE;
  typedef PlotUtils::MnvH1D MH1D;
  typedef PlotUtils::HistFolio<PlotUtils::MnvH1D> FOLIO;
  
  
public:
  //=======================================================================================
  // CTOR
  //=======================================================================================
  template <class... ARGS>
  Variable(ARGS... args) : PlotUtils::VariableBase<CVUniverse>(args...) {}
  
  //=======================================================================================
  // DECLARE NEW HISTOGRAMS
  //=======================================================================================
  // Histwrappers -- selected mc, selected data
  HW m_selected_mc_reco;
  HW m_selected_mc_truth;
  HW m_signal_mc_truth;
  HW m_selected_data;
  std::string m_units;
  bool hasData = false;
  bool hasMC = false;
  bool hasTruth = false;
  bool hasResponse = false;
  RESPONSE* m_response;
  // helpers for response
  std::map< std::string, std::vector<CVUniverse*> > m_universes;// map to get universe from name and index
  std::map< CVUniverse* , int> m_map; // map to get name and index from Universe;
  // Histofolio to categorize MC by interaction channel
  FOLIO m_selected_mc_by_channel;
  //FOLIO m_selected_mc_truth_by_channel;
  inline virtual std::string GetUnits(){return m_units;};
  inline virtual void SetUnits(std::string units){m_units=units;};
  
  //=======================================================================================
  // INITIALIZE ALL HISTOGRAMS
  //=======================================================================================
  template <typename T>
  void InitializeMCHistograms(T univs) {
    hasMC = true;
    std::vector<double> bins = GetBinVec();
    const char* name = GetName().c_str();
    
    m_selected_mc_reco = HW(Form("selected_mc_reco_%s", name), (GetName()+";"+m_xaxis_label).c_str(), GetNBins(), bins, univs);
    
    //Interaction types from https://nusoft.fnal.gov/minerva/minervadat/software_doxygen/HEAD/MINERVA/classMinerva_1_1GenMinInteraction.html
    std::vector<PlotUtils::NamedCategory<int>> genieCategories = {{1, "QE"}, {2, "RES"}, {3, "DIS"}, {4, "Coherent Pi"}, {8, "MEC"}};
    m_selected_mc_by_channel = FOLIO(genieCategories, std::string("selected_mc_in_") + name + "by_channel", GetNBins(), bins.data());
    m_selected_mc_truth = HW(Form("selected_mc_truth_%s", name), (GetName()+";"+m_xaxis_label).c_str(), GetNBins(), bins, univs);
  }
  
  template <typename T>
  void InitializeTruthHistograms(T univs) {
    hasTruth = true;
    std::vector<double> bins = GetBinVec();
    const char* name = GetName().c_str();
    
 
    m_signal_mc_truth = HW(Form("signal_mc_truth_%s", name), (GetName()+";"+m_xaxis_label).c_str(), GetNBins(), bins, univs);
    
  }
  
  //========== Add Response =================
  
  template <typename T>
  void InitializeMCResponse(T univs) {
    std::vector<double> bins = GetBinVec();
    const char* name = GetName().c_str();
    std::string resp_name=string("response_")+name;
    
    
    std::map<std::string, int> response_bands;
    // make response talk to universes
    
    for (auto band : univs){
      std::vector<CVUniverse*> error_band_universes = band.second;
      int size = band.second.size();
      assert(size>0);
      std::string realname = (band.second)[0]->ShortName();
      response_bands[realname] = size;
      // store the index for each universe in the error Band
      for (int i = 0; i < size; i++){
        m_map[band.second[i]] = i;
      }
    }
    
    hasResponse = true;
    m_response = new RESPONSE(resp_name.c_str(), resp_name.c_str(), GetNBins(), bins.data(), GetNBins(), bins.data(), response_bands);
  }
  
  template <typename T>
  void InitializeDataHistograms(T univs) {
    hasData = true;
    std::vector<double> bins = GetBinVec();
    const char* name = GetName().c_str();
    m_selected_data = HW(Form("selected_data_%s", name), (GetName()+";"+m_xaxis_label).c_str(), GetNBins(), bins, univs);
  }
  
  //=======================================================================================
  // WRITE ALL HISTOGRAMS
  //=======================================================================================
  void WriteAllHistogramsToFile(TFile& f) const {
    f.cd();
    
    
    // selected mc reco
    std::cout << " write out mc " << std::endl;
    if(hasMC) {
      m_selected_mc_reco.hist->Write();
          std::cout << " write out mc by channel" << std::endl;
      m_selected_mc_by_channel.WriteToFile(f);
      m_selected_mc_truth.hist->Write();
    }
    
    if(hasResponse){
      std::cout << " write out response " << std::endl;
      PlotUtils::MnvH2D* h_migration;
      PlotUtils::MnvH1D* h_reco;
      PlotUtils::MnvH1D* h_truth;
      std::cout << " GetMigrationObjects will now complain because I passed it pointers to uninitiated MnvH2D/1D to fill please ignore" << std::endl;
      m_response->GetMigrationObjects( h_migration, h_reco, h_truth);
      if (h_reco->GetEntries() > 0){
        h_migration->Write();
        h_reco->Write();
        h_truth->Write();
      }
    }
    
    if(hasTruth){
      m_signal_mc_truth.hist->Write();
    }
    
    if(hasData){
      std::cout << " write out data" << m_selected_data.hist->GetName() << std::endl;
      m_selected_data.hist->Print();
      m_selected_data.hist->Write();
    }
  }
  
  //=======================================================================================
  // SYNC ALL HISTOGRAMS
  //=======================================================================================
  void SyncAllHists() {
    if(hasMC){
      m_selected_mc_reco.SyncCVHistos();
      m_selected_mc_truth.SyncCVHistos();
    }
    if(hasData){
      m_selected_data.SyncCVHistos();
    }
    if(hasTruth){
      m_signal_mc_truth.SyncCVHistos();
    }
    
  }
  inline void FillResponse(CVUniverse* univ, const double value, const double truth, const double weight=1.0){
    m_response->Fill(value, truth, univ->ShortName(), UniverseIndex(univ), weight);
  }
  
  // helper to return the actual numeric index corresponding to a universe  ie, allows map from name,index space to pure universe space.
  
  inline int UniverseIndex(CVUniverse* univ){
    return m_map[univ];
  }
  
  //
};

}  // namespace Ben

#endif  // VARIABLE_H
