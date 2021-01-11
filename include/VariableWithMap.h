#ifndef VARIABLEWithMap_H
#define VARIABLEWithMap_H

#include "CVUniverse.h"
#include "PlotUtils/HistFolio.h"
#include "HistWrapperMap.h"
#include "MinervaUnfold/MnvResponse.h"


#ifndef __CINT__  // CINT doesn't know about std::function
#include "PlotUtils/VariableBase.h"
#endif  // __CINT__

namespace CCProtonPi0 {

class VariableWithMap : public PlotUtils::VariableBase<CVUniverse> {
private:
  typedef PlotUtils::HistWrapperMap<CVUniverse> HM;
  
  typedef PlotUtils::MnvH1D MH1D;
  typedef PlotUtils::HistFolio<PlotUtils::MnvH1D> FOLIO;
  
  
public:
  //=======================================================================================
  // CTOR
  //=======================================================================================
  template <class... ARGS>
  VariableWithMap(ARGS... args) : PlotUtils::VariableBase<CVUniverse>(args...) {}
  typedef std::map<std::string, std::vector<CVUniverse*>> UniverseMap;
  //=======================================================================================
  // DECLARE NEW HISTOGRAMS
  //=======================================================================================
  // Histwrappers -- selected mc, selected data
  
  HM m_selected_mc_reco;
  HM m_selected_mc_truth;
  HM m_signal_mc_truth;
  HM m_selected_data;
  UniverseMap m_universes;
  std::string m_units;
  std::map<const std::string, bool> hasData;
  std::map<const std::string, bool> hasMC;
  std::map<const std::string, bool> hasTruth;
  std::map<const std::string, bool> hasResponse;
  std::vector<std::string> m_tags;
  //RESPONSE* m_response;
  // helpers for response
  
  // std::map< CVUniverse* , int> m_map; // map to get name and index from Universe;
  // Histofolio to categorize MC by interaction channel
  //FOLIO m_selected_mc_by_channel;
  //FOLIO m_selected_mc_truth_by_channel;
  inline virtual std::string GetUnits(){return m_units;};
  
  inline virtual void SetUnits(std::string units){m_units=units;};
  
  //=======================================================================================
  // INITIALIZE ALL HISTOGRAMS
  //=======================================================================================
  
  
  void AddTags(const std::vector<std::string> tags){m_tags = tags;};
  
  template <typename T>
  
  void InitializeMCHistograms(T univs, const std::vector<std::string> tags) {
    for (auto tag:tags){
      hasMC[tag] = true;
    }
    std::vector<double> bins = GetBinVec();
    const char* name = GetName().c_str();
    
    m_selected_mc_reco = HM(Form("selected_mc_reco_%s", name), (GetName()+";"+m_xaxis_label).c_str(), GetNBins(), bins, univs, tags);
    //Interaction types from https://nusoft.fnal.gov/minerva/minervadat/software_doxygen/HEAD/MINERVA/classMinerva_1_1GenMinInteraction.html
    //  std::vector<PlotUtils::NamedCategory<int>> genieCategories = {{1, "QE"}, {2, "RES"}, {3, "DIS"}, {4, "Coherent Pi"}, {8, "MEC"}};
    //  m_selected_mc_by_channel = FOLIO(genieCategories, std::string("selected_mc_in_") + name + "by_channel", GetNBins(), bins.data());
    m_selected_mc_truth = HM(Form("selected_mc_truth_%s", name), (GetName()+";"+m_xaxis_label).c_str(), GetNBins(), bins, univs, tags);
  }
  
  
  
  
  
  template <typename T>
  void InitializeTruthHistograms(T univs, const std::vector<std::string> tags) {
    for (auto const tag:tags){
      hasTruth[tag] = true;
    }
    std::vector<double> bins = GetBinVec();
    const char* name = GetName().c_str();
    
    m_signal_mc_truth = HM(Form("signal_mc_truth_%s", name), (GetName()+";"+m_xaxis_label).c_str(), GetNBins(), bins, univs, tags);
    
  }
  
  //  void AddTruthWrapper(std::string tag){
  //     m_signal_mc_truth.AddResponse(tag);
  //  }
  
  //========== Add Response =================
  
  
  void AddMCResponse(const std::vector<std::string>  tags) {
    for (auto tag:tags){
      assert(hasMC[tag]);
    }
    m_selected_mc_reco.AddResponse(tags);
    for (auto tag:tags){
      hasResponse[tag] = true;
    }
    
  }
  
  
  
  
  template <typename T>
  void InitializeDataHistograms(T univs, const std::vector<std::string> tags) {
    for (auto tag:tags){
      hasData[tag] = true;
    }
    std::vector<double> bins = GetBinVec();
    const char* name = GetName().c_str();
    m_selected_data = HM(Form("selected_data_%s", name), (GetName()+";"+m_xaxis_label).c_str(), GetNBins(), bins, univs,tags);
  }
  
  
  
  //=======================================================================================
  // WRITE ALL HISTOGRAMS
  //=======================================================================================
  void WriteAllHistogramsToFile(TFile& f)  {
    std::cout << "should only be called once " << std::endl;
    f.cd();
    
    
    // selected mc reco
    
    for (auto tag:m_tags){
      std::cout << " write out flags " << hasMC[tag] << hasTruth[tag] << hasData[tag] <<  std::endl;
      if(hasMC[tag]) {
        m_selected_mc_reco.Write(tag);
        std::cout << " write out mc histogram " << m_selected_mc_reco.GetHist(tag)->GetName() << std::endl;
        // m_selected_mc_by_channel.WriteToFile(f);
        m_selected_mc_truth.Write(tag);
      }
      
      if(hasTruth[tag]){
        std::cout << " write out truth histogram " << m_selected_mc_truth.GetHist(tag)->GetName() << std::endl;
        m_signal_mc_truth.Write(tag);
      }
      
      if(hasData[tag]){
        std::cout << " write out data histogram " << m_selected_data.GetHist(tag)->GetName() << std::endl;
        // m_selected_data.hist->Print();
        m_selected_data.Write(tag);
      }
    }
  }
  
  //=======================================================================================
  // SYNC ALL HISTOGRAMS
  //=======================================================================================
  void SyncAllHists() {
    for (auto tag:m_tags){
      if(hasMC[tag]){
        m_selected_mc_reco.SyncCVHistos();
        m_selected_mc_truth.SyncCVHistos();
      }
      if(hasData[tag]){
        m_selected_data.SyncCVHistos();
      }
      if(hasTruth[tag]){
        m_signal_mc_truth.SyncCVHistos();
      }
    }
  }
  
  
  
  inline void FillResponse(std::string tag, CVUniverse* univ, const double value, const double truth, const double weight=1.0){
    m_selected_mc_reco.FillResponse(tag, univ, value, truth, weight);
  }
  
  // helper to return the actual numeric index corresponding to a universe  ie, allows map from name,index space to pure universe space.
  
  //  inline int UniverseIndex(CVUniverse* univ){
  //    return m_map[univ];
  //  }
  //
  //  //
};

}  // namespace Ben

#endif  // VARIABLE_H
