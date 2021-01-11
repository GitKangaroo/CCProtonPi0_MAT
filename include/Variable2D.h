#ifndef VARIABLE2D_H
#define VARIABLE2D_H

#include "CVUniverse.h"
#include "PlotUtils/Hist2DWrapper.h"
#include "PlotUtils/MnvH2D.h"

#ifndef __CINT__
#include "PlotUtils/Variable2DBase.h"
#endif

namespace CCProtonPi0 {
class Variable2D : public PlotUtils::Variable2DBase<CVUniverse> {
 private:
  //=======================================================================================
  // TYPEDEFS CONVENIENCE
  //=======================================================================================
  typedef PlotUtils::Hist2DWrapper<CVUniverse> HW2D;
  typedef PlotUtils::MnvH2D MH2D;

 public:
  //=======================================================================================
  // CTOR
  //=======================================================================================
  template <class... ARGS>
  Variable2D(ARGS... args) : Variable2DBase<CVUniverse>(args...) {}

  //=======================================================================================
  // DECLARE NEW HISTOGRAMS
  //=======================================================================================
  // HISTWRAPPER
  HW2D m_selected_mc_reco;

  //// HISTFOLIO
  // PlotUtils::HistFolio<MH2D> m_selected_mc_sb;

  //=======================================================================================
  // INITIALIZE ALL HISTOGRAMS
  //=======================================================================================
  template <typename T>
  void InitializeAllHistograms(T univs) {
    const bool clear_bands = true;  // we want empty histograms
    const char* name = GetName().c_str();

    // HISTWRAPPER
    // selected mc reco histwrapper
    MH2D* dummy_selected_mc_reco =
        new MH2D(Form("selected_mc_reco_%s", name), name, GetNBinsX(),
                 GetBinVecX().data(), GetNBinsY(), GetBinVecY().data());
    m_selected_mc_reco = HW2D(dummy_selected_mc_reco, univs, clear_bands);

    delete dummy_selected_mc_reco;
  }

  //=======================================================================================
  // WRITE ALL HISTOGRAMS
  //=======================================================================================
  void WriteAllHistogramsToFile(TFile& f) const {
    f.cd();

    // selected mc reco
    m_selected_mc_reco.hist->Write();
  }
};
}  // namespace Ben

#endif  // VARIABLE2D_H
