#ifndef plotting_functions_H
#define plotting_functions_H
#include <iostream>
#include <unordered_map>

#include "PlotUtils/HistogramUtils.h"
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvPlotter.h"
#include "PlotUtils/MnvVertErrorBand.h"
#include "Variable.h"
#include "Variable2D.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TStyle.h"

namespace plotting {
const double xmin = 0.;
const double xmax = 20.e3;
const int nbins = 30;

const bool do_fractional_uncertainty = true;
const bool do_cov_area_norm = false;
const bool include_stat_error = false;
}  // namespace plotting

void PlotTH1(TH1D* h, std::string filename) {
  TCanvas cF ("c4","c4"); 
  h->Draw();
  cF.Print(Form("%s.png", filename.c_str()));
}

const std::string do_fractional_uncertainty_str =
    plotting::do_fractional_uncertainty ? std::string("Frac")
                                        : std::string("Abs");
const std::string do_cov_area_norm_str =
    plotting::do_cov_area_norm ? std::string("CovAreaNorm") : std::string("");

void PlotErrorSummary(PlotUtils::MnvH1D* hist, std::string label);
void PlotVertBand(std::string band, std::string method_str,
                  PlotUtils::MnvH1D* hist);
void PlotVertBandCV(std::string band, std::string method_str,
                  PlotUtils::MnvH1D* h);
void PlotLatBand(std::string band, std::string method_str,
                 PlotUtils::MnvH1D* hist);
void PlotVertUniverse(std::string band, unsigned int universe,
                      std::string method_str, PlotUtils::MnvH1D* hist);
void PlotLatUniverse(std::string band, unsigned int universe,
                     std::string method_str, PlotUtils::MnvH1D* hist);
void PlotCVAndError(PlotUtils::MnvH1D* hist, std::string label);
void PlotTotalError(PlotUtils::MnvH1D* hist, std::string method_str);

void PlotTotalError(PlotUtils::MnvH1D* hist, std::string method_str) {
  TH1D* hTotalErr = (TH1D*)hist
                        ->GetTotalError(plotting::include_stat_error,
                                        plotting::do_fractional_uncertainty,
                                        plotting::do_cov_area_norm)
                        .Clone(Form("h_total_err_errSum_%d", __LINE__));
  TCanvas cF("c4", "c4");
  hTotalErr->SetTitle(
      Form("Total Uncertainty (%s); E_{#nu} (MeV)", method_str.c_str()));
  hTotalErr->Draw();
  cF.Print(Form("png/Enu_TotalUncertainty_%s_%s_%s.png",
                do_fractional_uncertainty_str.c_str(),
                do_cov_area_norm_str.c_str(), method_str.c_str()));
}

void PlotErrorSummary(PlotUtils::MnvH1D* h, std::string label) {
  PlotUtils::MnvH1D* hist = (PlotUtils::MnvH1D*)h->Clone("hist");
  PlotUtils::MnvPlotter mnvPlotter(PlotUtils::kCCNuPionIncStyle);
  TCanvas cE("c1", "c1");
  hist->GetXaxis()->SetTitle(label.c_str());


    mnvPlotter.error_summary_group_map.clear();
    mnvPlotter.error_summary_group_map["Genie_FSI"].push_back("GENIE_FrAbs_N");
    mnvPlotter.error_summary_group_map["Genie_FSI"].push_back("GENIE_FrAbs_pi");
    mnvPlotter.error_summary_group_map["Genie_FSI"].push_back("GENIE_FrCEx_N");
    mnvPlotter.error_summary_group_map["Genie_FSI"].push_back("GENIE_FrCEx_pi");
    mnvPlotter.error_summary_group_map["Genie_FSI"].push_back("GENIE_FrElas_N");
    mnvPlotter.error_summary_group_map["Genie_FSI"].push_back("GENIE_FrElas_pi");
    mnvPlotter.error_summary_group_map["Genie_FSI"].push_back("GENIE_FrInel_N");
    //mnvPlotter.error_summary_group_map["Genie_FSI"].push_back("GENIE_FrInel_pi");
    mnvPlotter.error_summary_group_map["Genie_FSI"].push_back("GENIE_FrPiProd_N");
    mnvPlotter.error_summary_group_map["Genie_FSI"].push_back(
        "GENIE_FrPiProd_pi");
    mnvPlotter.error_summary_group_map["Genie_FSI"].push_back("GENIE_MFP_N");
    mnvPlotter.error_summary_group_map["Genie_FSI"].push_back("GENIE_MFP_pi");
    mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_AGKYxF1pi");
    mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_AhtBY");
    mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_BhtBY");
    mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_CCQEPauliSupViaKF");
    mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_CV1uBY");
    mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_CV2uBY");
    mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_EtaNCEL");
    //mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
    //    "GENIE_MaCCQE");
    //mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
    //    "GENIE_MaCCQEshape");
    mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_MaNCEL");
    mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_MaRES");
    mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_MvRES");
    //mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
    //    "GENIE_NormCCQE");
    //mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
    //    "GENIE_NormCCRES");
    mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_NormDISCC");
    mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_NormNCRES");
    mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_RDecBR1gamma");
    mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_Rvn1pi");
    mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_Rvn2pi");
    mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_Rvn3pi");
    mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_Rvp1pi");
    mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_Rvp2pi");
    mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_Theta_Delta2Npi");
    mnvPlotter.error_summary_group_map["Genie_InteractionModel"].push_back(
        "GENIE_VecFFCCQEshape");

    mnvPlotter.error_summary_group_map["RPA"].push_back("RPA_HighQ2");
    mnvPlotter.error_summary_group_map["RPA"].push_back("RPA_LowQ2");
  
    mnvPlotter.error_summary_group_map["Muon_Energy"].push_back("Muon_Energy_MINOS");
    mnvPlotter.error_summary_group_map["Muon_Energy"].push_back("Muon_Energy_MINERvA");
    mnvPlotter.error_summary_group_map["Muon_Energy"].push_back("MinosEfficiency");
    mnvPlotter.error_summary_group_map["Muon_Energy"].push_back("Muon_Energy_Resolution");
    mnvPlotter.error_summary_group_map["Unfolding"].push_back("Unfolding");


  mnvPlotter.DrawErrorSummary(hist, "TR", plotting::include_stat_error, true,
                              0.0, plotting::do_cov_area_norm, "",
                              plotting::do_fractional_uncertainty);
  std::string plotname =
      Form("png/ErrorSummary_%s_%s_%s", do_fractional_uncertainty_str.c_str(),
           do_cov_area_norm_str.c_str(), label.c_str());
  mnvPlotter.MultiPrint(&cE, plotname, "png");
}

void PlotVertBandCV(std::string band, std::string method_str,
                  PlotUtils::MnvH1D* h) {
  PlotUtils::MnvH1D* hist = (PlotUtils::MnvH1D*)h->Clone("hist");
  TH1* h1 = (TH1*)hist->GetVertErrorBand(band.c_str())->Clone(Form("%s_%s", band.c_str(), method_str.c_str()));
  TCanvas cF("c4", "c4");
  h1->SetTitle(Form("%s Uncertainty; %s", band.c_str(), method_str.c_str()));
  h1->SetLineStyle(1);
  h1->Draw("hist");
  cF.Print(Form("png/%s_band_CV_%s.png", band.c_str(), method_str.c_str()));
}

void PlotVertBand(std::string band, std::string method_str,
                  PlotUtils::MnvH1D* h) {
  PlotUtils::MnvH1D* hist = (PlotUtils::MnvH1D*)h->Clone("hist");
  TH1* h1 = (TH1*)hist->GetVertErrorBand(band.c_str())
                ->GetErrorBand(plotting::do_fractional_uncertainty,
                               plotting::do_cov_area_norm)
                .Clone(Form("%s_%s", band.c_str(), method_str.c_str()));
  TCanvas cF("c4", "c4");
  h1->SetTitle(Form("%s Uncertainty; %s", band.c_str(), method_str.c_str()));
  h1->Draw("h");
  cF.Print(Form("png/%s_band_%s.png", band.c_str(), method_str.c_str()));
}

void PlotVertUniverse(std::string band, unsigned int universe,
                      std::string method_str, PlotUtils::MnvH1D* h) {
  PlotUtils::MnvH1D* hist = (PlotUtils::MnvH1D*)h->Clone("hist");
  TH1D* h1 = hist->GetVertErrorBand(band.c_str())->GetHist(universe);

  TCanvas cF("c1", "c1");
  h1->SetLineColor(kBlack);
  h1->SetLineStyle(1);
  h1->Draw("hist");
  cF.Print(Form("png/%s_band_universe%i_%s.png", band.c_str(), universe,
                method_str.c_str()));
}

void PlotLatUniverse(std::string band, unsigned int universe,
                     std::string method_str, PlotUtils::MnvH1D* hist) {
  TH1D* h1 = hist->GetLatErrorBand(band.c_str())->GetHist(universe);
  TCanvas cF("c1", "c1");
  h1->SetLineColor(kBlack);
  h1->SetLineStyle(1);
  h1->Draw("hist");
  cF.Print(Form("png/Enu_%s_band_universe%i_%s.png", band.c_str(), universe + 1,
                method_str.c_str()));
}

void PlotMCAndError(PlotUtils::MnvH1D* h, std::string label) {
  PlotUtils::MnvH1D* hist = (PlotUtils::MnvH1D*)h->Clone("hist");
  PlotUtils::MnvPlotter mnvPlotter(PlotUtils::kCCNuPionIncStyle);
  TCanvas cE("c1", "c1");
  PlotUtils::MnvH1D* datahist = new PlotUtils::MnvH1D(
      "adsf", "", plotting::nbins, plotting::xmin, plotting::xmax);
  bool statPlusSys = true;
  int mcScale = 1.;
  bool useHistTitles = false;
  const PlotUtils::MnvH1D* bkgdHist = NULL;
  const PlotUtils::MnvH1D* dataBkgdHist = NULL;
  mnvPlotter.DrawDataMCWithErrorBand(datahist, hist, mcScale, "TL",
                                     useHistTitles, NULL, NULL, false,
                                     statPlusSys);
  std::string plotname = Form("MCCV_w_err_%s", label.c_str());
  mnvPlotter.MultiPrint(&cE, plotname, "png");
  
  
  
  delete datahist;
}

void PlotDataMCAndError(PlotUtils::MnvH1D* d, PlotUtils::MnvH1D* m, double data_mc_pot_ratio, std::string label) {
  PlotUtils::MnvH1D* data = (PlotUtils::MnvH1D*)d->Clone("data");
  PlotUtils::MnvH1D* mc = (PlotUtils::MnvH1D*)m->Clone("mc");
  PlotUtils::MnvPlotter mnvPlotter(PlotUtils::kCCNuPionIncStyle);
  TCanvas cE("c1", "c1");
  bool do_bin_width_norm = true;
  bool statPlusSys = true;
  bool useHistTitles = false;
  const PlotUtils::MnvH1D* bkgdHist = NULL;
  const PlotUtils::MnvH1D* dataBkgdHist = NULL;

  if ( do_bin_width_norm ) {
    if (data) data->Scale(1., "width");
    mc->Scale(1., "width");
  }

  mnvPlotter.DrawDataMCWithErrorBand(data, mc, data_mc_pot_ratio, "TL",
                                     useHistTitles, NULL, NULL, false,
                                     statPlusSys);
  std::string plotname = Form("png/EventSelection_DataMC_%s", label.c_str());
  mnvPlotter.MultiPrint(&cE, plotname, "png");
  
  

   MnvH1D* d2 = (PlotUtils::MnvH1D*)d->Clone();
  MnvH1D* m2 = (PlotUtils::MnvH1D*)m->Clone();
 
   mc->SetDirectory(0);
   mc->ClearAllErrorBands();
   mc->AddMissingErrorBandsAndFillWithCV(*m);
   
   d2->Divide(d,mc,1.,1.);
  
   d2->GetYaxis()->SetTitle("Data/MC");
   m2->Divide(m,mc,1.,1.);
   // no idea why I should have to do this
   //m->Print("ALL");
   //datahist->Draw();
   //cov_area = false;
   //mnvPlotter.SetBinWidthNorm(false);
   mnvPlotter.DrawDataMCRatio(d2, m2,data_mc_pot_ratio,true,true);// "TL", useHistTitles, NULL, NULL,cov_area, statPlusSys);
   
   //%t->Draw();
   std::string plotname2 = Form("png/ratio_%s_CV_w_err_%s",d->GetName() ,label.c_str());
    mnvPlotter.MultiPrint(&cE, plotname2, "png");
   //cE.Print(cE.GetName(),plotname2.c_str());

}

void PlotStacked(const MnvH1D* data,const TObjArray& array_mc, const double pot_scale, std::string outfile_tag = "",
                 std::string plot_title = "", double ymax = -1) {
  // Never don't clone when plotting
  TObjArray array = *(TObjArray*)array_mc.Clone("mc");

  PlotUtils::MnvPlotter mnv_plotter;
  if (ymax > 0) mnv_plotter.axis_maximum = ymax;
  std::string outfile_name = Form("StackedBreakdown_%s", outfile_tag.c_str());
  TCanvas cE("c1", "c1");

  std::string x_label = ((PlotUtils::MnvH1D*)array.At(0))->GetXaxis()->GetTitle();
  std::string y_label = ((PlotUtils::MnvH1D*)array.At(0))->GetYaxis()->GetTitle();

  int nbins = ((PlotUtils::MnvH1D*)array.At(0))->GetNbinsX();
  double xmin = ((PlotUtils::MnvH1D*)array.At(0))->GetXaxis()->GetXmin();
  double xmax = ((PlotUtils::MnvH1D*)array.At(0))->GetXaxis()->GetXmax();

  //double pot_scale = 1.;
  //PlotUtils::MnvH1D* data = new PlotUtils::MnvH1D(
  //    "dummy", "dummy", nbins, xmin, xmax);

  // By default, this function uses the data x and y labels, instead of the mc
  // x and y labels. Given that this is much more commonly used with MC, this
  // seems like a bad choice, but who am I to break backwards compatibility?
  // Anyway, it means that we need to provide the x and y labels manually.
  mnv_plotter.DrawDataStackedMC(data, &array, pot_scale, "TR", "Data", -1, -1,
                                1001, x_label.c_str(), y_label.c_str());

  mnv_plotter.WritePreliminary("TL");
  mnv_plotter.AddHistoTitle(plot_title.c_str());
  mnv_plotter.MultiPrint(&cE, outfile_name, "png");

}

void Plot2D(PlotUtils::MnvH2D* hist, std::string name, std::string label_xaxis,
            std::string label_yaxis) {

    TCanvas c ("c1","c1"); 
    bool draw_as_matrix = false;
    hist->Draw("COLZ");
    hist->GetXaxis()->SetTitle(label_xaxis.c_str());
    hist->GetYaxis()->SetTitle(label_yaxis.c_str());
    c.Update();
    c.Print(Form("%s_2D.png", name.c_str()));
}

namespace
{
  std::string SafeROOTName(std::string copy)
  {
    //Replace decimal points with underscores
    if(copy.find_first_of(".") != std::string::npos) copy.replace(copy.find_first_of("."), 1, "_");

    //Replace spaces with underscores
    size_t spacePos = std::string::npos;
    while((spacePos = copy.find_first_of(" ")) != std::string::npos)
    {
      copy.replace(spacePos, 1, "_");
    }

    //Replace characters reserved in c++ with underscores
    size_t specialPos = std::string::npos;
    while((specialPos = copy.find_first_of("+-/.*&<>,{}()^|\\")) != std::string::npos)
    {
      copy.replace(specialPos, 1, "_");
    }

    return copy;
  }
}

void DrawStack(const std::string& cutName, CCProtonPi0::Variable* var, std::unordered_map<std::string, double>& histNameToMax, const std::string& histName)
{
  TCanvas canvas(cutName.c_str());
  THStack stack((histName + "_before_" + cutName).c_str(), (histName + "_before_" + cutName).c_str());
  auto histMap = var->m_selected_mc_by_channel.GetHistMap();
  assert(histMap.size() > 0 && "Drawing a HistFolio with 0 entries will result in a crash here.");
  stack.SetHistogram(histMap.begin()->second->Clone());
  stack.SetMaximum(histNameToMax[histName]);
  for(auto& pair: histMap) stack.Add(pair.second);
  stack.Draw("HIST");

  //Draw a legend that doesn't overlap the peak if possible
  const double maxBinCenter =  stack.GetHistogram()->GetBinCenter(((TH1*)(stack.GetStack()->Last()))->GetMaximumBin()),
               xMax = stack.GetHistogram()->GetXaxis()->GetXmax();
  const double xRange = 0.25;
  double xStart = 0.15;
  if(maxBinCenter < xMax/2.)
  {
    xStart = 0.85 - xRange;
  }

  canvas.BuildLegend(xStart, 0.45, xStart + xRange, 0.85); //x1, y1, x2, y2 in Normalized Device Coordinates ([0, 1])
  canvas.Print((::SafeROOTName(std::string("GENIE_breakdown_before_") + cutName + "_cut_in_" + histName) + ".png").c_str());
}

#endif
