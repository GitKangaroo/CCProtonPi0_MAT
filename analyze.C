bool DEBUG=0;

#include <iostream>
#include <string>
#include "PlotUtils/MnvH1D.h"
#include "PlotUtils/MnvH2D.h"
#include "MinervaUnfold/MnvResponse.h"
#include "MinervaUnfold/MnvUnfold.h"
#include "utils/RebinFlux.h"
#include "utils/SyncBands.h"
#include "TFile.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TError.h"
#include "TH2D.h"
#include "TF2.h"
#include <stdlib.h>
#include "PlotUtils/ChainWrapper.h"
#include "TCanvas.h"
#include "utils/NuConfig.h"
#include "PlotUtils/FluxReweighter.h"
//#include "utils/POTCounter.h"
//#define LOCAL  // this is for local  data samples


double targets = 3.23478E30; //nucleons - is this the same for MC or different

//bool neutrinoMode = false;
//
//std::string tag; //  = "NU" or "ANTINU"
//std::string datastring; // txt file with data list
//std::string mcstring;   // txt file with mc list

// ----- get all the usual code



#// header to make certain error bands don't lose their cV


#ifndef __CINT__
#include "include/plotting_pdf.h"
#endif

// pull in some code that defines cuts and setups up systematics


// Main  - this will do set up, loop over events and then make some plots


int main(const int argc, const char *argv[] ) {
  
  gROOT->ProcessLine("gErrorIgnoreLevel = kWarning;");
  //++++++++++++++++++=  Initialization +++++++++++++++++++++++++
  
  std::string pl = "1A";
  if (argc > 1){
    pl = std::string(argv[1]);
  }
  else{
    std::cout << " arguments are:\n analyze <config> [<prescale from config file>] " << std::endl;
    exit(0);
  }
  if (DEBUG) std::cout << " test NuConfig" << std::endl;
  NuConfig config;
  config.Read(pl+".json");
  config.Print();
  std::vector<std::string> variables = config.GetStringVector("AnalyzeVariables");
  bool LOCAL = config.GetBool("local");
  std::string playlist = config.GetString("playlist");
  std::string inputname = config.GetString("outRoot")+"_"+pl+".root";
  int m_fluxUniverses = config.GetInt("fluxUniverses");
  
  
  
  
  //========================================= Now do some analysis
  
  MnvH1D* h_flux_dewidthed = GetFlux(config);
  double flux = h_flux_dewidthed->Integral() ;
  // make containers for different analysis levels
  
  std::map<std::string, MnvH1D*> data;
  std::map<std::string, MnvH1D*> mc_sig;
  std::map<std::string, MnvH1D*> mc_bkg;
  std::map<std::string, MnvH1D*> mc; // sum of sig ang bkg
  std::map<std::string, MnvH1D*> selected;
  std::map<std::string, MnvH1D*> truth;
  std::map<std::string, MnvH2D*> migration;
  
  // derived
  std::map<std::string, MnvH1D*> signalFraction;
  std::map<std::string, MnvH1D*> unsmeared;
  std::map<std::string, MnvH1D*> bkgsub;
  std::map<std::string, MnvH1D*> efficiency;
  std::map<std::string, MnvH1D*> unsmearing; // just a temp when not doing unfolding
  std::map<std::string, MnvH1D*> effcorr;
  std::map<std::string, MnvH1D*> sigma;
  std::map<std::string, MnvH1D*> sigmaMC;
  
  std::map<std::string, MnvH1D*> h_flux_ebins;
  // tags for the histograms
  
  std::string sig = config.GetString("signal");
  std::string bkg = config.GetString("background");
  std::string dat = config.GetString("data");
  
  
  
  std::string truth_tag = sig+"_signal_mc_truth_";
  std::string selected_tag = sig+"_selected_mc_truth_";
  std::string mc_sig_tag = sig+"_selected_mc_reco_";
  std::string mc_bkg_tag = bkg+"_selected_mc_reco_";
  std::string data_tag = dat+"_selected_data_";
  std::string migration_tag = sig+"_response_selected_mc_reco_";
  std::string migration_end = "_migration";
  
  if(DEBUG) std::cout << " mc_bkg " << mc_bkg_tag << std::endl;
  
  
  // choose variables to analyze
  
  //std::vector<std::string> variables =  {"pT","enuQE","Log10Q2QE"};//, "thetaX","thetaY","theta","Log10Q2QE","recoil"};
  
  // get a normalization factor
  
  
  
  // now loop over histograms, unsmear, efficiency correct and normalize
  
  // unfolding setup
  MinervaUnfold::MnvUnfold unfold;
  MnvH2D* MigrationMatrix;
  double num_iter = 5;
  
  // input
  
  TFile* f = TFile::Open(inputname.c_str(),"READONLY");
  f->ls();
  
  TH1D* h_pot = (TH1D*)f->Get("POT_summary");
  double dataPOT = h_pot->GetBinContent(1);
  double mcPOTprescaled = h_pot->GetBinContent(3);
  double POTScale = dataPOT/mcPOTprescaled;
  
  double norm = 1./dataPOT/targets/flux;
  std::cout << "integrated luminosity is" << norm << std::endl;
  
  std::cout <<  "POT MC " <<  mcPOTprescaled  <<std::endl;
  std::cout <<  "POT DATA " << dataPOT << std::endl;
  
  
  
  for (auto name: variables){
    //  fdata->cd();
    std::cout << "tags" << data_tag+name << " " << mc_sig_tag+name << " " << mc_bkg_tag+name << std::endl;
    data[name] = (MnvH1D*)f->Get((data_tag + name).c_str());
    mc_sig[name] = (MnvH1D*)f->Get((mc_sig_tag + name).c_str());
    mc_bkg[name] = (MnvH1D*)f->Get((mc_bkg_tag + name).c_str());
    std::cout << selected_tag+name << " " << truth_tag+name << " " << migration_tag+name+migration_end << std::endl;
    selected[name] = (MnvH1D*)f->Get((selected_tag + name).c_str());
    truth[name] = (MnvH1D*)f->Get((truth_tag + name).c_str());
    migration[name] = (MnvH2D*)f->Get((migration_tag+name+migration_end).c_str());
    //std::cout << "migration " << migration_tag+name+migration_end << std::endl;
    //std::cout << " have filled histograms" << data[name] << " " <<  mc_sig[name] << " " << mc_bkg[name] << " " << truth[name] << " " << migration[name] << std::endl;
    
    
  }
  
  for (auto name: variables){
    data[name]->SetDirectory(0);
    mc_sig[name]->SetDirectory(0);
    mc_bkg[name]->SetDirectory(0);
    selected[name]->SetDirectory(0);
    truth[name]->SetDirectory(0);
    migration[name]->SetDirectory(0);
  }
  
  std::cout << " have filled histograms" << std::endl;
  
  // output
  std::string outroot = "analyze_"+inputname;
  TFile* o = TFile::Open(outroot.c_str(),"RECREATE");
  
  for (auto name: variables){
    mc_sig[name]->Scale(POTScale);
    mc_bkg[name]->Scale(POTScale);
    selected[name]->Scale(POTScale);
    truth[name]->Scale(POTScale);
    migration[name]->Scale(POTScale);
  }
  
  std::string pdfname = "analyze_"+inputname.replace(inputname.end()-5,inputname.end(),"");
  std::string pdffilename = pdfname + ".pdf";
  std::string pdfstart = pdfname + ".pdf(";
  std::string pdfend = pdfname + ".pdf)";
  TCanvas canvas(pdffilename.c_str());
  canvas.SetLeftMargin(0.15);
  canvas.SetBottomMargin(0.15);
  canvas.Print(pdfstart.c_str(),"pdf");
  
  bool binwid = true;  // flag you need if MnvPlotter does not do binwid correction.  If it does not, you need to set this to true
  
  o->cd();
  
  for (auto name: variables){  // only do this for a subset to save output time.
    int logscale = 0; // 0 for none, 1 for x, 2 for y, 3 for both
    
    if (name == "Q2QE") logscale = 1;
    if (DEBUG) std::cout << "logscale" << name << " " << logscale << std::endl;
    // look at the data and compare to mc
    std::cout << " starting on variable " << name << std::endl;
    data[name]->Write();
    
    mc_sig[name]->Write();
    mc_bkg[name]->Write();
    mc[name] = (MnvH1D*)mc_sig[name]->Clone((name+"_mc_tot").c_str());
    mc[name]->SetDirectory(0);
    mc[name]->Add(mc_bkg[name]);
    SyncBands(mc[name]);
    mc[name]->Write();
    
    signalFraction[name] = (MnvH1D*)mc_sig[name]->Clone((name+"_signalFraction").c_str());
    signalFraction[name]->SetDirectory(0);
    signalFraction[name]->Divide(mc_sig[name],mc[name],1.,1.,"B");
    SyncBands(signalFraction[name]);
    signalFraction[name]->Write();
    PlotCVAndError(canvas,data[name],mc[name], "DATA_vs_MC" ,true,logscale,binwid);
    PlotErrorSummary(canvas,mc[name],"MC errors" ,logscale );
    
    // here is where you would do background subtraction
    
    bkgsub[name] = data[name]->Clone((name+"_bkgsub").c_str());
    bkgsub[name]->SetDirectory(0);
    bkgsub[name]->AddMissingErrorBandsAndFillWithCV(*mc[name]);
    bkgsub[name]->Multiply(bkgsub[name],signalFraction[name]);
    SyncBands(bkgsub[name]);
    //bkgsub[name]->AddMissingErrorBandsAndFillWithCV(*mc[name]);
    bkgsub[name]->Write();
    
    PlotCVAndError(canvas,bkgsub[name],mc_sig[name], "BKGsub vs. MC signal" ,true,logscale,binwid);
    PlotErrorSummary(canvas,mc_sig[name],"MC signal  errors" ,logscale );
    
    // here is where you do unfolding - note that unsmeared gets filled by the unfolding but must NOT have error bands when fed in.
    // hence a clone of data and not bkgsub which might have them.  I found this out the hard way...
    
    unsmeared[name] = data[name]->Clone((name+"_unsmeared").c_str());
    unsmeared[name]->SetDirectory(0);
    
    // this section does a crude ratio unfolding on things that don't really need unfolding
    
    if (name == "recoil" || name == "recoilfull" || name == "multiplicity"){ // these guys don't have smearing really
      std::cout << " use ratio unfolding for " << name << std::endl;
      unsmearing[name] = bkgsub[name]->Clone((name+"_unsmearing").c_str());
      unsmearing[name]->SetDirectory(0);
      unsmearing[name]->Divide(selected[name],mc[name],1.,1.,"B");
      unsmeared[name]->Multiply(unsmeared[name],unsmearing[name],1.,1.);
      PlotCVAndError(canvas,mc[name],selected[name],"fraction_unsmearing" ,true,logscale,binwid);
      unsmearing[name]->Write();
    }
    else{ // real unfolding
      std::cout << " use Bayesian unfolding for " << name << std::endl;
      MigrationMatrix = migration[name];
      SyncBands(MigrationMatrix);
      
      // you can't the cv in as the unfolding code does not expect it.
      MigrationMatrix->PopVertErrorBand("cv");
      if (DEBUG)std::cout << " Migration matrix has size " << MigrationMatrix->GetErrorBandNames().size() << std::endl;
      if (DEBUG)std::cout << " Data has  size " << bkgsub[name]->GetErrorBandNames().size() << std::endl;
      
      // make an empty covariance matrix for the unfolding to give back to you
      TMatrixD covmatrix;
      bool data_unfolded = unfold.UnfoldHisto(unsmeared[name],covmatrix,MigrationMatrix,bkgsub[name],RooUnfold::kBayes, num_iter ,true,true);
      if (DEBUG)std::cout << " unfolding status " << data_unfolded << std::endl;
      
      // set the diagonal of the covmatrix to 0 as it is the same as the stat error on the central value.
      for (int i = 0; i < covmatrix.GetNrows(); i++){
        covmatrix[i][i] = 0.0;
      }
      unsmeared[name]->FillSysErrorMatrix("Unfolding",covmatrix);
      
    }
    
    SyncBands(unsmeared[name]);
    unsmeared[name]->Write();
    if (DEBUG) unsmeared[name]->Print();
    PlotCVAndError(canvas,bkgsub[name],unsmeared[name],"data before and after unsmearing", true,logscale,binwid);
    PlotCVAndError(canvas,unsmeared[name],selected[name], "unsmeared_data compared to selected MC" ,true,logscale,binwid);
    PlotErrorSummary(canvas,unsmeared[name],"unsmeared_data errors" ,logscale );
    
    // do an efficiency correction
    
    effcorr[name] = unsmeared[name]->Clone((name+"_effcorr").c_str());
    effcorr[name]->SetDirectory(0);
    
    // make the efficiency
    efficiency[name] = selected[name]->Clone((name+"_efficiency").c_str());
    efficiency[name]->AddMissingErrorBandsAndFillWithCV(*truth[name]);
    efficiency[name]->SetDirectory(0);
    efficiency[name]->Divide(efficiency[name],truth[name],1.,1.,"B");
    efficiency[name]->AddMissingErrorBandsAndFillWithCV(*unsmeared[name]);
    SyncBands(efficiency[name]);
    PlotCVAndError(canvas,selected[name],truth[name], "efficiency: selected and true" ,true,logscale,binwid);
    
    // apply the efficiency
    effcorr[name]->Divide(unsmeared[name],efficiency[name],1.,1.);
    SyncBands(effcorr[name]);
    effcorr[name]->Write();
    efficiency[name]->Write();
    if (DEBUG) effcorr[name]->Print();
    
    
    PlotCVAndError(canvas,effcorr[name],truth[name], "effcorr data vs truth" ,true,logscale,binwid);
    PlotErrorSummary(canvas,effcorr[name],"effcorr data/truth" ,logscale );
    
    
    // normalize to POT/targets/flux
    
    // if not an energy dependent cross section, normalize by total flux
    if(name != "enuQE" && name != "EnuHad"){
      binwid = true;
      sigma[name] = effcorr[name]->Clone((name+"_sigma").c_str());
      sigma[name]->SetDirectory(0);
      sigmaMC[name] = truth[name]->Clone((name+"_sigmaMC").c_str());
      sigmaMC[name]->SetDirectory(0);
      sigma[name]->Scale(norm);
      
      //sigma[name]->GetSysErrorMatrix("Unfolding").Print();
      sigmaMC[name]->Scale(norm);
      
    }
    else{  // this is energy dependent part which uses the full flux
      h_flux_ebins[name] = GetFluxEbins(h_flux_dewidthed,truth[name]);
      
      binwid = false;
      double Enorm = 1./targets/dataPOT;  // norm factor without flux
      sigmaMC[name] = truth[name]->Clone((name+"_sigmaMC").c_str());
      sigmaMC[name]->SetDirectory(0);
      sigma[name] = (MnvH1D*)effcorr[name]->Clone((name+"_sigma").c_str());
      sigma[name]->SetDirectory(0);
      if (DEBUG) std::cout << " about to divide sigma by flux " << Enorm <<  std::endl;
      //sigma[name]->Print("ALL");
      if (DEBUG)std::cout << " try to make error bands consistent" << std::endl;
      //sigma[name]->GetSysErrorMatrix("Unfolding").Print();
      h_flux_ebins[name]->AddMissingErrorBandsAndFillWithCV(*sigma[name]);
      sigma[name]->AddMissingErrorBandsAndFillWithCV(*h_flux_ebins[name]);
      if (DEBUG)std::cout << " after add missing error bands " << std::endl;
      //sigma[name]->GetSysErrorMatrix("Unfolding").Print();
      h_flux_ebins[name]->AddMissingErrorBandsAndFillWithCV(*sigmaMC[name]);
      sigmaMC[name]->AddMissingErrorBandsAndFillWithCV(*h_flux_ebins[name]);
      //h_flux_ebins->Print("ALL");
      sigma[name]->Scale(Enorm);
      if (DEBUG)std::cout << " after scale " << std::endl;
      if (DEBUG) sigma[name]->GetSysErrorMatrix("Unfolding").Print();
      sigmaMC[name]->Scale(Enorm);
      sigma[name]->Divide(sigma[name],h_flux_ebins[name]);
      if (DEBUG)std::cout << "just did the divide" << std::endl;
      sigmaMC[name]->Divide(sigmaMC[name],h_flux_ebins[name],1.0,1.0);
      
    }
    SyncBands(sigma[name]);
    SyncBands(sigmaMC[name]);
    sigma[name]->Write();
    sigma[name]->Print();
    if(DEBUG) sigma[name]->Print("ALL");
    sigmaMC[name]->Write();
    sigmaMC[name]->Print();
    double id = integrator(sigma[name],binwid);
    double im = integrator(sigmaMC[name],binwid);
    std::cout << "sigma " << name << " " << id << " " << sigma[name]->Integral() << " " << im << " " << sigmaMC[name]->Integral() << std::endl;
    PlotCVAndError(canvas,sigma[name],sigmaMC[name], "sigma" ,true,logscale,binwid);
    PlotErrorSummary(canvas,sigma[name],"sigma" ,logscale );
    
    std::cout << " end of loop over types " << std::endl;
   
  }
  
  std::cout << " end of loops" << std::endl;
  //
  canvas.Print(pdfend.c_str(),"pdf");
  
  o->Close();
  exit(1);
  //  for (auto h:histsets){
  //    h->DeleteResponse();
  //  }
  
}
