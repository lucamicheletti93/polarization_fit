#if !defined(__CINT__) || defined(__MAKECINT__)
#include <stdio.h>

#include <TMinuit.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TF2.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TPad.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TFitResult.h>
#include <TMatrixDSym.h>
#include <TPaveText.h>
#include <TGaxis.h>
#endif

void helicity_polarization_fit_2D(){

  gSystem -> CompileMacro("../signal_extraction/settings.h");
  gROOT -> ProcessLine(".x ../signal_extraction/binning.C");

  gStyle -> SetOptStat(0);
  gStyle -> SetOptFit(1);
  TGaxis::SetMaxDigits(2);

  //============================================================================
  // N_Jpsi HISTOGRAMS
  //============================================================================

  //----------------------------------------------------------------------------
  //RANGE 0 < phi < PI ; -1 < cost < 1 [BC = Bent&Cut]
  //----------------------------------------------------------------------------

  int matrix_N_Jpsi_HE[N_cost_bins][N_phi_bins];
  int matrix_stat_N_Jpsi_HE[N_cost_bins][N_phi_bins];

  TFile *N_Jpsi_file = new TFile("/home/luca/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/SIGNAL_EXTRACTION/HISTOS_FOR_SIGNAL_EXTRACTION/GIT_OUTPUT/N_Jpsi.root","READ");
  //TFile *N_Jpsi_file = new TFile("/home/luca/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/POLARIZED_DISTRIBUTIONS/sum_variable_Jpsi_polarization_LowStat.root");
  //TFile *N_Jpsi_file = new TFile("/home/luca/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/POLARIZED_DISTRIBUTIONS/variable_Jpsi_polarization_sample_LowStat_reduced.root");

  TH2D *hist_N_Jpsi_2pt6_HE = new TH2D("hist_N_Jpsi_2pt6_HE","hist_N_Jpsi_2pt6_HE",N_cost_bins,value_cost,N_phi_bins,value_phi);
  hist_N_Jpsi_2pt6_HE -> GetXaxis() -> SetLabelSize(0);
  hist_N_Jpsi_2pt6_HE -> GetYaxis() -> SetLabelSize(0);

  TH2D *hist_N_Jpsi_2pt6_HE = (TH2D*) N_Jpsi_file -> Get("hist_N_Jpsi_HE");
  //TH2D *hist_N_Jpsi_2pt6_HE = (TH2D*) N_Jpsi_file -> Get("hist_rec_polarization_Rebin14");
  hist_N_Jpsi_2pt6_HE -> Sumw2();
  //hist_N_Jpsi_2pt6_HE -> Draw();

  //int counter = 0;
  //TH1D *hist_rel_err = new TH1D("hist_rel_err","hist_rel_err",180,0,180);
  //for(int i = 0;i < 18;i++){
    //for(int j = 0;j < 10;j++){
      //hist_rel_err -> SetBinContent(counter+1,hist_N_Jpsi_2pt6_HE -> GetBinError(i+1,j+1)/hist_N_Jpsi_2pt6_HE -> GetBinContent(i+1,j+1));
      //counter++;
    //}
  //}

  //hist_rel_err -> Draw();

  /*TH2D *hist_N_Jpsi_area_2pt6_HE = new TH2D("hist_N_Jpsi_area_2pt6_HE","hist_N_Jpsi_area_2pt6_HE",N_cost_bins,value_cost,N_phi_bins,value_phi);

  for(int i = 0;i< N_cost_bins;i++){
    for(int j = 0;j < N_phi_bins;j++){
      hist_N_Jpsi_area_2pt6_HE -> SetBinContent(i+1,j+1,(hist_N_Jpsi_2pt6_HE -> GetBinContent(i+1,j+1))/bin_area[i][j]);
      hist_N_Jpsi_area_2pt6_HE -> SetBinError(i+1,j+1,(hist_N_Jpsi_2pt6_HE -> GetBinError(i+1,j+1))/bin_area[i][j]);
      //printf("%i +- %i \n",(hist_N_Jpsi_accxeff_corrected_2pt6_HE -> GetBinContent(i+1,j+1))/bin_area[i][j],(hist_N_Jpsi_accxeff_corrected_2pt6_HE -> GetBinError(i+1,j+1))/bin_area[i][j]);
    }
  }

  TCanvas *ccc = new TCanvas("ccc","ccc",20,20,600,600);
  hist_N_Jpsi_area_2pt6_HE -> Draw("COLZtext");*/

  printf("%f +- %f \n",hist_N_Jpsi_2pt6_HE -> GetBinContent(10,10),hist_N_Jpsi_2pt6_HE -> GetBinError(10,10));

  //============================================================================
  // ACCXEFF HISTOGRAMS
  //============================================================================

  TFile *accxeff_file = new TFile("/home/luca/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/ACCXEFF/HISTOS_FOR_ACCXEFF/GIT_OUTPUT/accxeff.root","READ");
  TH2D *hist_accxeff_2pt6_HE = (TH2D*) accxeff_file -> Get("hist_accxeff_HE_2pt6_rebin");

  TCanvas *cccc = new TCanvas("cccc","cccc",10,10,600,600);
  hist_accxeff_2pt6_HE -> Draw("COLZtext");

  printf("%f +- %f \n",hist_accxeff_2pt6_HE -> GetBinContent(10,10),hist_accxeff_2pt6_HE -> GetBinError(10,10));

  //============================================================================
  // N_Jpsi CORRECTED FOR ACCXEFF DISTRIBUTION
  //============================================================================

  TH2D *hist_N_Jpsi_accxeff_corrected_2pt6_HE = new TH2D("hist_N_Jpsi_accxeff_corrected_2pt6_HE","hist_N_Jpsi_accxeff_corrected_2pt6_HE",N_cost_bins,value_cost,N_phi_bins,value_phi);
  hist_N_Jpsi_accxeff_corrected_2pt6_HE -> Sumw2();
  hist_N_Jpsi_accxeff_corrected_2pt6_HE -> Divide(hist_N_Jpsi_2pt6_HE,hist_accxeff_2pt6_HE,1,1);

  printf("%f +- %f \n",hist_N_Jpsi_accxeff_corrected_2pt6_HE -> GetBinContent(10,10),hist_N_Jpsi_accxeff_corrected_2pt6_HE -> GetBinError(10,10));

  //============================================================================
  // N_Jpsi CORRECTED FOR ACCXEFF AND BIN-AREA DISTRIBUTION
  //============================================================================

  //----------------------------------------------------------------------------
  // RANGE 0 < phi < PI ; -1 < cost < 1
  //----------------------------------------------------------------------------
  TH2D *hist_N_Jpsi_accxeff_area_corrected_2pt6_HE = new TH2D("hist_N_Jpsi_accxeff_area_corrected_2pt6_HE","hist_N_Jpsi_accxeff_area_corrected_2pt6_HE",N_cost_bins,value_cost,N_phi_bins,value_phi);

  for(int i = 0;i< N_cost_bins;i++){
    for(int j = 0;j < N_phi_bins;j++){
      hist_N_Jpsi_accxeff_area_corrected_2pt6_HE -> SetBinContent(i+1,j+1,(hist_N_Jpsi_accxeff_corrected_2pt6_HE -> GetBinContent(i+1,j+1))/bin_area[i][j]);
      hist_N_Jpsi_accxeff_area_corrected_2pt6_HE -> SetBinError(i+1,j+1,(hist_N_Jpsi_accxeff_corrected_2pt6_HE -> GetBinError(i+1,j+1))/bin_area[i][j]);
      //printf("%i +- %i \n",(hist_N_Jpsi_accxeff_corrected_2pt6_HE -> GetBinContent(i+1,j+1))/bin_area[i][j],(hist_N_Jpsi_accxeff_corrected_2pt6_HE -> GetBinError(i+1,j+1))/bin_area[i][j]);
    }
  }

  //TH1D *projection_phi = (TH1D*) hist_N_Jpsi_accxeff_corrected_2pt6_HE -> ProjectionY("projection_phi");
  TH1D *projection_cost = (TH1D*) hist_N_Jpsi_accxeff_corrected_2pt6_HE -> ProjectionX("projection_cost");

  TCanvas *cccc = new TCanvas("cccc","cccc",10,10,600,600);
  //hist_N_Jpsi_accxeff_area_corrected_2pt6_HE -> Draw("COLZtext");
  //projection_phi -> Draw();
  projection_cost -> Draw();


  //============================================================================
  // POLARIZATION FIT
  //============================================================================

  //----------------------------------------------------------------------------
  // RANGE 0 < phi < PI ; -1 < cost < 1
  //----------------------------------------------------------------------------
  TH2D *hist_N_Jpsi_fit_HE = (TH2D*) hist_N_Jpsi_accxeff_area_corrected_2pt6_HE -> Clone("hist_N_Jpsi_fit_HE");

  TF2 *func_W_HE = new TF2("func_W_HE",Func_W,-0.6,0.6,0.502655,2.63894,4); // binning 2
  //TF2 *func_W_HE = new TF2("func_W_HE",Func_W,-0.8,0.8,0.942478,2.19911,4); // binning 1
  //TF2 *func_W_HE = new TF2("func_W_HE",Func_W,-0.8,0.8,0,PI,4);
  func_W_HE -> SetParameter(0,1000);
  func_W_HE -> SetParName(0,"N");
  func_W_HE -> SetParameter(1,0);
  func_W_HE -> SetParName(1,"#lambda_{#theta}");
  func_W_HE -> FixParameter(2,0);
  func_W_HE -> SetParName(2,"#lambda_{#phi}");
  func_W_HE -> FixParameter(3,0);
  func_W_HE -> SetParName(3,"#lambda_{#theta#phi}");
  //hist_N_Jpsi_accxeff_area_corrected_2pt6_HE -> Fit(func_W_HE,"RSI0LW");
  hist_N_Jpsi_accxeff_area_corrected_2pt6_HE -> Fit(func_W_HE,"RSI0");

  char title[100];
  TPaveText *t_fit_HE = new TPaveText(0.65,0.75,0.95,0.95,"brNDC");
  t_fit_HE -> SetFillColor(kWhite);
  sprintf(title,"N = %3.2f #pm %3.2f",func_W_HE -> GetParameter(0),func_W_HE -> GetParError(0));
  t_fit_HE -> AddText(title);
  sprintf(title,"#lambda_{#theta} = %f #pm %f",func_W_HE -> GetParameter(1),func_W_HE -> GetParError(1));
  t_fit_HE -> AddText(title);
  sprintf(title,"#lambda_{#phi} = %f #pm %f",func_W_HE -> GetParameter(2),func_W_HE -> GetParError(2));
  t_fit_HE -> AddText(title);
  sprintf(title,"#lambda_{#theta#phi} = %f #pm %f",func_W_HE -> GetParameter(3),func_W_HE -> GetParError(3));
  t_fit_HE -> AddText(title);

  TCanvas *c_fit_2pt6_HE = new TCanvas("c_fit_2pt6_HE","c_fit_2pt6_HE",4,132,1024,768);

  //TH3D *h_frame = new TH3D("h_frame","",100,-1,1,50,0,PI,100,20000,200000);
  //hist_N_Jpsi_fit_HE -> Draw("SURF1");
  //func_W_HE -> Draw("SURFsame");
  //h_frame -> Draw();
  //hist_N_Jpsi_fit_HE -> Draw("SURF1same");
  //func_W_HE -> Draw("SURF");
  //hist_N_Jpsi_fit_HE -> Draw("SURF1same");
  hist_N_Jpsi_fit_HE -> Draw("COLZtexterror");
  func_W_HE -> Draw("same");
  t_fit_HE -> Draw("same");

  cout << func_W_HE -> GetChisquare()/func_W_HE -> GetNDF() << endl;
  cout << func_W_HE -> Eval(0,0) << endl;

  double central_cost[N_cost_bins];
  double central_phi[N_cost_bins];

  for(int i = 0;i < N_cost_bins;i++){
    central_cost[i] = value_cost[i] + (TMath::Abs(value_cost[i+1] - value_cost[i])/2);
    for(int j = 0;j < N_phi_bins;j++){
      central_phi[j] = value_phi[j] + (TMath::Abs(value_phi[j+1] - value_phi[j])/2);
    }
  }

  TH2D *hist_ratio_func_histo = new TH2D("hist_ratio_func_histo","",N_cost_bins,value_cost,N_phi_bins,value_phi);

  for(int i = 0;i < N_cost_bins;i++){
    for(int j = 0;j < N_phi_bins;j++){
      hist_ratio_func_histo -> SetBinContent(i+1,j+1,hist_N_Jpsi_accxeff_area_corrected_2pt6_HE -> GetBinContent(i+1,j+1)/func_W_HE -> Eval(central_cost[i],central_phi[j]));
    }
  }

  TCanvas *c_histo_ratio = new TCanvas("c_histo_ratio","c_histo_ratio",20,20,600,600);
  hist_ratio_func_histo -> Draw("COLZtext");

  //============================================================================
  // CHECK THE TH1 PROJECTION
  //============================================================================

  //----------------------------------------------------------------------------
  // RANGE 0 < phi < PI ; -1 < cost < 1
  //----------------------------------------------------------------------------
  TH1D *projection_cost = (TH1D*) hist_N_Jpsi_fit_HE -> ProjectionX("projection_cost");

  TF1 *projection_cost_W = new TF1("projection_cost_W",Func_cost,-0.8,0.8,2);
  projection_cost_W -> SetParameter(0,func_W_HE -> GetParameter(0));
  projection_cost_W -> SetParameter(1,func_W_HE -> GetParameter(1));
  //projection_cost_W -> FixParameter(1,0.6);
  projection_cost -> Fit(projection_cost_W,"RSI0LW");

  TCanvas *c_projection_cost = new TCanvas("c_projection_cost","c_projection_cost",4,132,1024,768);
  projection_cost -> Draw();
  projection_cost_W -> Draw("same");
  //projection_cost_W -> Draw();
  //projection_cost -> Draw("same");

  TH1D *projection_phi = (TH1D*) hist_N_Jpsi_fit_HE -> ProjectionY("projection_phi");

  TF1 *projection_phi_W = new TF1("projection_phi_W",Func_phi,0.502655,2.63894,3);
  projection_phi_W -> SetParameter(0,func_W_HE -> GetParameter(0));
  projection_phi_W -> FixParameter(1,func_W_HE -> GetParameter(1)); // it has to be fixed
  projection_phi_W -> SetParameter(2,func_W_HE -> GetParameter(2));
  projection_phi -> Fit(projection_phi_W,"RSI0LW");

  TCanvas *c_projection_phi = new TCanvas("c_projection_phi","c_projection_phi",4,132,1024,768);
  projection_phi -> Draw();
  projection_phi_W -> Draw("same");

}
////////////////////////////////////////////////////////////////////////////////
double Func_W(double *x, double *par){

  double N = par[0];
  double L_th = par[1];
  double L_ph = par[2];
  double L_thph = par[3];

  double costh = x[0];
  double phi = x[1];

  double cosph = TMath::Cos(phi);
  double cos2ph = TMath::Cos(2*phi);

  double W =  (N/(3 + L_th))*(1 + (L_th + L_ph*cos2ph)*costh*costh + 2*L_thph*costh*cosph*TMath::Sqrt(1 - costh*costh) + L_ph*cos2ph);
  return W;
}
//------------------------------------------------------------------------------
double Func_cost(double *x, double *par){

  double N = par[0];
  double L_th = par[1];

  double costh = x[0];

  double W =  (N/(3 + L_th))*(1 + L_th*costh*costh);
  return W;
}
//------------------------------------------------------------------------------
double Func_phi(double *x, double *par){

  double N = par[0];
  double L_th = par[1];
  double L_ph = par[2];

  double phi = x[0];

  double cos2ph = TMath::Cos(2*phi);

  double W =  N*(1 + ((2*L_ph)/(3 + L_th))*cos2ph);
  return W;
}
