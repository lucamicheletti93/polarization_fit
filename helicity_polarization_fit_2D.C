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

  gSystem -> CompileMacro("../settings.h");
  gROOT -> ProcessLine(".x ../binning.C");

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
  TTree *N_Jpsi_tree = (TTree*) N_Jpsi_file -> Get("CB2_VWG");

  N_Jpsi_tree -> SetBranchAddress("N_Jpsi_HE",matrix_N_Jpsi_HE);
  N_Jpsi_tree -> SetBranchAddress("Stat_Jpsi_HE",matrix_stat_N_Jpsi_HE);

  for(int i = 0;i < N_Jpsi_tree -> GetEntries();i++){N_Jpsi_tree -> GetEntry(i);}

  /*for(int i = 0;i < N_cost_bins;i++){
    for(int j = 0;j < N_phi_bins;j++){
      cout << matrix_N_Jpsi_HE[i][j] << " ";
    }
    cout << endl;
  }*/

  TH2D *hist_N_Jpsi_2pt6_HE = new TH2D("hist_N_Jpsi_2pt6_HE","hist_N_Jpsi_2pt6_HE",N_cost_bins,value_cost,N_phi_bins,value_phi);

  TH2D *h_N_Jpsi_HE = new TH2D("h_N_Jpsi_HE","h_N_Jpsi_HE",N_cost_bins,value_cost,N_phi_bins,value_phi);
  h_N_Jpsi_HE -> GetXaxis() -> SetLabelSize(0);
  h_N_Jpsi_HE -> GetYaxis() -> SetLabelSize(0);

  for(int i = 0;i< N_cost_bins;i++){
    for(int j = 0;j < N_phi_bins;j++){
      hist_N_Jpsi_2pt6_HE -> SetBinContent(i+1,j+1,matrix_N_Jpsi_HE[i][j]);
      hist_N_Jpsi_2pt6_HE -> SetBinError(i+1,j+1,matrix_stat_N_Jpsi_HE[i][j]);
    }
  }

  printf("%f +- %f \n",hist_N_Jpsi_2pt6_HE -> GetBinContent(1,1),hist_N_Jpsi_2pt6_HE -> GetBinError(1,1));

  /*TCanvas *c_N_Jpsi_2pt6_HE = new TCanvas("c_N_Jpsi_2pt6_HE","c_N_Jpsi_2pt6_HE",4,132,1024,768);
  h_N_Jpsi_HE -> Draw();
  hist_N_Jpsi_2pt6_HE -> Draw("COLZtext");

  for(int i = 0;i < N_line_cost;i++){
    if(i < N_line_phi) line_phi[i] -> Draw("same");
    line_cost[i] -> Draw("same");
  }*/

  //============================================================================
  // ACCXEFF HISTOGRAMS
  //============================================================================

  TFile *accxeff_file = new TFile("/home/luca/cernbox/JPSI/JPSI_POLARIZATION/ANALYSIS/TWO_DIM_APPROACH/ACCXEFF/HISTOS_FOR_ACCXEFF/GIT_OUTPUT/accxeff.root","READ");
  TH2D *hist_accxeff_2pt6_HE = (TH2D*) accxeff_file -> Get("hist_accxeff_HE_2pt6_rebin");

  printf("%f +- %f \n",hist_accxeff_2pt6_HE -> GetBinContent(1,1),hist_accxeff_2pt6_HE -> GetBinError(1,1));

  //============================================================================
  // N_Jpsi CORRECTED FOR ACCXEFF DISTRIBUTION
  //============================================================================

  TH2D *hist_N_Jpsi_accxeff_corrected_2pt6_HE = new TH2D("hist_N_Jpsi_accxeff_corrected_2pt6_HE","hist_N_Jpsi_accxeff_corrected_2pt6_HE",N_cost_bins,value_cost,N_phi_bins,value_phi);
  hist_N_Jpsi_accxeff_corrected_2pt6_HE -> Divide(hist_N_Jpsi_2pt6_HE,hist_accxeff_2pt6_HE,1,1);

  printf("%f +- %f \n",hist_N_Jpsi_accxeff_corrected_2pt6_HE -> GetBinContent(1,1),hist_N_Jpsi_accxeff_corrected_2pt6_HE -> GetBinError(1,1));

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


  /*for(int i = 0;i< N_cost_bins;i++){
    cout << "{";
    for(int j = 0;j < N_phi_bins;j++){
      //cout << hist_N_Jpsi_accxeff_area_corrected_2pt6_HE -> GetBinContent(i+1,j+1) << " ";
      printf("%i,",hist_N_Jpsi_accxeff_area_corrected_2pt6_HE -> GetBinError(i+1,j+1));
    }
    cout << "}," << endl;
  }*/


  //============================================================================
  // POLARIZATION FIT
  //============================================================================

  //----------------------------------------------------------------------------
  // RANGE 0 < phi < PI ; -1 < cost < 1
  //----------------------------------------------------------------------------
  TH2D *hist_N_Jpsi_fit_HE = (TH2D*) hist_N_Jpsi_accxeff_area_corrected_2pt6_HE -> Clone("hist_N_Jpsi_fit_HE");

  TF2 *func_W_HE = new TF2("func_W_HE",Func_W,-0.5,0.5,0.95,2.2,4);
  func_W_HE -> SetParameter(0,1000);
  func_W_HE -> SetParName(0,"N");
  func_W_HE -> SetParameter(1,1);
  func_W_HE -> SetParName(1,"#lambda_{#theta}");
  func_W_HE -> SetParameter(2,1);
  func_W_HE -> SetParName(2,"#lambda_{#phi}");
  func_W_HE -> SetParameter(3,1);
  func_W_HE -> SetParName(3,"#lambda_{#theta#phi}");
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
  //hist_N_Jpsi_fit_HE -> Draw("SURF1");
  //func_W_HE -> Draw("SURFsame");
  hist_N_Jpsi_fit_HE -> Draw("COLZ");
  func_W_HE -> Draw("same");
  t_fit_HE -> Draw("same");

  return;

  //============================================================================
  // CHECK THE TH1 PROJECTION
  //============================================================================

  //----------------------------------------------------------------------------
  // RANGE 0 < phi < PI ; -1 < cost < 1
  //----------------------------------------------------------------------------
  TH1D *projection_cost = (TH1D*) hist_N_Jpsi_fit_HE -> ProjectionX("projection_cost");

  TF1 *projection_cost_W = new TF1("projection_cost_W",Func_cost,-1,1,2);
  projection_cost_W -> SetParameter(0,func_W_HE -> GetParameter(0));
  projection_cost_W -> SetParameter(1,func_W_HE -> GetParameter(1));
  projection_cost -> Fit(projection_cost_W,"RS0");

  TCanvas *c_projection_cost = new TCanvas("c_projection_cost","c_projection_cost",4,132,1024,768);
  projection_cost -> Draw();
  projection_cost_W -> Draw("same");
  //projection_cost_W -> Draw();
  //projection_cost -> Draw("same");

  TH1D *projection_phi = (TH1D*) hist_N_Jpsi_fit_HE -> ProjectionY("projection_phi");

  TF1 *projection_phi_W = new TF1("projection_phi_W",Func_phi,0,PI,3);
  projection_phi_W -> SetParameter(0,func_W_HE -> GetParameter(0));
  projection_phi_W -> FixParameter(1,func_W_HE -> GetParameter(1));
  projection_phi_W -> SetParameter(2,func_W_HE -> GetParameter(2));
  projection_phi -> Fit(projection_phi_W,"RS0");

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
