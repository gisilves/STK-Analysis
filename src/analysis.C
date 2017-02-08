#include <TROOT.h>
#include <TSystem.h>
#include <TF1.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TColor.h>
#include <TGraph.h>
#include <TString.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TRandom3.h>
#include <TPaveText.h>
#include <TGraphErrors.h>
#include <TClonesArray.h>
#include <TObject.h>
#include <TProfile.h>
#include <TMinuit.h>
#include <TLine.h>
#include <TVirtualFitter.h>

#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

#include <vector>
#include <string>
#include <sstream>
#include <utility>
#include <fstream>
#include <iostream>
#include <algorithm>

using namespace std;


#include "/home/gigi/Root/DmpSoftware/trunk/Event/include/DmpStkSiCluster.h"
#include "/home/gigi/Root/DmpSoftware/trunk/Event/include/DmpStkLadderAdc.hh"
#include "/home/gigi/Root/DmpSoftware/trunk/Event/include/Scintillator.h"
#include "/home/gigi/Root/DmpSoftware/trunk/Event/include/ScintillatorCluster.h"

bool debug = false;

int run_id = 0;
TString ion_name[16]	 = {	"run"	,"H"	,"He"	,"Li"	,"Be"	,"B"	,"C"	,"N"	,"O"	,"F"	,"Ne"	,"Na"	,"Mg"	,"Al"	,"Si"	,"P"    }; // run_id
double ion_mean_pos[32] = {	00	,-1	,-1	,-1	,-1	,-1	,-1	,-1	,-1	,-1	,-1	,-1	,-1	,-1	,-1	,-1     ,  // 0
				155     ,-1     ,20.89  ,29.97  ,40.4   ,49.78  ,58.63  ,66     ,73.48  ,79.76  ,86.51  ,92.2   ,97.55  ,103.27 ,107.83 ,111.91 }; // 40

//---------------------------------------- track fitting -------------------------------------//

// define a track function - straigt line
double func(float z, double *par)
{
  double value = par[0]*z + par[1];
  return value;
}

// define global data points
int ndatapoints = 6;
vector<double> pos_xy, pos_xy_error, pos_z;

void init_data_points(int n)
{
  ndatapoints = n;
  pos_xy.clear();
  pos_xy.assign(ndatapoints,0);
  pos_z.clear();
  pos_z.assign(ndatapoints,0);
  pos_xy_error.clear();
  pos_xy_error.assign(ndatapoints,0);
}

// define FCN
double chi2;
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  Double_t chisq = 0;
  Double_t delta;
  Int_t i;
  //calculate chisquare 
  for (int i=0;i<ndatapoints; i++) {
    //delta  = (pos_xy[i]-func(pos_z[i],par))/pos_xy_error[i];
    delta  = (pos_xy.at(i)-func(pos_z.at(i),par))/pos_xy_error.at(i);
    chisq += delta*delta; }
  f = chisq;
  chi2 = chisq;
}

// execute the fit
vector<double> fit_track()
{
  if (debug) printf("\n\n@@@@@@@@@@@@@@@ FIT TRACK @@@@@@@@@@@@@@ \n\n ");

  //initialize TMinuit with a maximum of 5 params
  TMinuit *gMinuit = new TMinuit(5);  
  gMinuit->SetFCN(fcn);

  if (debug) gMinuit->SetPrintLevel(1);
  else  gMinuit->SetPrintLevel(-1);

  Double_t arglist[10];
  Int_t ierflg = 0;
  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

  // Set starting values and step sizes for parameters
  static Double_t vstart[2] = {   -0.1    ,  99.0 };
  static Double_t vlow[2]   = {  -10.0    ,   0.0 };
  static Double_t vup[2]    = {   10.0    , 100.0 };
  static Double_t step[2]   = {    0.0001 ,   0.1 };
  gMinuit->mnparm(0, "a1", vstart[0], step[0], vlow[0], vup[0], ierflg);
  gMinuit->mnparm(1, "a2", vstart[1], step[1], vlow[1], vup[1], ierflg);

  // Now ready for minimization step
  arglist[0] = 500;
  arglist[1] = 1.;
  gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

  // Get fit results
  double slope;
  double slope_err;
  double intersec;
  double intersec_err;
  gMinuit->GetParameter(0,slope,slope_err);
  gMinuit->GetParameter(1,intersec,intersec_err);

  // Return results
  vector<double> results;
  results.push_back(slope);        // 0
  results.push_back(slope_err);    // 1
  results.push_back(intersec);     // 2
  results.push_back(intersec_err); // 3
  results.push_back(chi2);         // 4
  return results;
}

//------------------------------------------- main function -------------------------------------//

void analysis(
	      int run_id = 0,
	      int neventsforalignment = 30000,
	      int neventstoprocess = 100000,
	      int firsteventtodisplay = 0,
	      int neventstodisplay = 30,
	      int stop = 0,
	      bool dbg = false
	      )
{  
  TString run_id_str = TString::Format("%04d",(int)ion_mean_pos[run_id*16]);
  TString filename = "/home/gigi/tmp/processing0/run_"+run_id_str+"_MERGED_v2.root";

  if (dbg) debug = true;

  if (!debug) {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(11); 
  } else {
    gStyle->SetOptStat(11111111);
    gStyle->SetOptFit(1111); 
  }
  gSystem->Load("/home/gigi/Root/DmpSoftware/trunk/Event/libDmpEvent.so");

  TString term_cmd = "mkdir -pv "+TString::Format("plots/run%d",(int)ion_mean_pos[run_id*16]);
  system(term_cmd.Data());

  TFile* f = new TFile(filename,"READ");
  TTree* t = (TTree*)f->Get("CollectionTree");

  // Register STK collections
  TClonesArray* stkclusters = new TClonesArray("DmpStkSiCluster");
  t->SetBranchAddress("StkClusterCollection",&stkclusters);

  TClonesArray* stkladderadc = new TClonesArray("DmpStkLadderAdc");
  t->SetBranchAddress("DmpStkLadderAdcCollection", &stkladderadc);

  // Register SCINTILLATOR collections
  TClonesArray* scintillators = new TClonesArray("Scintillator");  
  t->SetBranchAddress("ScintillatorAdcCollection",&scintillators);
  TClonesArray* scintclusters = new TClonesArray("ScintillatorCluster");   
  t->SetBranchAddress("ScintillatorClusterCollection",&scintclusters);       
    
  //----------------------------------------- telescope and DUT planes positions
  double stk_clusters_z[20] = { 14.2, 14.2, 81.6, 81.6, 87.0, 87.0, 112.0, 113.0, 116.0, 117.0, 120.0, 121.0, 124.0, 125.0, 146.4, 146.4, 152.5, 152.5, 200.0, 200.0};
  int stk_ladder_id[20]    =  {   31,   30,   25,   24,   72,   73,   103,   102,    97,    96,   127,   126,   121,   120,    78,    79,     0,     1,     6,     7};
  //                               0     1     2     3     4     5      6      7      8      9     10     11     12     13     14     15     16     17     18     19
  //                               Y     X     Y     X     Y     X      Y      Y      Y      Y      Y      Y      Y      Y      Y      X      Y      X      Y      X
  //                           ------------------------------------                                                        ------------------------------------------
  //                           |            telescope             | <<<<<<<<<<<<<<<<<<<<<<<<<< DUT >>>>>>>>>>>>>>>>>>>>>>> |                telescope               |
  //                           ------------------------------------                                                        ------------------------------------------
  int x_plane_id[6 ] = {1,3,5,15,17,19};
  int y_plane_id[14] = {0,2,4,6,7,8,9,10,11,12,13,14,16,18};

  /////////////////////////////////////////// alignment ///////////////////////////////////////////
  
  TH2F *non_aligned_corr[18];
  for (int i=0; i<18; i++) non_aligned_corr[i] = new TH2F("non_aligned_corr_"+TString::Format("%d",i), "non_aligned_corr_"+TString::Format("%d",i), 200,0,100,200,0,100);

  // Event LOOP
  if (neventsforalignment==-1 || neventsforalignment > t->GetEntries()) neventsforalignment = t->GetEntries();
  for(int entry=0; entry<neventsforalignment; entry++){
   
    t->GetEntry(entry);

    double stkcluster_maxehitposition[20];
    double stkcluster_maxe[20];
    for (int i=0; i<20; i++) {
      stkcluster_maxehitposition[i] = -999.9;
      stkcluster_maxe[i] = -999.9; }

    // Loop over STK clusters in event
    int nstkclusters = stkclusters->GetLast()+1;
    // printf("nstkclusters = %d\n",nstkclusters);
    for(int i=0; i<nstkclusters; i++){
      DmpStkSiCluster* cluster = (DmpStkSiCluster*)stkclusters->ConstructedAt(i);
      double stkclusterenergy = cluster->getEnergy();
      double stkclusterposition = cluster->getLadderCentroid();
      int stkludderid = cluster->getLadderHardware();

      //----------------------------- find maximum energy cluster position in every ladder
      for (int i=0; i<20; i++) {
      if (stkludderid==stk_ladder_id[i] && stkclusterenergy > stkcluster_maxe[i]) { 
	stkcluster_maxe[i] = stkclusterenergy; 
	stkcluster_maxehitposition[i] = stkclusterposition; }}
    }// loop over STK clusters

    // Y telescope planes
    non_aligned_corr[0 ]->Fill( 93-stkcluster_maxehitposition[0 ], 93-stkcluster_maxehitposition[2 ]); // 31  -  25
    non_aligned_corr[1 ]->Fill( 93-stkcluster_maxehitposition[2 ],    stkcluster_maxehitposition[4 ]); // 25  -  72
    non_aligned_corr[2 ]->Fill(    stkcluster_maxehitposition[4 ],    stkcluster_maxehitposition[14]); // 72  -  78
    non_aligned_corr[3 ]->Fill(    stkcluster_maxehitposition[14],    stkcluster_maxehitposition[16]); // 78  -   0
    non_aligned_corr[4 ]->Fill(    stkcluster_maxehitposition[16],    stkcluster_maxehitposition[18]); //  0  -   6 
    // X telescope planes
    non_aligned_corr[5 ]->Fill( 93-stkcluster_maxehitposition[1 ], 93-stkcluster_maxehitposition[3 ]); // 30  -  24
    non_aligned_corr[6 ]->Fill( 93-stkcluster_maxehitposition[3 ],    stkcluster_maxehitposition[5 ]); // 24  -  73
    non_aligned_corr[7 ]->Fill(    stkcluster_maxehitposition[5 ],    stkcluster_maxehitposition[15]); // 73  -  79
    non_aligned_corr[8 ]->Fill(    stkcluster_maxehitposition[15],    stkcluster_maxehitposition[17]); // 79  -   1
    non_aligned_corr[9 ]->Fill(    stkcluster_maxehitposition[17],    stkcluster_maxehitposition[19]); //  1  -   7
    // Y DUT planes //NOTE: ladder 73 as reference
    non_aligned_corr[10]->Fill(    stkcluster_maxehitposition[5],     stkcluster_maxehitposition[6 ]); // 73  - 103
    non_aligned_corr[11]->Fill(    stkcluster_maxehitposition[5],     stkcluster_maxehitposition[7 ]); // 73  - 102
    non_aligned_corr[12]->Fill(    stkcluster_maxehitposition[5],     stkcluster_maxehitposition[8 ]); // 73  -  97
    non_aligned_corr[13]->Fill(    stkcluster_maxehitposition[5],     stkcluster_maxehitposition[9 ]); // 73  -  96
    non_aligned_corr[14]->Fill(    stkcluster_maxehitposition[5],     stkcluster_maxehitposition[10]); // 73  - 127
    non_aligned_corr[15]->Fill(    stkcluster_maxehitposition[5],     stkcluster_maxehitposition[11]); // 73  - 126
    non_aligned_corr[16]->Fill(    stkcluster_maxehitposition[5],     stkcluster_maxehitposition[12]); // 73  - 121
    non_aligned_corr[17]->Fill(    stkcluster_maxehitposition[5],     stkcluster_maxehitposition[13]); // 73  - 120
    
  }// end of event loop
  
  // plot correlation in telescope planes
  TCanvas *non_aligned_corr_telescope_canv = new TCanvas("non_aligned_corr_telescope_canv","non_aligned_corr_telescope_canv",1800,1000);
  non_aligned_corr_telescope_canv->Divide(5,2);
  for (int i=0; i<10; i++) { non_aligned_corr_telescope_canv->cd(i+1); non_aligned_corr[i]->Draw("COLZ"); }
  non_aligned_corr_telescope_canv->SaveAs(TString::Format("plots/run%d/non_aligned_corr_telescope",(int)ion_mean_pos[run_id*16])+"_canv.png");

  // plot correlations in DUT planes
  TCanvas *non_aligned_corr_dut_canv = new TCanvas("non_aligned_corr_dut_canv","non_aligned_corr_dut_canv",1800,1000);
  non_aligned_corr_dut_canv->Divide(4,2);
  for (int i=10; i<18; i++) { non_aligned_corr_dut_canv->cd(i-9); non_aligned_corr[i]->Draw("COLZ"); }
  non_aligned_corr_telescope_canv->SaveAs(TString::Format("plots/run%d/non_aligned_corr_telescope",(int)ion_mean_pos[run_id*16])+"_canv.png");
  //_canv->SaveAs(TString::Format("plots/run%d/",(int)ion_mean_pos[run_id*16])+"_canv.png");

  if (stop==1) return;

  //-------- correlations
  for (int i=0; i<18; i++) {
    non_aligned_corr[i]->Sumw2();
    non_aligned_corr[i]->Scale(1.0/non_aligned_corr[i]->Integral()); }

  // remove background hits
  for (int i=1; i<=200; i++) {
    for (int j=1; j<=200; j++) {
      for (int k=0; k<18; k++) {
	if (non_aligned_corr[k]->GetBinContent(i,j)<0.00015) non_aligned_corr[k]->SetBinContent(i,j,0); }}}
  
  // extract alignment parameters
  TF1 *fitf = new TF1("fitf","[0]+[1]*x");
  double p0[18];
  double p1[18];
  TProfile *non_aligned_corr_prof[18];
  for (int k=0; k<18; k++) {
    non_aligned_corr_prof[k] = non_aligned_corr[k]->ProfileX();    
    non_aligned_corr_prof[k]->Fit("fitf","NQ0","",40,50);
    p0[k] = fitf->GetParameter(0); p1[k] = fitf->GetParameter(1); }

  //----------------------------------- profile plots
  TCanvas *non_aligned_tlp_prof_canv;
  TCanvas *non_aligned_dut_prof_canv;
  
  //-------- plot correlations in telelscope planes
  non_aligned_tlp_prof_canv = new TCanvas("non_aligned_tlp_prof_canv","non_aligned_tlp_prof_canv",1800,1000);
  non_aligned_tlp_prof_canv->Divide(5,2);
  for (int k=0; k<10; k++) {
    non_aligned_tlp_prof_canv->cd(k+1);  
    non_aligned_corr_prof[k] = non_aligned_corr[k]->ProfileX();    
    non_aligned_corr_prof[k]->Fit("pol1","Q","",35,50); 
    non_aligned_corr_prof[k]->GetYaxis()->SetRangeUser(0,100);     
    non_aligned_tlp_prof_canv->cd(k+1)->SetGridx();  
    non_aligned_tlp_prof_canv->cd(k+1)->SetGridy(); }
  
  //-------- plot correlations in DUT planes
  non_aligned_dut_prof_canv = new TCanvas("non_aligned_dut_prof_canv","non_aligned_dut_prof_canv",1800,1000);
  non_aligned_dut_prof_canv->Divide(4,2);
  for (int k=10; k<18; k++) {
    non_aligned_dut_prof_canv->cd(k-9);  
    non_aligned_corr_prof[k] = non_aligned_corr[k]->ProfileX();    
    non_aligned_corr_prof[k]->Fit("pol1","Q","",30,60); 
    non_aligned_corr_prof[k]->GetYaxis()->SetRangeUser(0,100);     
    non_aligned_dut_prof_canv->cd(k-9)->SetGridx();  
    non_aligned_dut_prof_canv->cd(k-9)->SetGridy(); }
  non_aligned_tlp_prof_canv->SaveAs(TString::Format("plots/run%d/non_aligned_tlp_prof",(int)ion_mean_pos[run_id*16])+"_canv.png");
  non_aligned_dut_prof_canv->SaveAs(TString::Format("plots/run%d/non_aligned_dut_prof",(int)ion_mean_pos[run_id*16])+"_canv.png");

  if (stop==2) return;

  //////////////////////////////////////// end of alignment ///////////////////////////////////////////
  

  // //////////////////////////////////////////////// analysis ///////////////////////////////////////////

  // Correlation plots after alignment
  TH2F *aligned_corr[18];
  for (int i=0; i<18; i++) aligned_corr[i] = new TH2F("aligned_corr_"+TString::Format("%d",i), "aligned_corr_"+TString::Format("%d",i), 200,0,100,200,0,100);

  //----------------------------------- event displays  
  TH2F *stk_hits_XZ[neventstodisplay];
  TH2F *stk_hits_YZ[neventstodisplay];
  for (int i=0; i<neventstodisplay; i++){
    stk_hits_XZ[i] = new TH2F("stk_hits_XZ_event_"+TString::Format("%d",firsteventtodisplay+i),
  			      "stk_hits_XZ_event_"+TString::Format("%d",firsteventtodisplay+i),
  			      200,0,220,200,0,100);
    stk_hits_YZ[i] = new TH2F("stk_hits_YZ_event_"+TString::Format("%d",firsteventtodisplay+i),
  			      "stk_hits_YZ_event_"+TString::Format("%d",firsteventtodisplay+i),
  			      200,0,220,200,0,100); }


  int event_display_i = 0;
  TLine *x_track_6_line[neventstodisplay];
  TLine *x_track_3a_line[neventstodisplay];
  TLine *x_track_3b_line[neventstodisplay];
  TLine *y_track_6_line[neventstodisplay];
  TLine *y_track_3a_line[neventstodisplay];
  TLine *y_track_3b_line[neventstodisplay];
  TCanvas *evtdisp_canv[neventstodisplay];
  for (int i=0; i<neventstodisplay; i++) {
    evtdisp_canv[i] = new TCanvas("evtdisp_canv"+TString::Format("%d",firsteventtodisplay+i),
  			  "evtdisp_canv"+TString::Format("%d",firsteventtodisplay+i),
  			  1800,1000);
    evtdisp_canv[i]->Divide(1,2); 
    x_track_6_line[i] = new TLine();
    x_track_3a_line[i] = new TLine();
    x_track_3b_line[i] = new TLine();
    y_track_6_line[i] = new TLine();
    y_track_3a_line[i] = new TLine();
    y_track_3b_line[i] = new TLine();
  }

  //----------------------------------------- Chi2 plots
  // X tracks
  TH1F *chi2_x_track_3a = new TH1F("chi2_x_track_3a","chi2_x_track_3a",300,-8,4);
  TH1F *chi2_x_track_3b = new TH1F("chi2_x_track_3b","chi2_x_track_3b",300,-8,4);
  TH1F *chi2_x_track_6  = new TH1F("chi2_x_track_6", "chi2_x_track_6", 300,-8,4);
  TH2F *chi2_x_track_3a_vs_3b = new TH2F("chi2_x_track_3a_vs_3b","chi2_x_track_3a_vs_3b",300,-8,4,300,-8,4);
  TH2F *chi2_x_track_3a_vs_6  = new TH2F("chi2_x_track_3a_vs_6", "chi2_x_track_3a_vs_6", 300,-8,4,300,-8,4);
  TH2F *chi2_x_track_3b_vs_6  = new TH2F("chi2_x_track_3b_vs_6", "chi2_x_track_3b_vs_6", 300,-8,4,300,-8,4);
  // Y tracks
  TH1F *chi2_y_track_3a = new TH1F("chi2_y_track_3a","chi2_y_track_3a",300,-8,4);
  TH1F *chi2_y_track_3b = new TH1F("chi2_y_track_3b","chi2_y_track_3b",300,-8,4);
  TH1F *chi2_y_track_6  = new TH1F("chi2_y_track_6", "chi2_y_track_6", 300,-8,4);
  TH2F *chi2_y_track_3a_vs_3b = new TH2F("chi2_y_track_3a_vs_3b","chi2_y_track_3a_vs_3b",300,-8,4,300,-8,4);
  TH2F *chi2_y_track_3a_vs_6  = new TH2F("chi2_y_track_3a_vs_6", "chi2_y_track_3a_vs_6", 300,-8,4,300,-8,4);
  TH2F *chi2_y_track_3b_vs_6  = new TH2F("chi2_y_track_3b_vs_6", "chi2_y_track_3b_vs_6", 300,-8,4,300,-8,4);

  //----------------------------------------- Cluster energy plots
  /*TH1F *cluster_q_x_track[6];
  TH1F *cluster_q_y_track[6];
  for (int i=0; i<6; i++) {
    cluster_q_x_track[i] = new TH1F("cluster_q_x_track_lad"+TString::Format("%d",i),
  				    "cluster_q_x_track_lad"+TString::Format("%d",i),300,0,200);
    cluster_q_y_track[i] = new TH1F("cluster_q_y_track_lad"+TString::Format("%d",i),
    "cluster_q_y_track_lad"+TString::Format("%d",i),300,0,200); }*/
  TH1F *cluster_q_track[20];
  for (int i=0; i<20; i++) {
    cluster_q_track[i] = new TH1F("cluster_q_track_lad"+TString::Format("%d",i),
  				  "cluster_q_track_lad"+TString::Format("%d",i),300,0,200);}

  TH1F *cluster_z_dut_track = new TH1F("cluster_z_dut_track_lad",
  				       "cluster_z_dut_track_lad",600,0,20);
  
  TH1F *cluster_adc_track[20];
  for (int i=0; i<20; i++) {
    cluster_adc_track[i] = new TH1F("cluster_adc_track_lad"+TString::Format("%d",i),
  				    "cluster_adc_track_lad"+TString::Format("%d",i),300,0,10000);}
  
  
  //----------------------------------------- Occupancy plots
  TH1F *occupancy_plane[20];
  for (int i=0; i<20; i++) occupancy_plane[i] = new TH1F("occupancy_plane_lad"+TString::Format("%d",i),
  							 "occupancy_plane_lad"+TString::Format("%d",i),
  							 2000,-15,100);
  TH2F *energy_occupancy_plane[20];
  for (int i=0; i<20; i++) energy_occupancy_plane[i] = new TH2F("energy_occupancy_plane_lad"+TString::Format("%d",i),
  								"energy_occupancy_plane_lad"+TString::Format("%d",i),
  								2000,-15,100,100,1,5);
  TH2F *nstr_occupancy_plane[20];
  for (int i=0; i<20; i++) nstr_occupancy_plane[i] = new TH2F("nstr_occupancy_plane_lad"+TString::Format("%d",i),
  							      "nstr_occupancy_plane_lad"+TString::Format("%d",i),
  							      2000,-15,100,40,0.5,40.5);
  
  //------------------------------------------- cluster multiplicities
  TH1F *cluster_mult_plane[20];
  for (int i=0; i<20; i++) cluster_mult_plane[i] = new TH1F("cluster_mult_plane_lad"+TString::Format("%d",i),
  							    "cluster_mult_plane_lad"+TString::Format("%d",i),
  							    50,0.5,50.5);

  TH2F *maxq_vs_cluster_mult_plane[20];
  for (int i=0; i<20; i++) maxq_vs_cluster_mult_plane[i] = new TH2F("maxq_vs_cluster_mult_plane_lad"+TString::Format("%d",i),
  								    "maxq_vs_cluster_mult_plane_lad"+TString::Format("%d",i),
  								    300,0,200,50,0.5,50.5);
  
  //------------------------------------------- secondary clusters
  TH1F *secondary_clusters_energy_ratio_plane[20];
  for (int i=0; i<20; i++) secondary_clusters_energy_ratio_plane[i] = new TH1F("secondary_clusters_energy_ratio_plane_lad"+TString::Format("%d",i),
  									       "secondary_clusters_energy_ratio_plane_lad"+TString::Format("%d",i),
  									       200,-3,2);
  
  //------------------------------------------- cluster position residuals
  double track_ladder_pos[20];
  TH1F *cluster_track_residuals_plane[20];
  for (int i=0; i<20; i++) {
    track_ladder_pos[i] = -999.99;
    cluster_track_residuals_plane[i] = new TH1F("cluster_track_residuals_plane_lad"+TString::Format("%d",i),
  						"cluster_track_residuals_plane_lad"+TString::Format("%d",i),
  						300,-1,1); }
  
  //------------------------------------------- counters
  int n_single_particle_trigger_events = 0;
  int n_good_events = 0;
  int n_fragmentation_events = 0;
  int n_good_upsteam_events = 0;

  //------------------------------------------- scintillator clusters
  TH1F *scint_cluster_q[7];
  for (int i=0; i<7; i++) scint_cluster_q[i] = new TH1F("scint_cluster_q_ch"+TString::Format("%d",i),
  							"scint_cluster_q_ch"+TString::Format("%d",i),
  							600,0,400);
  TH1F *scint_cluster_pos[7];
  for (int i=0; i<7; i++) scint_cluster_pos[i] = new TH1F("scint_cluster_pos_ch"+TString::Format("%d",i),
  							  "scint_cluster_pos_ch"+TString::Format("%d",i),
  							  300,0,4000);
  TH1F *scint_cluster_width[7];
  for (int i=0; i<7; i++) scint_cluster_width[i] = new TH1F("scint_cluster_width_ch"+TString::Format("%d",i),
  							    "scint_cluster_width_ch"+TString::Format("%d",i),
  							    300,0,600);
  TH2F *stkmaxq_vs_scint_cluster_q[7];
  for (int i=0; i<7; i++) stkmaxq_vs_scint_cluster_q[i] = new TH2F("stkmaxq_vs_scint_cluster_q_ch"+TString::Format("%d",i),
  								   "stkmaxq_vs_scint_cluster_q_ch"+TString::Format("%d",i),
  								   300,0,200,600,0,400);
  //-------------------------------------------- tracking plots
  TH1F *delta_phi_upstream_downstream_x_tracks = new TH1F("delta_phi_upstream_downstream_x_tracks","delta_phi_upstream_downstream_x_tracks",200,-2,2);
  TH1F *delta_phi_upstream_downstream_y_tracks = new TH1F("delta_phi_upstream_downstream_y_tracks","delta_phi_upstream_downstream_y_tracks",200,-2,2);
  
  /////////////////////////////////// Event LOOP //////////////////////////////////

  if (neventstoprocess==-1 || neventstoprocess > t->GetEntries()) neventstoprocess = t->GetEntries();
  for(int entry=0; entry<neventstoprocess; entry++){
   
    t->GetEntry(entry);

    double stkcluster_maxehitposition[20];
    double stkcluster_maxe[20];
    int stknclusters[20];
    for (int i=0; i<20; i++) {
      stkcluster_maxehitposition[i] = -999.9;
      stkcluster_maxe[i] = -999.9; 
      stknclusters[i] = 0;}

    vector<double> stk_clusters_pos[20];
    vector<double> stk_clusters_e[20];
    
    //////////////////////////// Loop over STK clusters ////////////////////////////

    int nstkclusters = stkclusters->GetLast()+1;
    // printf("nstkclusters = %d\n",nstkclusters);
    for(int i=0; i<nstkclusters; i++){

      DmpStkSiCluster* cluster = (DmpStkSiCluster*)stkclusters->ConstructedAt(i);
      double stkclusterenergy = cluster->getEnergy();
      int stkclusternstr = cluster->getNstrip();         
      double stkclusterposition = cluster->getLadderCentroid();
      int stkludderid = cluster->getLadderHardware();


      //------------------------- reclculate cluster positions
      // XZ telescope planes
      if (stkludderid==stk_ladder_id[0 ]) stkclusterposition =  93.-stkclusterposition; 
      if (stkludderid==stk_ladder_id[2 ]) stkclusterposition = (93.-stkclusterposition-p0[0])/p1[0];
      if (stkludderid==stk_ladder_id[4 ]) stkclusterposition =    ((stkclusterposition-p0[1])/p1[1]-p0[0])/p1[0];
      if (stkludderid==stk_ladder_id[14]) stkclusterposition =   (((stkclusterposition-p0[2])/p1[2]-p0[1])/p1[1]-p0[0])/p1[0];
      if (stkludderid==stk_ladder_id[16]) stkclusterposition =  ((((stkclusterposition-p0[3])/p1[3]-p0[2])/p1[2]-p0[1])/p1[1]-p0[0])/p1[0];
      if (stkludderid==stk_ladder_id[18]) stkclusterposition = (((((stkclusterposition-p0[4])/p1[4]-p0[3])/p1[3]-p0[2])/p1[2]-p0[1])/p1[1]-p0[0])/p1[0];
      // YZ telescope planes                     
      if (stkludderid==stk_ladder_id[1 ]) stkclusterposition =  93.-stkclusterposition;
      if (stkludderid==stk_ladder_id[3 ]) stkclusterposition = (93.-stkclusterposition-p0[5])/p1[5];
      if (stkludderid==stk_ladder_id[5 ]) stkclusterposition =    ((stkclusterposition-p0[6])/p1[6]-p0[5])/p1[5];
      if (stkludderid==stk_ladder_id[15]) stkclusterposition =   (((stkclusterposition-p0[7])/p1[7]-p0[6])/p1[6]-p0[5])/p1[5];
      if (stkludderid==stk_ladder_id[17]) stkclusterposition =  ((((stkclusterposition-p0[8])/p1[8]-p0[7])/p1[7]-p0[6])/p1[6]-p0[5])/p1[5];
      if (stkludderid==stk_ladder_id[19]) stkclusterposition = (((((stkclusterposition-p0[9])/p1[9]-p0[8])/p1[8]-p0[7])/p1[7]-p0[6])/p1[6]-p0[5])/p1[5];
      // YZ DUT planes                        
      if (stkludderid==stk_ladder_id[6 ]) stkclusterposition =   (((stkclusterposition-p0[10])/p1[10]-p0[6])/p1[6]-p0[5])/p1[5];
      if (stkludderid==stk_ladder_id[7 ]) stkclusterposition =   (((stkclusterposition-p0[11])/p1[11]-p0[6])/p1[6]-p0[5])/p1[5];
      if (stkludderid==stk_ladder_id[8 ]) stkclusterposition =   (((stkclusterposition-p0[12])/p1[12]-p0[6])/p1[6]-p0[5])/p1[5];
      if (stkludderid==stk_ladder_id[9 ]) stkclusterposition =   (((stkclusterposition-p0[13])/p1[13]-p0[6])/p1[6]-p0[5])/p1[5];
      if (stkludderid==stk_ladder_id[10]) stkclusterposition =   (((stkclusterposition-p0[14])/p1[14]-p0[6])/p1[6]-p0[5])/p1[5];
      if (stkludderid==stk_ladder_id[11]) stkclusterposition =   (((stkclusterposition-p0[15])/p1[15]-p0[6])/p1[6]-p0[5])/p1[5];
      if (stkludderid==stk_ladder_id[12]) stkclusterposition =   (((stkclusterposition-p0[16])/p1[16]-p0[6])/p1[6]-p0[5])/p1[5];
      if (stkludderid==stk_ladder_id[13]) stkclusterposition =   (((stkclusterposition-p0[17])/p1[17]-p0[6])/p1[6]-p0[5])/p1[5];

      //------------------------- fill occupancy plots
      for (int i=0; i<20; i++) {
  	if (stkludderid==stk_ladder_id[i]) {
  	  occupancy_plane[i]->Fill(stkclusterposition);  
  	  energy_occupancy_plane[i]->Fill(stkclusterposition,log10(stkclusterenergy)); 
  	  nstr_occupancy_plane[i]->Fill(stkclusterposition,stkclusternstr); }}
      
      //------------------------- fill vectors
      for (int i=0; i<20; i++) {
  	if (stkludderid==stk_ladder_id[i]) {
  	  stk_clusters_pos[i].push_back(stkclusterposition); }}
      
      //-------------------------- fill cluster energy
      for (int i=0; i<20; i++) {
  	if (stkludderid==stk_ladder_id[i]) {
  	  stk_clusters_e[i].push_back(stkclusterenergy); }}

      //----------------------------- find maximum energy cluster position in every ladder
      for (int i=0; i<20; i++) {
  	if (stkludderid==stk_ladder_id[i] && stkclusterenergy > stkcluster_maxe[i]) { 
  	  stkcluster_maxe[i] = stkclusterenergy;
  	  stkcluster_maxehitposition[i] = stkclusterposition; }}

      //------------------------- fill cluster multiplicities arrays
      for (int i=0; i<20; i++) {
  	if (stkludderid==stk_ladder_id[i]){ stknclusters[i]++; }}

    }// Loop over STK clusters

    //-------------------- find events with particles before or after a trigger
    int n_scint_clusters = scintclusters->GetLast()+1;
    bool non_trigger_particle = false;
    double scint_cluster_charge[4];
    bool saturated_signal[4] = {false,false,false,false};
    if (n_scint_clusters!=2) continue;
    // if (non_trigger_particle) continue;
    
    // Loop over reconstructed SCINTILLATOR Clusters
    for (int i=0; i<n_scint_clusters; i++) {
      
      ScintillatorCluster* cluster = (ScintillatorCluster*)scintclusters->ConstructedAt(i);
      int scintid = cluster->GetScintillatorID(); 
      int scintclusterfirstsampleid = cluster->GetPosition(0);
      int scintclustersize = cluster->GetSize();
      scint_cluster_charge[scintid] = sqrt(cluster->GetSignalTotal());

      scint_cluster_pos[scintid]->Fill(scintclusterfirstsampleid);

    }// loop over scintillator clusters
    

    //---------------------- calculate secondary clusters energy ratio
    double secondary_clusters_energy_ratio[20];
    for (int i=0; i<20; i++){
      double total_e = 0;
      for (uint j=0; j<stk_clusters_e[i].size(); j++) {
  	total_e += stk_clusters_e[i].at(j); }
      double e_frac = (total_e-stkcluster_maxe[i])/stkcluster_maxe[i];///(stk_clusters_e[i].size()-1);
      secondary_clusters_energy_ratio[i] = e_frac; }
    
    //----------------------- fill aligned correlation histograms 
    if (neventstoprocess >= neventsforalignment) {

      // X telescope planes
      aligned_corr[0 ]->Fill( stkcluster_maxehitposition[0 ], stkcluster_maxehitposition[2 ]); // 0   -   6
      aligned_corr[1 ]->Fill( stkcluster_maxehitposition[2 ], stkcluster_maxehitposition[4 ]); // 6   -  72
      aligned_corr[2 ]->Fill( stkcluster_maxehitposition[4 ], stkcluster_maxehitposition[14]); // 72  -  96
      aligned_corr[3 ]->Fill( stkcluster_maxehitposition[14], stkcluster_maxehitposition[16]); // 96  - 126
      aligned_corr[4 ]->Fill( stkcluster_maxehitposition[16], stkcluster_maxehitposition[18]); // 126 - 120
      // Y telescope planes
      aligned_corr[5 ]->Fill( stkcluster_maxehitposition[1 ], stkcluster_maxehitposition[3 ]); // 1   -   7
      aligned_corr[6 ]->Fill( stkcluster_maxehitposition[3 ], stkcluster_maxehitposition[5 ]); // 7   -  73
      aligned_corr[7 ]->Fill( stkcluster_maxehitposition[5 ], stkcluster_maxehitposition[15]); // 73  -  97
      aligned_corr[8 ]->Fill( stkcluster_maxehitposition[15], stkcluster_maxehitposition[17]); // 97  - 127
      aligned_corr[9 ]->Fill( stkcluster_maxehitposition[17], stkcluster_maxehitposition[19]); // 127 - 121
      // Y DUT planes
      aligned_corr[10]->Fill( stkcluster_maxehitposition[5],  stkcluster_maxehitposition[6 ]); // 73  -  78
      aligned_corr[11]->Fill( stkcluster_maxehitposition[5],  stkcluster_maxehitposition[7 ]); // 73  -  79
      aligned_corr[12]->Fill( stkcluster_maxehitposition[5],  stkcluster_maxehitposition[8 ]); // 73  - 175
      aligned_corr[13]->Fill( stkcluster_maxehitposition[5],  stkcluster_maxehitposition[9 ]); // 73  - 102
      aligned_corr[14]->Fill( stkcluster_maxehitposition[5],  stkcluster_maxehitposition[10]); // 73  - 103
      aligned_corr[15]->Fill( stkcluster_maxehitposition[5],  stkcluster_maxehitposition[11]); // 73  - 169
      aligned_corr[16]->Fill( stkcluster_maxehitposition[5],  stkcluster_maxehitposition[12]); // 73  - 168
      aligned_corr[17]->Fill( stkcluster_maxehitposition[5],  stkcluster_maxehitposition[13]);}// 73  - 174

    //------------------------- fit tracks based on maximum energy clusters

  if (stop==3) return;

    //// XZ telescope planes
  
    // upstream telescope planes
    bool missing_x_hit_3a = false;
    init_data_points(3);
    if (stkcluster_maxe[0 ]==-999.9 || stkcluster_maxehitposition[0 ]==-999.9) missing_x_hit_3a = true; else {pos_xy.at(0) = stkcluster_maxehitposition[0 ]; pos_z.at(0) = stk_clusters_z[0]; pos_xy_error.at(0) = 2.0;}
    if (stkcluster_maxe[2 ]==-999.9 || stkcluster_maxehitposition[2 ]==-999.9) missing_x_hit_3a = true; else {pos_xy.at(1) = stkcluster_maxehitposition[2 ]; pos_z.at(1) = stk_clusters_z[2]; pos_xy_error.at(1) = 2.0;}
    if (stkcluster_maxe[4 ]==-999.9 || stkcluster_maxehitposition[4 ]==-999.9) missing_x_hit_3a = true; else {pos_xy.at(2) = stkcluster_maxehitposition[4 ]; pos_z.at(2) = stk_clusters_z[4]; pos_xy_error.at(2) = 2.0;}
    vector<double> x_track_fit_results_3a;
    if (!missing_x_hit_3a) {
      x_track_fit_results_3a = fit_track();
      chi2_x_track_3a->Fill(log10(x_track_fit_results_3a.at(4)));}
    
    // downstream telescope planes
    bool missing_x_hit_3b = false;
    init_data_points(3);
    if (stkcluster_maxe[14]==-999.9 || stkcluster_maxehitposition[14]==-999.9) missing_x_hit_3b = true; else {pos_xy.at(0) = stkcluster_maxehitposition[14]; pos_z.at(0) = stk_clusters_z[14]; pos_xy_error.at(0) = 2.0;}
    if (stkcluster_maxe[16]==-999.9 || stkcluster_maxehitposition[16]==-999.9) missing_x_hit_3b = true; else {pos_xy.at(1) = stkcluster_maxehitposition[16]; pos_z.at(1) = stk_clusters_z[16]; pos_xy_error.at(1) = 2.0;}
    if (stkcluster_maxe[18]==-999.9 || stkcluster_maxehitposition[18]==-999.9) missing_x_hit_3b = true; else {pos_xy.at(2) = stkcluster_maxehitposition[18]; pos_z.at(2) = stk_clusters_z[18]; pos_xy_error.at(2) = 2.0;}
    vector<double> x_track_fit_results_3b;
    if (!missing_x_hit_3b) {
      x_track_fit_results_3b = fit_track();
      chi2_x_track_3b->Fill(log10(x_track_fit_results_3b.at(4)));}
    
    // all telescope planes
    bool missing_x_hit_6 = false;
    init_data_points(6);
    if (stkcluster_maxe[0 ]==-999.9 || stkcluster_maxehitposition[0 ]==-999.9) missing_x_hit_6 = true; else {pos_xy.at(0) = stkcluster_maxehitposition[0 ]; pos_z.at(0) = stk_clusters_z[0 ]; pos_xy_error.at(0) = 2.0;}
    if (stkcluster_maxe[2 ]==-999.9 || stkcluster_maxehitposition[2 ]==-999.9) missing_x_hit_6 = true; else {pos_xy.at(1) = stkcluster_maxehitposition[2 ]; pos_z.at(1) = stk_clusters_z[2 ]; pos_xy_error.at(1) = 2.0;}
    if (stkcluster_maxe[4 ]==-999.9 || stkcluster_maxehitposition[4 ]==-999.9) missing_x_hit_6 = true; else {pos_xy.at(2) = stkcluster_maxehitposition[4 ]; pos_z.at(2) = stk_clusters_z[4 ]; pos_xy_error.at(2) = 2.0;}
    if (stkcluster_maxe[14]==-999.9 || stkcluster_maxehitposition[14]==-999.9) missing_x_hit_6 = true; else {pos_xy.at(3) = stkcluster_maxehitposition[14]; pos_z.at(3) = stk_clusters_z[14]; pos_xy_error.at(3) = 2.0;}
    if (stkcluster_maxe[16]==-999.9 || stkcluster_maxehitposition[16]==-999.9) missing_x_hit_6 = true; else {pos_xy.at(4) = stkcluster_maxehitposition[16]; pos_z.at(4) = stk_clusters_z[16]; pos_xy_error.at(4) = 2.0;}
    if (stkcluster_maxe[18]==-999.9 || stkcluster_maxehitposition[18]==-999.9) missing_x_hit_6 = true; else {pos_xy.at(5) = stkcluster_maxehitposition[18]; pos_z.at(5) = stk_clusters_z[18]; pos_xy_error.at(5) = 2.0;}
    vector<double> x_track_fit_results_6;
    if (!missing_x_hit_6) {
      x_track_fit_results_6 = fit_track();
      chi2_x_track_6->Fill(log10(x_track_fit_results_6.at(4)));}
    // correlations
    if (!missing_x_hit_3a && !missing_x_hit_3b) chi2_x_track_3a_vs_3b->Fill(log10(x_track_fit_results_3a.at(4)),log10(x_track_fit_results_3b.at(4)));
    if (!missing_x_hit_3a && !missing_x_hit_6)  chi2_x_track_3a_vs_6 ->Fill(log10(x_track_fit_results_3a.at(4)),log10(x_track_fit_results_6 .at(4)));
    if (!missing_x_hit_3b && !missing_x_hit_6)  chi2_x_track_3b_vs_6 ->Fill(log10(x_track_fit_results_3b.at(4)),log10(x_track_fit_results_6 .at(4)));



    
    //// YZ telescope planes
    
    // upstream telescope planes
    bool missing_y_hit_3a = false;
    init_data_points(3);
    if (stkcluster_maxe[1 ]==-999.9 || stkcluster_maxehitposition[1 ]==-999.9) missing_y_hit_3a = true; else {pos_xy.at(0) = stkcluster_maxehitposition[1 ]; pos_z.at(0) = stk_clusters_z[1]; pos_xy_error.at(0) = 2.0;}
    if (stkcluster_maxe[3 ]==-999.9 || stkcluster_maxehitposition[3 ]==-999.9) missing_y_hit_3a = true; else {pos_xy.at(1) = stkcluster_maxehitposition[3 ]; pos_z.at(1) = stk_clusters_z[3]; pos_xy_error.at(1) = 2.0;}
    if (stkcluster_maxe[5 ]==-999.9 || stkcluster_maxehitposition[5 ]==-999.9) missing_y_hit_3a = true; else {pos_xy.at(2) = stkcluster_maxehitposition[5 ]; pos_z.at(2) = stk_clusters_z[5]; pos_xy_error.at(2) = 2.0;}
    vector<double> y_track_fit_results_3a;
    if (!missing_y_hit_3a) {
      y_track_fit_results_3a = fit_track();
      chi2_y_track_3a->Fill(log10(y_track_fit_results_3a.at(4)));}
    
    // downstream telescope planes
    bool missing_y_hit_3b = false;
    init_data_points(3);
    if (stkcluster_maxe[15]==-999.9 || stkcluster_maxehitposition[15]==-999.9) missing_y_hit_3b = true; else {pos_xy.at(0) = stkcluster_maxehitposition[15]; pos_z.at(0) = stk_clusters_z[15]; pos_xy_error.at(0) = 2.0;}
    if (stkcluster_maxe[17]==-999.9 || stkcluster_maxehitposition[17]==-999.9) missing_y_hit_3b = true; else {pos_xy.at(1) = stkcluster_maxehitposition[17]; pos_z.at(1) = stk_clusters_z[17]; pos_xy_error.at(1) = 2.0;}
    if (stkcluster_maxe[19]==-999.9 || stkcluster_maxehitposition[19]==-999.9) missing_y_hit_3b = true; else {pos_xy.at(2) = stkcluster_maxehitposition[19]; pos_z.at(2) = stk_clusters_z[19]; pos_xy_error.at(2) = 2.0;}
    vector<double> y_track_fit_results_3b;
    if (!missing_y_hit_3b) {
      y_track_fit_results_3b = fit_track();
      chi2_y_track_3b->Fill(log10(y_track_fit_results_3b.at(4)));}
    
    // all telescope planes
    bool missing_y_hit_6 = false;
    init_data_points(6);
    if (stkcluster_maxe[1 ]==-999.9 || stkcluster_maxehitposition[1 ]==-999.9) missing_y_hit_6 = true; else {pos_xy.at(0) = stkcluster_maxehitposition[1 ]; pos_z.at(0) = stk_clusters_z[1 ]; pos_xy_error.at(0) = 2.0;}
    if (stkcluster_maxe[3 ]==-999.9 || stkcluster_maxehitposition[3 ]==-999.9) missing_y_hit_6 = true; else {pos_xy.at(1) = stkcluster_maxehitposition[3 ]; pos_z.at(1) = stk_clusters_z[3 ]; pos_xy_error.at(1) = 2.0;}
    if (stkcluster_maxe[5 ]==-999.9 || stkcluster_maxehitposition[5 ]==-999.9) missing_y_hit_6 = true; else {pos_xy.at(2) = stkcluster_maxehitposition[5 ]; pos_z.at(2) = stk_clusters_z[5 ]; pos_xy_error.at(2) = 2.0;}
    if (stkcluster_maxe[15]==-999.9 || stkcluster_maxehitposition[15]==-999.9) missing_y_hit_6 = true; else {pos_xy.at(3) = stkcluster_maxehitposition[15]; pos_z.at(3) = stk_clusters_z[15]; pos_xy_error.at(3) = 2.0;}
    if (stkcluster_maxe[17]==-999.9 || stkcluster_maxehitposition[17]==-999.9) missing_y_hit_6 = true; else {pos_xy.at(4) = stkcluster_maxehitposition[17]; pos_z.at(4) = stk_clusters_z[17]; pos_xy_error.at(4) = 2.0;}
    if (stkcluster_maxe[19]==-999.9 || stkcluster_maxehitposition[19]==-999.9) missing_y_hit_6 = true; else {pos_xy.at(5) = stkcluster_maxehitposition[19]; pos_z.at(5) = stk_clusters_z[19]; pos_xy_error.at(5) = 2.0;}
    vector<double> y_track_fit_results_6;
    if (!missing_y_hit_6) {
      y_track_fit_results_6 = fit_track();
      chi2_y_track_6->Fill(log10(y_track_fit_results_6.at(4)));
    }

    
    // correlations
    if (!missing_y_hit_3a && !missing_y_hit_3b) chi2_y_track_3a_vs_3b->Fill(log10(y_track_fit_results_3a.at(4)),log10(y_track_fit_results_3b.at(4)));
    if (!missing_y_hit_3a && !missing_y_hit_6)  chi2_y_track_3a_vs_6 ->Fill(log10(y_track_fit_results_3a.at(4)),log10(y_track_fit_results_6 .at(4)));
    if (!missing_y_hit_3b && !missing_y_hit_6)  chi2_y_track_3b_vs_6 ->Fill(log10(y_track_fit_results_3b.at(4)),log10(y_track_fit_results_6 .at(4)));


  if (stop==4) return;

    //-------------------------- track based event selection ----------------//

   
    double max_chi2 = -1.5;
    bool x_track_3ap_found = false;
    if (!missing_x_hit_3a && log10(x_track_fit_results_3a.at(4)) < max_chi2) x_track_3ap_found = true;
    bool x_track_3bp_found = false;
    if (!missing_x_hit_3b && log10(x_track_fit_results_3b.at(4)) < max_chi2) x_track_3bp_found = true;
    bool x_track_6p_found = false;
    if (!missing_x_hit_6 && log10(x_track_fit_results_6.at(4))   < max_chi2+0.5) x_track_6p_found  = true;
    bool y_track_3ap_found = false;
    if (!missing_y_hit_3a && log10(y_track_fit_results_3a.at(4)) < max_chi2) y_track_3ap_found = true;
    bool y_track_3bp_found = false;
    if (!missing_y_hit_3b && log10(y_track_fit_results_3b.at(4)) < max_chi2) y_track_3bp_found = true;
    bool y_track_6p_found = false;
    if (!missing_y_hit_6 && log10(y_track_fit_results_6.at(4))   < max_chi2+0.5) y_track_6p_found  = true;

    bool x_upstream_track_is_clean = false;
    bool y_upstream_track_is_clean = false;
    //int max_n_cluster = 10;
    double max_logeratio = -0.5;
    if (
    	log10(secondary_clusters_energy_ratio[0])<max_logeratio+0.2 &&
    	log10(secondary_clusters_energy_ratio[2])<max_logeratio &&
    	log10(secondary_clusters_energy_ratio[4])<max_logeratio
    	) x_upstream_track_is_clean = true;
    if (
    	log10(secondary_clusters_energy_ratio[1])<max_logeratio+0.2 &&
    	log10(secondary_clusters_energy_ratio[3])<max_logeratio &&
    	log10(secondary_clusters_energy_ratio[5])<max_logeratio
    	) y_upstream_track_is_clean = true;

    bool good_event = false;
    bool fragmentation = false;
    bool good_upstream_event = false;

    if (x_track_3ap_found &&
    	x_track_3bp_found &&
    	x_track_6p_found &&
    	y_track_3ap_found &&
    	y_track_3bp_found &&
    	y_track_6p_found &&
    	x_upstream_track_is_clean &&
    	y_upstream_track_is_clean
    	) {
      n_good_events++;
      good_event = true; }

    if (x_track_3ap_found &&
    	y_track_3ap_found &&
    	(!x_track_3bp_found ||
    	 !y_track_3bp_found) &&
    	x_upstream_track_is_clean &&
    	y_upstream_track_is_clean
    	) {
      n_fragmentation_events++;      
      fragmentation = true; }
  
    if (x_track_3ap_found &&
    	y_track_3ap_found &&
    	x_upstream_track_is_clean &&
    	y_upstream_track_is_clean
    	) {
      n_good_upsteam_events++;
      good_upstream_event = true; }


    //--------- fill histograms of cluster residuals
    if (good_event) {
      for (int i=0; i<20; i++) {
      	double slope;
    	double intersec;
    	for (int j=0; j<6; j++) {
    	  if (i == x_plane_id[j] ) {
    	    slope = x_track_fit_results_6.at(0);
    	    intersec = x_track_fit_results_6.at(2); }}
    	for (int j=0; j<14; j++) {
    	  if (i == y_plane_id[j] ) {
    	    slope = y_track_fit_results_6.at(0);
    	    intersec = y_track_fit_results_6.at(2); }}
    	  track_ladder_pos[i] = slope*stk_clusters_z[i] + intersec; 
    	cluster_track_residuals_plane[i]->Fill(stkcluster_maxehitposition[i]-track_ladder_pos[i]);
      }
    }

    //------------------- fill histogram of delta phi between upstream and downstream tracks
    if (good_event) {
      double phi_x_3a = atan(x_track_fit_results_3a.at(0))*180.0/3.14;
      double phi_x_3b = atan(x_track_fit_results_3b.at(0))*180.0/3.14;
      delta_phi_upstream_downstream_x_tracks->Fill(phi_x_3a-phi_x_3b);
      double phi_y_3a = atan(y_track_fit_results_3a.at(0))*180.0/3.14;
      double phi_y_3b = atan(y_track_fit_results_3b.at(0))*180.0/3.14;
      delta_phi_upstream_downstream_y_tracks->Fill(phi_y_3a-phi_y_3b); } 

    //------------------------------- cuts based on track fitting
    //if (good_upstream_event) {
    if (good_event) {
    //if (fragmentation) {
    //{
      for (int i=0; i<20; i++){ cluster_q_track[i]->Fill(sqrt(stkcluster_maxe[i])); }
      for (int i=0; i<20; i++){ cluster_adc_track[i]->Fill(stkcluster_maxe[i]); }

      double stkcluster_mean_Z = 0;
      int n_good_clusters = 0;
      for (int i=6; i<14; i++) { stkcluster_mean_Z += stkcluster_maxe[i]; n_good_clusters++; } 
      if (n_good_clusters>0) stkcluster_mean_Z = sqrt(stkcluster_mean_Z/n_good_clusters/60.0);
      cluster_z_dut_track->Fill(stkcluster_mean_Z);

      //--------- fill cluster multiplicities histograms
      //if (x_track_6p_found && y_track_6p_found) 
      for (int i=0; i<20; i++){ cluster_mult_plane[i]->Fill(stknclusters[i]); }
    
      //if (x_track_6p_found && y_track_6p_found)
      for (int i=0; i<20; i++){ maxq_vs_cluster_mult_plane[i]->Fill(sqrt(stkcluster_maxe[i]),stknclusters[i]); }

      //--------- fill secondary clusters histograms
      for (int i=0; i<20; i++){
    	secondary_clusters_energy_ratio_plane[i]->Fill(log10(secondary_clusters_energy_ratio[i])); }
      
      //-------------------------- draw event display -------------------------//
      if (entry>=firsteventtodisplay && event_display_i<neventstodisplay) {
	
      	printf("\n@@@@@@@@@@@@@@@ Event %d [%d] @@@@@@@@@@@@@@@\n",event_display_i,entry);

      	for (int i=0;  i<6; i++) { 
      	  for (uint j=0; j<stk_clusters_pos[x_plane_id[i]].size(); j++) {
      	    stk_hits_XZ[event_display_i]->Fill( stk_clusters_z[x_plane_id[i]], stk_clusters_pos[x_plane_id[i]].at(j), stk_clusters_e[x_plane_id[i]].at(j)); }}

      	for (int i=0;  i<14; i++) { 
      	  for (uint j=0; j<stk_clusters_pos[y_plane_id[i]].size(); j++) {
      	    stk_hits_YZ[event_display_i]->Fill( stk_clusters_z[y_plane_id[i]], stk_clusters_pos[y_plane_id[i]].at(j), stk_clusters_e[y_plane_id[i]].at(j)); }}

      	evtdisp_canv[event_display_i]->cd(1);
      	evtdisp_canv[event_display_i]->cd(1)->SetLogz();
      	stk_hits_XZ[event_display_i]->Draw("COLZ");
      	if (!missing_x_hit_3a) {
      	  x_track_3a_line[event_display_i]->SetX1(stk_clusters_z[0]);
      	  x_track_3a_line[event_display_i]->SetY1(x_track_fit_results_3a.at(0)*stk_clusters_z[0]+x_track_fit_results_3a.at(2));
      	  x_track_3a_line[event_display_i]->SetX2(stk_clusters_z[4]);
      	  x_track_3a_line[event_display_i]->SetY2(x_track_fit_results_3a.at(0)*stk_clusters_z[4]+x_track_fit_results_3a.at(2));
      	  x_track_3a_line[event_display_i]->SetLineColor(kGreen+2);
      	  x_track_3a_line[event_display_i]->SetLineWidth(4);
      	  x_track_3a_line[event_display_i]->SetLineStyle(2);
      	  x_track_3a_line[event_display_i]->Draw("SAME"); }
      	if (!missing_x_hit_3b) {
      	  x_track_3b_line[event_display_i]->SetX1(stk_clusters_z[14]);
      	  x_track_3b_line[event_display_i]->SetY1(x_track_fit_results_3b.at(0)*stk_clusters_z[14]+x_track_fit_results_3b.at(2));
      	  x_track_3b_line[event_display_i]->SetX2(stk_clusters_z[18]);
      	  x_track_3b_line[event_display_i]->SetY2(x_track_fit_results_3b.at(0)*stk_clusters_z[18]+x_track_fit_results_3b.at(2));
      	  x_track_3b_line[event_display_i]->SetLineColor(kGreen+2);
      	  x_track_3b_line[event_display_i]->SetLineWidth(4);
      	  x_track_3b_line[event_display_i]->SetLineStyle(2);
      	  x_track_3b_line[event_display_i]->Draw("SAME"); }
      	if (!missing_x_hit_6) {
      	  x_track_6_line[event_display_i]->SetX1(stk_clusters_z[0]);
      	  x_track_6_line[event_display_i]->SetY1(x_track_fit_results_6.at(0)*stk_clusters_z[0]+x_track_fit_results_6.at(2));
      	  x_track_6_line[event_display_i]->SetX2(stk_clusters_z[18]);
      	  x_track_6_line[event_display_i]->SetY2(x_track_fit_results_6.at(0)*stk_clusters_z[18]+x_track_fit_results_6.at(2));
      	  x_track_6_line[event_display_i]->SetLineColor(kRed);
      	  x_track_6_line[event_display_i]->SetLineWidth(1);
      	  x_track_6_line[event_display_i]->SetLineStyle(1);
      	  x_track_6_line[event_display_i]->Draw("SAME"); }
      	evtdisp_canv[event_display_i]->cd(2);
      	evtdisp_canv[event_display_i]->cd(2)->SetLogz();
      	stk_hits_YZ[event_display_i]->Draw("COLZ");
      	if (!missing_y_hit_3a) {
      	  y_track_3a_line[event_display_i]->SetX1(stk_clusters_z[1]);
      	  y_track_3a_line[event_display_i]->SetY1(y_track_fit_results_3a.at(0)*stk_clusters_z[1]+y_track_fit_results_3a.at(2));
      	  y_track_3a_line[event_display_i]->SetX2(stk_clusters_z[5]);
      	  y_track_3a_line[event_display_i]->SetY2(y_track_fit_results_3a.at(0)*stk_clusters_z[5]+y_track_fit_results_3a.at(2));
      	  y_track_3a_line[event_display_i]->SetLineColor(kGreen+2);
      	  y_track_3a_line[event_display_i]->SetLineWidth(4);
      	  y_track_3a_line[event_display_i]->SetLineStyle(2);
      	  y_track_3a_line[event_display_i]->Draw("SAME"); }
      	if (!missing_y_hit_3b) {
      	  y_track_3b_line[event_display_i]->SetX1(stk_clusters_z[15]);
      	  y_track_3b_line[event_display_i]->SetY1(y_track_fit_results_3b.at(0)*stk_clusters_z[15]+y_track_fit_results_3b.at(2));
      	  y_track_3b_line[event_display_i]->SetX2(stk_clusters_z[19]);
      	  y_track_3b_line[event_display_i]->SetY2(y_track_fit_results_3b.at(0)*stk_clusters_z[19]+y_track_fit_results_3b.at(2));
      	  y_track_3b_line[event_display_i]->SetLineColor(kGreen+2);
      	  y_track_3b_line[event_display_i]->SetLineWidth(4);
      	  y_track_3b_line[event_display_i]->SetLineStyle(2);
      	  y_track_3b_line[event_display_i]->Draw("SAME"); }
      	if (!missing_y_hit_6) {
      	  y_track_6_line[event_display_i]->SetX1(stk_clusters_z[1]);
      	  y_track_6_line[event_display_i]->SetY1(y_track_fit_results_6.at(0)*stk_clusters_z[1]+y_track_fit_results_6.at(2));
      	  y_track_6_line[event_display_i]->SetX2(stk_clusters_z[19]);
      	  y_track_6_line[event_display_i]->SetY2(y_track_fit_results_6.at(0)*stk_clusters_z[19]+y_track_fit_results_6.at(2));
      	  y_track_6_line[event_display_i]->SetLineColor(kRed);
      	  y_track_6_line[event_display_i]->SetLineWidth(1);
      	  y_track_6_line[event_display_i]->SetLineStyle(1);
      	  y_track_6_line[event_display_i]->Draw("SAME"); }
      	evtdisp_canv[event_display_i]->SaveAs(TString::Format("plots/run%d/evt_dsp_entry%d.png",(int)ion_mean_pos[run_id*16],entry));
      	//------------------------------------------------------------------------//
      	event_display_i++; }

      //------------------------- scintillators ----------------------------//
      for(int i=0; i<scintillators->GetLast()+1; i++){
        Scintillator* scintillator = (Scintillator*)scintillators->ConstructedAt(i);
        int nsamples = scintillator->GetSize();
        // printf("scint%d: %d samples\n",i,nsamples);
        for(int j=0; j<nsamples; j++){
      	scintillator->GetValue(j);
        }
      }
      // ... there should be 4 scintillators per event, corresponding to wave0, wave1, wave2, wave3
      //  scintillator class is defined here: Scintillator.h

      // Loop over reconstructed SCINTILLATOR Clusters
      for(int i=0; i<scintclusters->GetLast()+1; i++){
    	ScintillatorCluster* cluster = (ScintillatorCluster*)scintclusters->ConstructedAt(i);
    	int scintid = cluster->GetScintillatorID();                   // 0 - wave0,  1 - wave1,  2 - wave2,  3 - wave3
    	int scintclustersize = cluster->GetSize();                    // Number of samples in the cluster
    	int scintclusterenergy = cluster->GetSignalTotal();           // Total energy of the cluster
    	int scintclusterfirstsampleid = cluster->GetPosition(0);

    	// //printf("scint.cluster %d: wave%d size=%d energy=%d firstsampleid=%d\n",i,scintid,scintclustersize,scintclusterenergy,scintclusterfirstsampleid);
    	// for(int j=0; j<scintclustersize; j++){
    	//   //cluster->GetPosition(j);                     // position of j-th sample in the cluster; 
    	//   //   there are 4000 samples per event (may be more, depending on the run)
    	//   //   so, for example, given the cluster starts at 100-th sample and lasts for 20 
    	//   //   samples. the positions will be 100,101, 102,...119
	
    	//   // signal of j-th sample in the cluster (pedestal is subtracted)
    	//   if (cluster->GetSignal(j) > 3954) suturated_signal[scintid] = true;
    	// }


    	scint_cluster_q[scintid]->Fill(sqrt(scintclusterenergy));
    	scint_cluster_width[scintid]->Fill(scintclustersize);
    	stkmaxq_vs_scint_cluster_q[scintid]->Fill(sqrt(stkcluster_maxe[0]),sqrt(scintclusterenergy));

      }// loop over scintillator clusters

    }// cuts based on track fitting
  
  }
  // END OF EVENT LOOP
  
  //----------------------------------- correlations
  TCanvas *aligned_tlp_prof_canv;
  TCanvas *aligned_dut_prof_canv;
  if (neventstoprocess >= neventsforalignment) {
    //-------- correlations
    for (int i=0; i<18; i++) {
      aligned_corr[i]->Sumw2();
      aligned_corr[i]->Scale(1.0/aligned_corr[i]->Integral()); }
    //-------- remove background hits
    for (int i=1; i<=200; i++) {
      for (int j=1; j<=200; j++) {
  	for (int k=0; k<18; k++) {
  	  if (aligned_corr[k]->GetBinContent(i,j)<0.00015) aligned_corr[k]->SetBinContent(i,j,0); }}}
    //-------- define correlation profiles
    TProfile *aligned_corr_prof[18];
    //-------- plot correlatoins in telelscope planes
    aligned_tlp_prof_canv = new TCanvas("aligned_tlp_prof_canv","aligned_tlp_prof_canv",1800,1000);
    aligned_tlp_prof_canv->Divide(5,2);
    for (int k=0; k<10; k++) {
      aligned_tlp_prof_canv->cd(k+1);  
      aligned_corr_prof[k] = aligned_corr[k]->ProfileX();    
      aligned_corr_prof[k]->Fit("pol1","Q","",35,50); 
      aligned_corr_prof[k]->GetYaxis()->SetRangeUser(0,100);     
      aligned_tlp_prof_canv->cd(k+1)->SetGridx();  
      aligned_tlp_prof_canv->cd(k+1)->SetGridy(); }
    //-------- plot correlations in DUT planes
    aligned_dut_prof_canv = new TCanvas("aligned_dut_prof_canv","aligned_dut_prof_canv",1800,1000);
    aligned_dut_prof_canv->Divide(4,2);
    for (int k=10; k<18; k++) {
      aligned_dut_prof_canv->cd(k-9);  
      aligned_corr_prof[k] = aligned_corr[k]->ProfileX();    
      aligned_corr_prof[k]->Fit("pol1","Q","",30,60); 
      aligned_corr_prof[k]->GetYaxis()->SetRangeUser(0,100);     
      aligned_dut_prof_canv->cd(k-9)->SetGridx();  
      aligned_dut_prof_canv->cd(k-9)->SetGridy(); }}
  aligned_tlp_prof_canv->SaveAs(TString::Format("plots/run%d/aligned_tlp_prof",(int)ion_mean_pos[run_id*16])+"_canv.png");
  aligned_dut_prof_canv->SaveAs(TString::Format("plots/run%d/aligned_dut_prof",(int)ion_mean_pos[run_id*16])+"_canv.png");

  //----------------------------------- plot chi2 of tracks
  TCanvas *chi2_tracks_canv = new TCanvas("chi2_tracks_canv","chi2_tracks_canv",1800,1000);
  chi2_tracks_canv->Divide(3,2);
  chi2_tracks_canv->cd(1); chi2_x_track_3a->Draw();
  chi2_tracks_canv->cd(2); chi2_x_track_3b->Draw();
  chi2_tracks_canv->cd(3); chi2_x_track_6 ->Draw();
  chi2_tracks_canv->cd(4); chi2_y_track_3a->Draw();
  chi2_tracks_canv->cd(5); chi2_y_track_3b->Draw();
  chi2_tracks_canv->cd(6); chi2_y_track_6 ->Draw();
  chi2_tracks_canv->SaveAs(TString::Format("plots/run%d/chi2_tracks",(int)ion_mean_pos[run_id*16])+"_canv.png");

  TCanvas *chi2_corr_tracks_canv = new TCanvas("chi2_corr_tracks_canv","chi2_corr_tracks_canv",1800,1000);
  chi2_corr_tracks_canv->Divide(3,2);
  chi2_corr_tracks_canv->cd(1); chi2_x_track_3a_vs_3b->Draw("COLZ");
  chi2_corr_tracks_canv->cd(2); chi2_x_track_3a_vs_6 ->Draw("COLZ");
  chi2_corr_tracks_canv->cd(3); chi2_x_track_3b_vs_6 ->Draw("COLZ");
  chi2_corr_tracks_canv->cd(4); chi2_y_track_3a_vs_3b->Draw("COLZ"); 
  chi2_corr_tracks_canv->cd(5); chi2_y_track_3a_vs_6 ->Draw("COLZ"); 
  chi2_corr_tracks_canv->cd(6); chi2_y_track_3b_vs_6 ->Draw("COLZ"); 
  chi2_corr_tracks_canv->SaveAs(TString::Format("plots/run%d/chi2_corr_tracks",(int)ion_mean_pos[run_id*16])+"_canv.png");

  //---------------------------------- plot cluster charges
  TCanvas *cluster_q_x_telescope_plane_canv = new TCanvas("cluster_q_x_telescope_plane_canv","cluster_q_x_telescope_plane_canv",1800,1000);
  cluster_q_x_telescope_plane_canv->Divide(3,2);
  for (int i=0; i<6; i++) {
    cluster_q_x_telescope_plane_canv->cd(i+1); /*cluster_q_track[x_plane_id[i]]->Rebin(5);*/ cluster_q_track[x_plane_id[i]]->Draw(); }
  cluster_q_x_telescope_plane_canv->SaveAs(TString::Format("plots/run%d/cluster_q_x_telescope_plane",(int)ion_mean_pos[run_id*16])+"_canv.png");

  TCanvas *cluster_q_y_telescope_plane_canv = new TCanvas("cluster_q_y_telescope_plane_canv","cluster_q_y_telescope_plane_canv",1800,1000);
  cluster_q_y_telescope_plane_canv->Divide(3,2);
  for (int i=0; i<6; i++) {
    cluster_q_y_telescope_plane_canv->cd(i+1); /*cluster_q_track[y_plane_id[i]]->Rebin(5);*/ cluster_q_track[y_plane_id[i]]->Draw(); }
  cluster_q_y_telescope_plane_canv->SaveAs(TString::Format("plots/run%d/cluster_q_y_telescope_plane",(int)ion_mean_pos[run_id*16])+"_canv.png");

  TCanvas *cluster_q_y_dut_plane_canv = new TCanvas("cluster_q_y_dut_plane_canv","cluster_q_y_dut_plane_canv",1800,1000);
  cluster_q_y_dut_plane_canv->Divide(4,2);
  for (int i=6; i<14; i++) { cluster_q_y_dut_plane_canv->cd(i-5); /*cluster_q_track[i]->Rebin(5);*/ cluster_q_track[i]->Draw();}
  cluster_q_y_dut_plane_canv->SaveAs(TString::Format("plots/run%d/cluster_q_y_dut_plane",(int)ion_mean_pos[run_id*16])+"_canv.png");

  //-----

  TCanvas *cluster_z_dut_canv = new TCanvas("cluster_z_dut_canv","cluster_z_dut_canv",800,600);
  cluster_z_dut_canv->Divide(4,2);
  cluster_z_dut_track->Draw();
  cluster_z_dut_canv->SaveAs(TString::Format("plots/run%d/cluster_z_dut_",(int)ion_mean_pos[run_id*16])+"_canv.png");
  cluster_z_dut_canv->SaveAs(TString::Format("plots/run%d/cluster_z_dut_",(int)ion_mean_pos[run_id*16])+"_canv.C");
  
  //-----

  TCanvas *cluster_adc_x_telescope_plane_canv = new TCanvas("cluster_adc_x_telescope_plane_canv","cluster_adc_x_telescope_plane_canv",1800,1000);
  cluster_adc_x_telescope_plane_canv->Divide(3,2);
  for (int i=0; i<6; i++) {
    cluster_adc_x_telescope_plane_canv->cd(i+1); /*cluster_adc_track[x_plane_id[i]]->Rebin(5);*/ cluster_adc_track[x_plane_id[i]]->Draw(); }
  cluster_adc_x_telescope_plane_canv->SaveAs(TString::Format("plots/run%d/cluster_adc_x_telescope_plane",(int)ion_mean_pos[run_id*16])+"_canv.png");

  TCanvas *cluster_adc_y_telescope_plane_canv = new TCanvas("cluster_adc_y_telescope_plane_canv","cluster_adc_y_telescope_plane_canv",1800,1000);
  cluster_adc_y_telescope_plane_canv->Divide(3,2);
  for (int i=0; i<6; i++) {
    cluster_adc_y_telescope_plane_canv->cd(i+1); /*cluster_adc_track[y_plane_id[i]]->Rebin(5);*/ cluster_adc_track[y_plane_id[i]]->Draw(); }
  cluster_adc_y_telescope_plane_canv->SaveAs(TString::Format("plots/run%d/cluster_adc_y_telescope_plane",(int)ion_mean_pos[run_id*16])+"_canv.png");

  TCanvas *cluster_adc_y_dut_plane_canv = new TCanvas("cluster_adc_y_dut_plane_canv","cluster_adc_y_dut_plane_canv",1800,1000);
  cluster_adc_y_dut_plane_canv->Divide(4,2);
  for (int i=6; i<14; i++) { cluster_adc_y_dut_plane_canv->cd(i-5); /*cluster_adc_track[i]->Rebin(5);*/ cluster_adc_track[i]->Draw();}
  cluster_adc_y_dut_plane_canv->SaveAs(TString::Format("plots/run%d/cluster_adc_y_dut_plane",(int)ion_mean_pos[run_id*16])+"_canv.png");

  //---------------------------------- plot occupancies

  TCanvas *occupancy_x_telescope_plane_canv = new TCanvas("occupancy_x_telescope_plane_canv","occupancy_x_telescope_plane_canv",1800,1000);
  occupancy_x_telescope_plane_canv->Divide(3,2);
  for (int i=0; i<6; i++) {
    occupancy_x_telescope_plane_canv->cd(i+1); /*occupancy_plane[x_plane_id[i]]->Rebin(5);*/ occupancy_plane[x_plane_id[i]]->Draw(); }
  occupancy_x_telescope_plane_canv->SaveAs(TString::Format("plots/run%d/occupancy_x_telescope_plane",(int)ion_mean_pos[run_id*16])+"_canv.png");

  TCanvas *occupancy_y_telescope_plane_canv = new TCanvas("occupancy_y_telescope_plane_canv","occupancy_y_telescope_plane_canv",1800,1000);
  occupancy_y_telescope_plane_canv->Divide(3,2);
  for (int i=0; i<6; i++) {
    occupancy_y_telescope_plane_canv->cd(i+1); /*occupancy_plane[y_plane_id[i]]->Rebin(5);*/ occupancy_plane[y_plane_id[i]]->Draw(); }
  occupancy_y_telescope_plane_canv->SaveAs(TString::Format("plots/run%d/occupancy_y_telescope_plane",(int)ion_mean_pos[run_id*16])+"_canv.png");

  TCanvas *occupancy_y_dut_plane_canv = new TCanvas("occupancy_y_dut_plane_canv","occupancy_y_dut_plane_canv",1800,1000);
  occupancy_y_dut_plane_canv->Divide(4,2);
  for (int i=6; i<14; i++) { occupancy_y_dut_plane_canv->cd(i-5); /*occupancy_plane[i]->Rebin(5);*/ occupancy_plane[i]->Draw();}
  occupancy_y_dut_plane_canv->SaveAs(TString::Format("plots/run%d/occupancy_y_dut_plane",(int)ion_mean_pos[run_id*16])+"_canv.png");

  TCanvas *energy_occupancy_x_telescope_plane_canv = new TCanvas("energy_occupancy_x_telescope_plane_canv","energy_occupancy_x_telescope_plane_canv",1800,1000);
  energy_occupancy_x_telescope_plane_canv->Divide(3,2);
  for (int i=0; i<6; i++) {
    energy_occupancy_x_telescope_plane_canv->cd(i+1); energy_occupancy_plane[x_plane_id[i]]->Rebin2D(5,1); energy_occupancy_plane[x_plane_id[i]]->Draw("COLZ"); energy_occupancy_x_telescope_plane_canv->cd(i+1)->SetLogz(); }
  energy_occupancy_x_telescope_plane_canv->SaveAs(TString::Format("plots/run%d/energy_occupancy_x_telescope_plane",(int)ion_mean_pos[run_id*16])+"_canv.png");

  TCanvas *energy_occupancy_y_telescope_plane_canv = new TCanvas("energy_occupancy_y_telescope_plane_canv","energy_occupancy_y_telescope_plane_canv",1800,1000);
  energy_occupancy_y_telescope_plane_canv->Divide(3,2);
  for (int i=0; i<6; i++) {
    energy_occupancy_y_telescope_plane_canv->cd(i+1); energy_occupancy_plane[y_plane_id[i]]->Rebin2D(5,1); energy_occupancy_plane[y_plane_id[i]]->Draw("COLZ"); energy_occupancy_y_telescope_plane_canv->cd(i+1)->SetLogz(); }
  energy_occupancy_y_telescope_plane_canv->SaveAs(TString::Format("plots/run%d/energy_occupancy_y_telescope_plane",(int)ion_mean_pos[run_id*16])+"_canv.png");

  TCanvas *energy_occupancy_y_dut_plane_canv = new TCanvas("energy_occupancy_y_dut_plane_canv","energy_occupancy_y_dut_plane_canv",1800,1000);
  energy_occupancy_y_dut_plane_canv->Divide(4,2);
  for (int i=6; i<14; i++) { energy_occupancy_y_dut_plane_canv->cd(i-5); energy_occupancy_plane[i]->Rebin2D(5,1); energy_occupancy_plane[i]->Draw("COLZ"); energy_occupancy_y_dut_plane_canv->cd(i-5)->SetLogz();}
  energy_occupancy_y_dut_plane_canv->SaveAs(TString::Format("plots/run%d/energy_occupancy_y_dut_plane",(int)ion_mean_pos[run_id*16])+"_canv.png");

  TCanvas *nstr_occupancy_x_telescope_plane_canv = new TCanvas("nstr_occupancy_x_telescope_plane_canv","nstr_occupancy_x_telescope_plane_canv",1800,1000);
  nstr_occupancy_x_telescope_plane_canv->Divide(3,2);
  for (int i=0; i<6; i++) {
    nstr_occupancy_x_telescope_plane_canv->cd(i+1); nstr_occupancy_plane[x_plane_id[i]]->Rebin2D(5,1); nstr_occupancy_plane[x_plane_id[i]]->Draw("COLZ"); nstr_occupancy_x_telescope_plane_canv->cd(i+1)->SetLogz(); }
  nstr_occupancy_x_telescope_plane_canv->SaveAs(TString::Format("plots/run%d/nstr_occupancy_x_telescope_plane",(int)ion_mean_pos[run_id*16])+"_canv.png");

  TCanvas *nstr_occupancy_y_telescope_plane_canv = new TCanvas("nstr_occupancy_y_telescope_plane_canv","nstr_occupancy_y_telescope_plane_canv",1800,1000);
  nstr_occupancy_y_telescope_plane_canv->Divide(3,2);
  for (int i=0; i<6; i++) {
    nstr_occupancy_y_telescope_plane_canv->cd(i+1); nstr_occupancy_plane[y_plane_id[i]]->Rebin2D(5,1); nstr_occupancy_plane[y_plane_id[i]]->Draw("COLZ"); nstr_occupancy_y_telescope_plane_canv->cd(i+1)->SetLogz(); }
  nstr_occupancy_y_telescope_plane_canv->SaveAs(TString::Format("plots/run%d/nstr_occupancy_y_telescope_plane",(int)ion_mean_pos[run_id*16])+"_canv.png");

  TCanvas *nstr_occupancy_y_dut_plane_canv = new TCanvas("nstr_occupancy_y_dut_plane_canv","nstr_occupancy_y_dut_plane_canv",1800,1000);
  nstr_occupancy_y_dut_plane_canv->Divide(4,2);
  for (int i=6; i<14; i++) { nstr_occupancy_y_dut_plane_canv->cd(i-5); nstr_occupancy_plane[i]->Rebin2D(5,1); nstr_occupancy_plane[i]->Draw("COLZ"); nstr_occupancy_y_dut_plane_canv->cd(i-5)->SetLogz();}
  nstr_occupancy_y_dut_plane_canv->SaveAs(TString::Format("plots/run%d/nstr_occupancy_y_dut_plane",(int)ion_mean_pos[run_id*16])+"_canv.png");

  //------------------------------------------ plot cluster multiplicities
  TCanvas *cluster_mult_x_telescope_plane_canv = new TCanvas("cluster_mult_x_telescope_plane_canv","cluster_mult_x_telescope_plane_canv",1800,1000);
  cluster_mult_x_telescope_plane_canv->Divide(3,2);
  for (int i=0; i<6; i++) {
    cluster_mult_x_telescope_plane_canv->cd(i+1); cluster_mult_plane[x_plane_id[i]]->Draw(); }
  cluster_mult_x_telescope_plane_canv->SaveAs(TString::Format("plots/run%d/cluster_mult_x_telescope_plane",(int)ion_mean_pos[run_id*16])+"_canv.png");

  TCanvas *cluster_mult_y_telescope_plane_canv = new TCanvas("cluster_mult_y_telescope_plane_canv","cluster_mult_y_telescope_plane_canv",1800,1000);
  cluster_mult_y_telescope_plane_canv->Divide(3,2);
  for (int i=0; i<6; i++) {
    cluster_mult_y_telescope_plane_canv->cd(i+1); cluster_mult_plane[y_plane_id[i]]->Draw(); }
  cluster_mult_y_telescope_plane_canv->SaveAs(TString::Format("plots/run%d/cluster_mult_y_telescope_plane",(int)ion_mean_pos[run_id*16])+"_canv.png");

  TCanvas *cluster_mult_y_dut_plane_canv = new TCanvas("cluster_mult_y_dut_plane_canv","cluster_mult_y_dut_plane_canv",1800,1000);
  cluster_mult_y_dut_plane_canv->Divide(4,2);
  for (int i=6; i<14; i++) { cluster_mult_y_dut_plane_canv->cd(i-5); cluster_mult_plane[i]->Draw();}
  cluster_mult_y_dut_plane_canv->SaveAs(TString::Format("plots/run%d/cluster_mult_y_dut_plane",(int)ion_mean_pos[run_id*16])+"_canv.png");

  TCanvas *maxq_vs_cluster_mult_x_telescope_plane_canv = new TCanvas("maxq_vs_cluster_mult_x_telescope_plane_canv","maxq_vs_cluster_mult_x_telescope_plane_canv",1800,1000);
  maxq_vs_cluster_mult_x_telescope_plane_canv->Divide(3,2);
  for (int i=0; i<6; i++) {
    maxq_vs_cluster_mult_x_telescope_plane_canv->cd(i+1); maxq_vs_cluster_mult_plane[x_plane_id[i]]->Draw("COLZ"); maxq_vs_cluster_mult_x_telescope_plane_canv->cd(i+1)->SetLogz(); }
  maxq_vs_cluster_mult_x_telescope_plane_canv->SaveAs(TString::Format("plots/run%d/maxq_vs_cluster_mult_x_telescope_plane",(int)ion_mean_pos[run_id*16])+"_canv.png");

  TCanvas *maxq_vs_cluster_mult_y_telescope_plane_canv = new TCanvas("maxq_vs_cluster_mult_y_telescope_plane_canv","maxq_vs_cluster_mult_y_telescope_plane_canv",1800,1000);
  maxq_vs_cluster_mult_y_telescope_plane_canv->Divide(3,2);
  for (int i=0; i<6; i++) {
    maxq_vs_cluster_mult_y_telescope_plane_canv->cd(i+1); maxq_vs_cluster_mult_plane[y_plane_id[i]]->Draw("COLZ"); maxq_vs_cluster_mult_y_telescope_plane_canv->cd(i+1)->SetLogz(); }
  maxq_vs_cluster_mult_y_telescope_plane_canv->SaveAs(TString::Format("plots/run%d/maxq_vs_cluster_mult_y_telescope_plane",(int)ion_mean_pos[run_id*16])+"_canv.png");

  TCanvas *maxq_vs_cluster_mult_y_dut_plane_canv = new TCanvas("maxq_vs_cluster_mult_y_dut_plane_canv","maxq_vs_cluster_mult_y_dut_plane_canv",1800,1000);
  maxq_vs_cluster_mult_y_dut_plane_canv->Divide(4,2);
  for (int i=6; i<14; i++) { maxq_vs_cluster_mult_y_dut_plane_canv->cd(i-5); maxq_vs_cluster_mult_plane[i]->Draw("COLZ"); maxq_vs_cluster_mult_y_dut_plane_canv->cd(i-5)->SetLogz(); }
  maxq_vs_cluster_mult_y_dut_plane_canv->SaveAs(TString::Format("plots/run%d/maxq_vs_cluster_mult_y_dut_plane",(int)ion_mean_pos[run_id*16])+"_canv.png");

  //----------------------------------- secondary clusters
  TCanvas *secondary_clusters_energy_ratio_x_telescope_plane_canv = new TCanvas("secondary_clusters_energy_ratio_x_telescope_plane_canv","secondary_clusters_energy_ratio_x_telescope_plane_canv",1800,1000);
  secondary_clusters_energy_ratio_x_telescope_plane_canv->Divide(3,2);
  for (int i=0; i<6; i++) {
    secondary_clusters_energy_ratio_x_telescope_plane_canv->cd(i+1); secondary_clusters_energy_ratio_plane[x_plane_id[i]]->Draw(); /*secondary_clusters_energy_ratio_x_telescope_plane_canv->cd(i+1)->SetLogy();*/ }
  secondary_clusters_energy_ratio_x_telescope_plane_canv->SaveAs(TString::Format("plots/run%d/secondary_clusters_energy_ratio_x_telescope_plane",(int)ion_mean_pos[run_id*16])+"_canv.png");

  TCanvas *secondary_clusters_energy_ratio_y_telescope_plane_canv = new TCanvas("secondary_clusters_energy_ratio_y_telescope_plane_canv","secondary_clusters_energy_ratio_y_telescope_plane_canv",1800,1000);
  secondary_clusters_energy_ratio_y_telescope_plane_canv->Divide(3,2);
  for (int i=0; i<6; i++) {
    secondary_clusters_energy_ratio_y_telescope_plane_canv->cd(i+1); secondary_clusters_energy_ratio_plane[y_plane_id[i]]->Draw(); /*secondary_clusters_energy_ratio_y_telescope_plane_canv->cd(i+1)->SetLogy();*/ }
  secondary_clusters_energy_ratio_y_telescope_plane_canv->SaveAs(TString::Format("plots/run%d/secondary_clusters_energy_ratio_y_telescope_plane",(int)ion_mean_pos[run_id*16])+"_canv.png");

  TCanvas *secondary_clusters_energy_ratio_y_dut_plane_canv = new TCanvas("secondary_clusters_energy_ratio_y_dut_plane_canv","secondary_clusters_energy_ratio_y_dut_plane_canv",1800,1000);
  secondary_clusters_energy_ratio_y_dut_plane_canv->Divide(4,2);
  for (int i=6; i<14; i++) { secondary_clusters_energy_ratio_y_dut_plane_canv->cd(i-5); secondary_clusters_energy_ratio_plane[i]->Draw(); /*secondary_clusters_energy_ratio_y_dut_plane_canv->cd(i-5)->SetLogy();*/ }
  secondary_clusters_energy_ratio_y_dut_plane_canv->SaveAs(TString::Format("plots/run%d/secondary_clusters_energy_ratio_y_dut_plane",(int)ion_mean_pos[run_id*16])+"_canv.png");

  //----------------------------------- scintillator clusters
  // TCanvas *scint_cluster_q_canv = new TCanvas("scint_cluster_q_canv","scint_cluster_q_canv",1800,1000);
  // scint_cluster_q_canv->Divide(2,2);
  // for (int i=0; i<4; i++){
  //   scint_cluster_q_canv->cd(i+1);
  //   if (i==0 || i==2) scint_cluster_q[i]->GetXaxis()->SetRangeUser(0,200);
  //   scint_cluster_q[i]->Draw();
    
  //   if (i==2) { // scint. ch2 - non amp., upstream
      
  //     TString multi_gaus_fitf_str = "";
  //     int gaus_first_par_id = 0;
  //     for (int j=1; j<16; j++) {
  // 	if (ion_mean_pos[run_id*16+j]!=-1) {
  // 	  if (multi_gaus_fitf_str!="") multi_gaus_fitf_str += "+";
  // 	  multi_gaus_fitf_str += "gaus("+TString::Format("%d",gaus_first_par_id)+")";
  // 	  gaus_first_par_id += 3; }}
  //     TF1 *multi_gaus_fitf = new TF1("multi_gaus_fitf",multi_gaus_fitf_str);
  //     double upper_limit = -999.99;
  //     double lower_limit = -999.99;
  //     gaus_first_par_id = 0;
  //     for (int j=1; j<16; j++) {
  // 	if (ion_mean_pos[run_id*16+j]!=-1) {
  // 	  // set initial values for parameters
  // 	  multi_gaus_fitf->SetParameter(gaus_first_par_id,1);
  // 	  multi_gaus_fitf->SetParameter(gaus_first_par_id+1,ion_mean_pos[run_id*16+j]-4);
  // 	  multi_gaus_fitf->SetParameter(gaus_first_par_id+2,0.6);
  // 	  // set limits
  // 	  multi_gaus_fitf->SetParLimits(gaus_first_par_id+1,ion_mean_pos[run_id*16+j]-4,ion_mean_pos[run_id*16+j]+1.5);
  // 	  multi_gaus_fitf->SetParLimits(gaus_first_par_id+2,0.6,2.5);
  // 	  if (lower_limit==-999.99) lower_limit = ion_mean_pos[run_id*16+j]-3;
  // 	  upper_limit = ion_mean_pos[run_id*16+j]+1;
  // 	  gaus_first_par_id += 3; }}

  //     multi_gaus_fitf->SetRange(lower_limit,upper_limit);
  //     TVirtualFitter::SetMaxIterations(50000);
  //     scint_cluster_q[i]->Fit(multi_gaus_fitf,"RLM");

  //     printf("\n");
  //     printf("%4s %10s %10s %10s\n","ion","const","mean","sigma");
  //     gaus_first_par_id = 0;
  //     for (int j=1; j<16; j++) {
  // 	if (ion_mean_pos[run_id*16+j]!=-1) {
  // 	  TString ion = ion_name[j];
  // 	  double norm = multi_gaus_fitf->GetParameter(gaus_first_par_id);
  // 	  double mean = multi_gaus_fitf->GetParameter(gaus_first_par_id+1);
  // 	  double sigma = multi_gaus_fitf->GetParameter(gaus_first_par_id+2);
  // 	  gaus_first_par_id += 3;
  // 	  printf("%4s %10.2f %10.2f %10.2f\n",ion.Data(),norm,mean,sigma); }}
  //     printf("\n");
  //   }
  // }
  // scint_cluster_q_canv->SaveAs(TString::Format("plots/run%d/scint_cluster_q",(int)ion_mean_pos[run_id*16])+"_canv.png");

  TCanvas *scint_cluster_pos_canv = new TCanvas("scint_cluster_pos_canv","scint_cluster_pos_canv",1800,1000);
  scint_cluster_pos_canv->Divide(4,2);
  for (int i=0; i<7; i++){
    scint_cluster_pos_canv->cd(i+1);
    scint_cluster_pos_canv->cd(i+1)->SetLogy();
    scint_cluster_pos[i]->Draw(); }  
  scint_cluster_pos_canv->SaveAs(TString::Format("plots/run%d/scint_cluster_pos",(int)ion_mean_pos[run_id*16])+"_canv.png");

  TCanvas *scint_cluster_width_canv = new TCanvas("scint_cluster_width_canv","scint_cluster_width_canv",1800,1000);
  scint_cluster_width_canv->Divide(4,2);
  for (int i=0; i<7; i++){
    scint_cluster_width_canv->cd(i+1);
    scint_cluster_width[i]->Draw(); }  
  scint_cluster_width_canv->SaveAs(TString::Format("plots/run%d/scint_cluster_width",(int)ion_mean_pos[run_id*16])+"_canv.png");

  TCanvas *stkmaxq_vs_scint_cluster_q_canv = new TCanvas("stkmaxq_vs_scint_cluster_q_canv","stkmaxq_vs_scint_cluster_q_canv",1800,1000);
  stkmaxq_vs_scint_cluster_q_canv->Divide(4,2);
  for (int i=0; i<7; i++){
    stkmaxq_vs_scint_cluster_q_canv->cd(i+1);
    stkmaxq_vs_scint_cluster_q_canv->cd(i+1)->SetLogz();
    if (i==0 || i==2) stkmaxq_vs_scint_cluster_q[i]->GetYaxis()->SetRangeUser(0,200);
    stkmaxq_vs_scint_cluster_q[i]->Draw("COLZ"); }
  stkmaxq_vs_scint_cluster_q_canv->SaveAs(TString::Format("plots/run%d/stkmaxq_vs_scint_cluster_q",(int)ion_mean_pos[run_id*16])+"_canv.png");

  //----------------------------------- cluster track residuals
  TCanvas *cluster_track_residuals_x_telescope_plane_canv = new TCanvas("cluster_track_residuals_x_telescope_plane_canv","cluster_track_residuals_x_telescope_plane_canv",1800,1000);
  cluster_track_residuals_x_telescope_plane_canv->Divide(3,2);
  for (int i=0; i<6; i++) {
    cluster_track_residuals_x_telescope_plane_canv->cd(i+1); cluster_track_residuals_plane[x_plane_id[i]]->Draw(); 
    cluster_track_residuals_x_telescope_plane_canv->cd(i+1)->SetLogy(); cluster_track_residuals_plane[x_plane_id[i]]->Fit("gaus","Q"); }
  cluster_track_residuals_x_telescope_plane_canv->SaveAs(TString::Format("plots/run%d/cluster_track_residuals_x_telescope_plane",(int)ion_mean_pos[run_id*16])+"_canv.png");

  TCanvas *cluster_track_residuals_y_telescope_plane_canv = new TCanvas("cluster_track_residuals_y_telescope_plane_canv","cluster_track_residuals_y_telescope_plane_canv",1800,1000);
  cluster_track_residuals_y_telescope_plane_canv->Divide(3,2);
  for (int i=0; i<6; i++) {
    cluster_track_residuals_y_telescope_plane_canv->cd(i+1); cluster_track_residuals_plane[y_plane_id[i]]->Draw(); 
    cluster_track_residuals_y_telescope_plane_canv->cd(i+1)->SetLogy(); cluster_track_residuals_plane[y_plane_id[i]]->Fit("gaus","Q"); }
  cluster_track_residuals_y_telescope_plane_canv->SaveAs(TString::Format("plots/run%d/cluster_track_residuals_y_telescope_plane",(int)ion_mean_pos[run_id*16])+"_canv.png");

  TCanvas *cluster_track_residuals_y_dut_plane_canv = new TCanvas("cluster_track_residuals_y_dut_plane_canv","cluster_track_residuals_y_dut_plane_canv",1800,1000);
  cluster_track_residuals_y_dut_plane_canv->Divide(4,2);
  for (int i=6; i<14; i++) { cluster_track_residuals_y_dut_plane_canv->cd(i-5); cluster_track_residuals_plane[i]->Draw(); 
    cluster_track_residuals_y_dut_plane_canv->cd(i-5)->SetLogy(); cluster_track_residuals_plane[i]->Fit("gaus","Q"); }
  cluster_track_residuals_y_dut_plane_canv->SaveAs(TString::Format("plots/run%d/cluster_track_residuals_y_dut_plane",(int)ion_mean_pos[run_id*16])+"_canv.png");

  //----------------------------------- delta phi between upstream and downstream tracks
  TCanvas *delta_phi_upstream_downstream_canv = new TCanvas("delta_phi_upstream_downstream_canv","delta_phi_upstream_downstream_canv",1800,1000);
  delta_phi_upstream_downstream_canv->Divide(2,1);
  delta_phi_upstream_downstream_canv->cd(1); delta_phi_upstream_downstream_canv->cd(1)->SetLogy(); delta_phi_upstream_downstream_x_tracks->Draw(); delta_phi_upstream_downstream_x_tracks->Fit("gaus","Q");
 delta_phi_upstream_downstream_canv->cd(2); delta_phi_upstream_downstream_canv->cd(2)->SetLogy(); delta_phi_upstream_downstream_y_tracks->Draw(); delta_phi_upstream_downstream_y_tracks->Fit("gaus","Q");
  delta_phi_upstream_downstream_canv->SaveAs(TString::Format("plots/run%d/delta_phi_upstream_downstream",(int)ion_mean_pos[run_id*16])+"_canv.png");

  //----------------------------------- statistics
  printf("\n[0] neventstoprocess: %d\n",neventstoprocess);
  printf("[1] n_good_upsteam_events: %d\n",n_good_upsteam_events);
  printf("[2] n_good_events: %d\n",n_good_events);
  printf("[3] n_fragmentation_events: %d\n",n_fragmentation_events);
  printf("\n");

  //------------------------------------------ END of the program -------------------------------//
  printf("\n\nAll done!!!\n\n");
  return;
}
