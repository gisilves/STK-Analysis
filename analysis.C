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

int run_id = 43;
TString ion_name[16]	 = {	"run"	,"H"	,"He"	,"Li"	,"Be"	,"B"	,"C"	,"N"	,"O"	,"F"	,"Ne"	,"Na"	,"Mg"	,"Al"	,"Si"	,"P"    }; // run_id
double ion_mean_pos[656] = {	00	,75.4	,-1	,-1	,-1	,-1	,-1	,-1	,-1	,-1	,-1	,-1	,-1	,-1	,-1	,-1     ,  // 0
       			   	48	,-1	,44.3	,-1	,84.95	,104.37	,122.06	,138.33	,153.54	,167.11	,178	,-1	,-1	,-1	,-1	,-1     ,  // 1
				50	,-1	,18.99	,26.2	,35.7	,44.82	,52.58	,59.6	,66.12	,71.8	,78.56	,83.5	,-1	,-1	,-1	,-1     ,  // 2 
				52	,-1	,19.41	,28.6	,36.6	,45.2	,53.07	,60.29	,67	,73.07	,78.48	,-1	,-1	,-1	,-1	,-1     ,  // 3 
				53	,-1	,-1	,-1	,36.7	,45.38	,53.23	,60.44	,67.6	,-1	,-1	,-1	,-1	,-1	,-1	,-1     ,  // 4 
				55	,-1	,19.8	,27.22	,36.7	,45.45	,53.39	,60.58	,67.74	,73.85	,79.75	,-1	,-1	,-1	,-1	,-1     ,  // 5 
				59	,-1	,19.24	,27.6	,36.8	,45.41	,53.26	,60.39	,67.9	,-1	,-1	,-1	,-1	,-1	,-1	,-1     ,  // 6 
				61	,-1	,-1	,-1	,36.93	,45.79	,53.64	,61.09	,67.62	,72.8	,-1	,-1	,-1	,-1	,-1	,-1     ,  // 7 
				71	,-1	,-1	,-1	,-1	,71.9	,86.6	,97.3	,107.28	,116.5	,123.57	,131.41	,-1	,-1	,-1	,-1     ,  // 8
				72	,-1	,-1	,-1	,60.6	,74.26	,86.11	,97.21	,107.11	,114.23	,125.32	,-1	,-1	,-1	,-1	,-1     ,  // 9
				74	,-1	,-1	,-1	,38.79	,47.44	,55.53	,63.05	,69.88	,76.13	,82.31	,87.8	,92.9	,97.98	,102.52	,106.22	,  // 10
				77	,-1	,-1	,28.3	,38.7	,47.77	,55.94	,63.43	,70.5	,76.69	,83.09	,88.5	,93.68	,98.38	,103.29	,109.94	,  // 11
				81	,-1	,-1	,-1	,38.1	,48.21	,56.47	,63.97	,70.91	,77.5	,83.43	,88.96	,94.28	,98.7	,-1	,-1     ,  // 12
				85	,-1	,-1	,-1	,37.6	,48.05	,56.06	,63.82	,70.86	,77.13	,83.56	,88.8	,94.42	,99.65	,104.11	,-1     ,  // 13
				86	,-1	,-1	,28	,39.21	,47.84	,56.03	,63.59	,70.59	,76.78	,83.01	,88.5	,93.78	,98.37	,-1	,-1     ,  // 14
				87	,-1	,-1	,-1	,39.2	,48.97	,56.2	,63.46	,70.45	,76.62	,82.77	,88.28	,93.45	,98.46	,102.96	,-1     ,  // 15
				89	,-1	,-1	,-1	,39	,47.89	,56	,63.38	,70.24	,76.76	,82.79	,88.5	,93.52	,98.34	,102.95	,107.25	,  // 16
				90	,-1	,-1	,-1	,-1	,-1	,54.3	,63.3	,70.79	,76.93	,83.2	,88.5	,94.2	,98.8	,103.35	,-1     ,  // 17
				91	,-1	,-1	,-1	,-1	,-1	,-1	,62.7	,71.08	,77.5	,83.52	,89.07	,94.36	,99.5	,103.76	,107.85	,  // 18
				93	,-1	,-1	,-1	,-1	,-1	,-1	,64.3	,71.49	,77.68	,83.94	,89.5	,94.88	,99.8	,104.23	,-1     ,  // 19
				96      ,-1     ,-1     ,-1     ,-1     ,-1     ,-1     ,65.6   ,72.79  ,79.24  ,85.45  ,91.2   ,96.51  ,101.74 ,105.97 ,-1     ,  // 20
				97      ,-1     ,-1     ,-1     ,-1     ,-1     ,-1     ,64.7   ,72.73  ,79.02  ,85.6   ,91.1   ,96.33  ,101.47 ,105.92 ,-1     ,  // 21
				100     ,-1     ,-1     ,-1     ,-1     ,-1     ,-1     ,64.6   ,71.7   ,78.13  ,84.44  ,89.8   ,95.08  ,100.5  ,104.96 ,-1     ,  // 22
				101     ,-1     ,-1     ,-1     ,-1     ,-1     ,-1     ,64     ,71.84  ,78.5   ,83.96  ,89.5   ,96.84  ,-1     ,104.82 ,-1     ,  // 23
				106     ,-1     ,-1     ,-1     ,-1     ,-1     ,-1     ,64.3   ,71.67  ,78.03  ,84.26  ,89.94  ,95.09  ,99.93  ,104.75 ,109.29 ,  // 24
				107     ,-1     ,-1     ,-1     ,39.55  ,48.5   ,56.9   ,64.46  ,71.58  ,78.24  ,83.9   ,88.8   ,94.74  ,99.8   ,104.4  ,109.04 ,  // 25
				109     ,-1     ,-1     ,-1     ,39.6   ,48.71  ,57.02  ,64.69  ,71.71  ,78.1   ,84.44  ,90.01  ,95.25  ,99.8   ,104.79 ,109.22 ,  // 26
				111     ,-1     ,-1     ,-1     ,39.6   ,49.2   ,57.68  ,65.35  ,72.58  ,79.2   ,85.54  ,91.5   ,96.5   ,101.5  ,105.86 ,-1     ,  // 27
				115     ,-1     ,-1     ,-1     ,-1     ,-1     ,56.59  ,65.3   ,73.13  ,79.16  ,85.72  ,91.34  ,96.82  ,101.79 ,106.44 ,109.63 ,  // 28
				117     ,-1     ,-1     ,-1     ,-1     ,-1     ,59.2   ,66.22  ,73.18  ,79.67  ,86.2   ,91.5   ,97.2   ,102.32 ,106.96 ,-1     ,  // 29
				118     ,-1     ,-1     ,-1     ,40.61  ,49.63  ,58.21  ,66.08  ,73.15  ,79.8   ,86.19  ,91.74  ,97.3   ,101.99 ,106.93 ,-1     ,  // 30
				120     ,-1     ,-1     ,-1     ,40     ,49.84  ,58.32  ,66.17  ,73.49  ,80.13  ,86.29  ,92.12  ,97.5   ,101.8  ,105.8  ,-1     ,  // 31
				121     ,-1     ,-1     ,-1     ,-1     ,-1     ,59.5   ,66.59  ,73.83  ,80.4   ,86.58  ,92.37  ,96.8   ,99.8   ,106.83 ,-1     ,  // 32
				128     ,-1     ,-1     ,-1     ,-1     ,-1     ,58.6   ,67.07  ,74.21  ,80.69  ,87.2   ,92.79  ,98.34  ,103.5  ,108.07 ,-1     ,  // 33
				139     ,-1     ,-1     ,-1     ,-1     ,-1     ,55.8   ,64.7   ,72.23  ,78.66  ,84.8   ,90.49  ,95.87  ,100.66 ,105.25 ,-1     ,  // 34
				141     ,-1     ,-1     ,30.6   ,39.83  ,48.8   ,57.22  ,64.93  ,72.03  ,78.49  ,84.79  ,90.23  ,95.68  ,100.8  ,105.33 ,109.58 ,  // 35
				143     ,-1     ,-1     ,29.3   ,39     ,49.05  ,57.47  ,65.31  ,72.4   ,78.98  ,85.05  ,90.5   ,96.2   ,101.03 ,106.62 ,109.5  ,  // 36
				145     ,-1     ,-1     ,-1     ,-1     ,-1     ,57.9   ,65.6   ,72.5   ,79.2   ,85.2   ,91.16  ,96.39  ,101.37 ,106.01 ,110.24 ,  // 37
				148     ,-1     ,-1     ,-1     ,40.47  ,49.48  ,57.88  ,65.67  ,72.84  ,79.23  ,85.84  ,91.29  ,96.81  ,102.18 ,106.67 ,111.53 ,  // 38
				149     ,-1     ,-1     ,-1     ,41.04  ,49.59  ,58.04  ,65.69  ,72.88  ,79.5   ,85.9   ,91.47  ,97.08  ,102.37 ,107.06 ,-1     ,  // 49
				155     ,-1     ,20.89  ,29.97  ,40.4   ,49.78  ,58.63  ,66     ,73.48  ,79.76  ,86.51  ,92.2   ,97.55  ,103.27 ,107.83 ,111.91 }; // 40

//---------------------------------------- track fitting -------------------------------------//

// define a track function - straigt line
double func(float z, double *par)
{
  double value = par[0]*z + par[1];
  return value;
}

// define global data points
//const Int_t ndatapoints = 6;
//double pos_xy[ndatapoints],pos_z[ndatapoints],pos_xy_error[ndatapoints];
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

  // // Print results
  // Double_t amin,edm,errdef;
  // Int_t nvpar,nparx,icstat;
  // if (debug) {
  //   gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  //   gMinuit->mnprin(3,amin); }

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
	      bool dbg = true
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
  //DmpStkEventMetadata* stkMetadata = new DmpStkEventMetadata();
  //t->SetBranchAddress("DmpStkEventMetadata", &stkMetadata);
  TClonesArray* stkladderadc = new TClonesArray("DmpStkLadderAdc");
  t->SetBranchAddress("DmpStkLadderAdcCollection", &stkladderadc);

  // // Register ANCILLARY container
  // AncillaryEventIons* ancevent = new AncillaryEventIons();
  // t->SetBranchAddress("AncillaryEventIons",&ancevent);

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

    // Y telescope planes //to do: change id numbers
    non_aligned_corr[0 ]->Fill( 93-stkcluster_maxehitposition[0 ], 93-stkcluster_maxehitposition[2 ]); // 0   -   6
    non_aligned_corr[1 ]->Fill( 93-stkcluster_maxehitposition[2 ],    stkcluster_maxehitposition[4 ]); // 6   -  72
    non_aligned_corr[2 ]->Fill(    stkcluster_maxehitposition[4 ],    stkcluster_maxehitposition[14]); // 72  -  96
    non_aligned_corr[3 ]->Fill(    stkcluster_maxehitposition[14],    stkcluster_maxehitposition[16]); // 96  - 126
    non_aligned_corr[4 ]->Fill(    stkcluster_maxehitposition[16],    stkcluster_maxehitposition[18]); // 126 - 120
    // X telescope planes
    non_aligned_corr[5 ]->Fill( 93-stkcluster_maxehitposition[1 ], 93-stkcluster_maxehitposition[3 ]); // 1   -   7
    non_aligned_corr[6 ]->Fill( 93-stkcluster_maxehitposition[3 ],    stkcluster_maxehitposition[5 ]); // 7   -  73
    non_aligned_corr[7 ]->Fill(    stkcluster_maxehitposition[5 ],    stkcluster_maxehitposition[15]); // 73  -  97
    non_aligned_corr[8 ]->Fill(    stkcluster_maxehitposition[15],    stkcluster_maxehitposition[17]); // 97  - 127
    non_aligned_corr[9 ]->Fill(    stkcluster_maxehitposition[17],    stkcluster_maxehitposition[19]); // 127 - 121
    // Y DUT planes
    non_aligned_corr[10]->Fill(    stkcluster_maxehitposition[5],     stkcluster_maxehitposition[6 ]); // 73  -  78
    non_aligned_corr[11]->Fill(    stkcluster_maxehitposition[5],     stkcluster_maxehitposition[7 ]); // 73  -  79
    non_aligned_corr[12]->Fill(    stkcluster_maxehitposition[5],     stkcluster_maxehitposition[8 ]); // 73  - 175
    non_aligned_corr[13]->Fill(    stkcluster_maxehitposition[5],     stkcluster_maxehitposition[9 ]); // 73  - 102
    non_aligned_corr[14]->Fill(    stkcluster_maxehitposition[5],     stkcluster_maxehitposition[10]); // 73  - 103
    non_aligned_corr[15]->Fill(    stkcluster_maxehitposition[5],     stkcluster_maxehitposition[11]); // 73  - 169
    non_aligned_corr[16]->Fill(    stkcluster_maxehitposition[5],     stkcluster_maxehitposition[12]); // 73  - 168
    non_aligned_corr[17]->Fill(    stkcluster_maxehitposition[5],     stkcluster_maxehitposition[13]); // 73  - 174
    
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
  TH1F *scint_cluster_q[4];
  for (int i=0; i<4; i++) scint_cluster_q[i] = new TH1F("scint_cluster_q_ch"+TString::Format("%d",i),
  							"scint_cluster_q_ch"+TString::Format("%d",i),
  							600,0,400);
  TH1F *scint_cluster_pos[4];
  for (int i=0; i<4; i++) scint_cluster_pos[i] = new TH1F("scint_cluster_pos_ch"+TString::Format("%d",i),
  							  "scint_cluster_pos_ch"+TString::Format("%d",i),
  							  300,0,4000);
  TH1F *scint_cluster_width[4];
  for (int i=0; i<4; i++) scint_cluster_width[i] = new TH1F("scint_cluster_width_ch"+TString::Format("%d",i),
  							    "scint_cluster_width_ch"+TString::Format("%d",i),
  							    300,0,600);
  TH2F *stkmaxq_vs_scint_cluster_q[4];
  for (int i=0; i<4; i++) stkmaxq_vs_scint_cluster_q[i] = new TH2F("stkmaxq_vs_scint_cluster_q_ch"+TString::Format("%d",i),
  								   "stkmaxq_vs_scint_cluster_q_ch"+TString::Format("%d",i),
  								   300,0,200,600,0,400);
  //-------------------------------------------- tracking plots
  TH1F *delta_phi_upstream_downstream_x_tracks = new TH1F("delta_phi_upstream_downstream_x_tracks","delta_phi_upstream_downstream_x_tracks",200,-2,2);
  TH1F *delta_phi_upstream_downstream_y_tracks = new TH1F("delta_phi_upstream_downstream_y_tracks","delta_phi_upstream_downstream_y_tracks",200,-2,2);
  
  /////////////////////////////////// Event LOOP //////////////////////////////////

  if (neventstoprocess==-1 || neventstoprocess > t->GetEntries()) neventstoprocess = t->GetEntries();
  for(int entry=0; entry<neventstoprocess; entry++){
   
    t->GetEntry(entry);

    // printf("@@@@@@@@@@ Event %d\n", entry);
    // // STK metadata
    // printf("STK data mode = %d\n", stkMetadata->fRunMode );
    // // STK modes: 2 (COMPRESSED), 3(RAW), 5(ULD),  6(DLD)
    // // ... for more information see: DmpStkEventMetadata.h

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

      // printf("\nSTK cluster info:\n");
      // printf("   total ADC counts = %f\n",cluster->getEnergy());
      // printf("   number of strips = %d\n",cluster->getNstrip());         
      // printf("   signal in the first strip = %f\n", cluster->GetSignal(0));
      // printf("   signal in the last  strip = %f\n", cluster->GetSignal(cluster->getNstrip()-1));
      // printf("   noise  in the first strip = %f\n", cluster->GetNoise(0));
      // printf("   noise  in the last  strip = %f\n", cluster->GetNoise(cluster->getNstrip()-1));
      // printf("   center of gravity = %f\n",cluster->getLadderCentroid());
      // printf("   ladder ID (hardware) = %d\n",cluster->getLadderHardware());
      //etc.
      // For more details on STK clusters see DmpStkSiCluster.h

      //------------------------- remove clusters with bad channels
      int firststripid = cluster->getFirstStrip();
      for (int i=firststripid; i<firststripid+stkclusternstr; i++) {
  	//!!!!!!! TO DO !!!!!!!!!!!
  	// stkludderid
  	;
      }

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
    printf("nscint %d\n", scintclusters->GetLast()+1);
    for(int icl=0; icl<scintclusters->GetLast()+1; icl++)
      {
	ScintillatorCluster* cluster = (ScintillatorCluster*)scintclusters->ConstructedAt(icl);
	cout<<cluster->GetScintillatorID()<<endl; 
      }


    // Loop over reconstructed SCINTILLATOR Clusters
    for (int i=0; i<n_scint_clusters; i++) {
      
      ScintillatorCluster* cluster = (ScintillatorCluster*)scintclusters->ConstructedAt(i);
      int scintid = cluster->GetScintillatorID(); 
      int scintclusterfirstsampleid = cluster->GetPosition(0);
      int scintclustersize = cluster->GetSize();
      scint_cluster_charge[scintid] = sqrt(cluster->GetSignalTotal());

      scint_cluster_pos[scintid]->Fill(scintclusterfirstsampleid);

      for(int j=0; j<scintclustersize; j++){
  	if (cluster->GetSignal(j) > 3954) saturated_signal[scintid] = true; }

      if (scintclusterfirstsampleid<1240 || scintclusterfirstsampleid>1300) {
  	non_trigger_particle = true;
  	break; }
    }// loop over scintillator clusters
	  	  
    //--------------------- skipp events if there is no signal in either of the scintillators
    if (n_scint_clusters!=4) continue;

    //--------------------- skipp events with extra particles
    if (non_trigger_particle) continue;
   
    n_single_particle_trigger_events++;

    //--------------------- skipp events if scintilator signals in ch2 are suturated
    if (saturated_signal[2]) continue;

    // //--------------------- choose ions based on energy of scintillator clusters
    //if (scint_cluster_charge[2]<55 || scint_cluster_charge[2]>63) continue;
    
    // double sum;
    // double average_q_x_except_max[6];
    // double average_q_y_except_max[14];
    // sum = 0.0; for (uint i=0; i<stk_clusters_x_e[0].size(); i++)  sum += stk_clusters_x_e[0].at(i);  average_q_x_except_max[0]  = (sum-stkcluster_maxe[0 ])/(stk_clusters_x_e[0 ].size()-1);

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
      chi2_y_track_6->Fill(log10(y_track_fit_results_6.at(4)));}
    // correlations
    if (!missing_y_hit_3a && !missing_y_hit_3b) chi2_y_track_3a_vs_3b->Fill(log10(y_track_fit_results_3a.at(4)),log10(y_track_fit_results_3b.at(4)));
    if (!missing_y_hit_3a && !missing_y_hit_6)  chi2_y_track_3a_vs_6 ->Fill(log10(y_track_fit_results_3a.at(4)),log10(y_track_fit_results_6 .at(4)));
    if (!missing_y_hit_3b && !missing_y_hit_6)  chi2_y_track_3b_vs_6 ->Fill(log10(y_track_fit_results_3b.at(4)),log10(y_track_fit_results_6 .at(4)));
    
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
    if (//stknclusters[0]<max_n_cluster*2 && 
  	//stknclusters[2]<max_n_cluster && 
  	//stknclusters[4]<max_n_cluster &&
  	log10(secondary_clusters_energy_ratio[0])<max_logeratio+0.2 &&
  	log10(secondary_clusters_energy_ratio[2])<max_logeratio &&
  	log10(secondary_clusters_energy_ratio[4])<max_logeratio
  	) x_upstream_track_is_clean = true;
    if (//stknclusters[1]<max_n_cluster*2 && 
  	//stknclusters[3]<max_n_cluster && 
  	//stknclusters[5]<max_n_cluster &&
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
      /*
      // cluster charge in X telescope planes
      cluster_q_x_track[0]->Fill(sqrt(stkcluster_maxe[0 ])); 
      cluster_q_x_track[1]->Fill(sqrt(stkcluster_maxe[2 ])); 
      cluster_q_x_track[2]->Fill(sqrt(stkcluster_maxe[4 ])); 
      cluster_q_x_track[3]->Fill(sqrt(stkcluster_maxe[14])); 
      cluster_q_x_track[4]->Fill(sqrt(stkcluster_maxe[16])); 
      cluster_q_x_track[5]->Fill(sqrt(stkcluster_maxe[18]));
      // cluster charge in Y telescope planes
      cluster_q_y_track[0]->Fill(sqrt(stkcluster_maxe[1 ])); 
      cluster_q_y_track[1]->Fill(sqrt(stkcluster_maxe[3 ])); 
      cluster_q_y_track[2]->Fill(sqrt(stkcluster_maxe[5 ])); 
      cluster_q_y_track[3]->Fill(sqrt(stkcluster_maxe[15])); 
      cluster_q_y_track[4]->Fill(sqrt(stkcluster_maxe[17])); 
      cluster_q_y_track[5]->Fill(sqrt(stkcluster_maxe[19])); 
      */
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
      //  For each wave, clusters are sorted in time, the first one in the collection comes first in the event
      //  There are 2 scintillators, for each of them you have amplified and non-amplified signals:
      //         wave0  (cluster->GetScintillatorID()==0)   is a non-amplified signal of first scintillator
      //         wave1  (cluster->GetScintillatorID()==1)   is an amplified signal of first scintillator
      //         wave2  (cluster->GetScintillatorID()==2)   is a non-amplified signal of second scintillator
      //         wave3  (cluster->GetScintillatorID()==3)   is an amplified signal of second scintillator
      //       
      //  scintillator class is defined here: http://dpnc.unige.ch/SVNDAMPE/DAMPE/DmpSoftware/DmpSoftware-4-5-4/Event/Ams/include/ScintillatorCluster.h      

    }// cuts based on track fitting


    // // Loop over STK adc 
    // for(int i=0; i<stkladderadc->GetLast()+1; i++){
    //   DmpStkLadderAdc* laddderadc = (DmpStkLadderAdc*)stkladderadc->ConstructedAt(i);
    //   int ladderid = laddderadc->GetLadderID();
    //   for(int j=0; j<384; j++){
    // 	double adc = laddderadc->GetChannelAdc(j);
    //   }
    // }
  
    //   // Ancillary detector
    //   printf("\nAnclillary detector information:\n");
    //   printf("   ADC s1 low gain   = %f\n",ancevent->s1_6db);
    //   printf("   ADC s1 high gain (x10)  = %f\n",ancevent->s1_12db);
    //   printf("   ADC s2 low gain  = %f\n",ancevent->s2_6db);
    //   printf("   ADC s2 high gain (x10) = %f\n",ancevent->s2_12db);
    //   printf("   ADC si0 low gain = %f\n",ancevent->si0_sg);
    //   printf("   ADC si0 high gain (x10) = %f\n",ancevent->si0_lg);
    //   printf("   ADC si1 low gain = %f\n",ancevent->si1_sg);
    //   printf("   ADC si1 high gain () = %f\n",ancevent->si1_lg);
    //   // etc.
    //   // For more information see http://dpnc.unige.ch/SVNDAMPE/DAMPE/DmpSoftware/DmpSoftware-4-5-4/Event/Ams/include/AncillaryEventIons.hh
  
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
  /*
  TCanvas *cluster_q_x_track_canv = new TCanvas("cluster_q_x_track_canv","cluster_q_x_track_canv",1800,1000);
  cluster_q_x_track_canv->Divide(3,2);
  for (int i=0; i<6; i++) {cluster_q_x_track_canv->cd(i+1); cluster_q_x_track[i]->Draw();}
  cluster_q_x_track_canv->SaveAs(TString::Format("plots/run%d/cluster_q_x_track",(int)ion_mean_pos[run_id*16])+"_canv.png");

  TCanvas *cluster_q_y_track_canv = new TCanvas("cluster_q_y_track_canv","cluster_q_y_track_canv",1800,1000);
  cluster_q_y_track_canv->Divide(3,2);
  for (int i=0; i<6; i++) {cluster_q_y_track_canv->cd(i+1); cluster_q_y_track[i]->Draw();}
  cluster_q_y_track_canv->SaveAs(TString::Format("plots/run%d/cluster_q_y_track",(int)ion_mean_pos[run_id*16])+"_canv.png");
  */
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
  TCanvas *scint_cluster_q_canv = new TCanvas("scint_cluster_q_canv","scint_cluster_q_canv",1800,1000);
  scint_cluster_q_canv->Divide(2,2);
  for (int i=0; i<4; i++){
    scint_cluster_q_canv->cd(i+1);
    if (i==0 || i==2) scint_cluster_q[i]->GetXaxis()->SetRangeUser(0,200);
    scint_cluster_q[i]->Draw();
    
    if (i==2) { // scint. ch2 - non amp., upstream
      
      TString multi_gaus_fitf_str = "";
      int gaus_first_par_id = 0;
      for (int j=1; j<16; j++) {
  	if (ion_mean_pos[run_id*16+j]!=-1) {
  	  if (multi_gaus_fitf_str!="") multi_gaus_fitf_str += "+";
  	  multi_gaus_fitf_str += "gaus("+TString::Format("%d",gaus_first_par_id)+")";
  	  gaus_first_par_id += 3; }}
      TF1 *multi_gaus_fitf = new TF1("multi_gaus_fitf",multi_gaus_fitf_str);
      double upper_limit = -999.99;
      double lower_limit = -999.99;
      gaus_first_par_id = 0;
      for (int j=1; j<16; j++) {
  	if (ion_mean_pos[run_id*16+j]!=-1) {
  	  // set initial values for parameters
  	  multi_gaus_fitf->SetParameter(gaus_first_par_id,1);
  	  multi_gaus_fitf->SetParameter(gaus_first_par_id+1,ion_mean_pos[run_id*16+j]-4);
  	  multi_gaus_fitf->SetParameter(gaus_first_par_id+2,0.6);
  	  // set limits
  	  multi_gaus_fitf->SetParLimits(gaus_first_par_id+1,ion_mean_pos[run_id*16+j]-4,ion_mean_pos[run_id*16+j]+1.5);
  	  multi_gaus_fitf->SetParLimits(gaus_first_par_id+2,0.6,2.5);
  	  if (lower_limit==-999.99) lower_limit = ion_mean_pos[run_id*16+j]-3;
  	  upper_limit = ion_mean_pos[run_id*16+j]+1;
  	  gaus_first_par_id += 3; }}

      multi_gaus_fitf->SetRange(lower_limit,upper_limit);
      TVirtualFitter::SetMaxIterations(50000);
      scint_cluster_q[i]->Fit(multi_gaus_fitf,"RLM");

      printf("\n");
      printf("%4s %10s %10s %10s\n","ion","const","mean","sigma");
      gaus_first_par_id = 0;
      for (int j=1; j<16; j++) {
  	if (ion_mean_pos[run_id*16+j]!=-1) {
  	  TString ion = ion_name[j];
  	  double norm = multi_gaus_fitf->GetParameter(gaus_first_par_id);
  	  double mean = multi_gaus_fitf->GetParameter(gaus_first_par_id+1);
  	  double sigma = multi_gaus_fitf->GetParameter(gaus_first_par_id+2);
  	  gaus_first_par_id += 3;
  	  printf("%4s %10.2f %10.2f %10.2f\n",ion.Data(),norm,mean,sigma); }}
      printf("\n");
    }
  }
  scint_cluster_q_canv->SaveAs(TString::Format("plots/run%d/scint_cluster_q",(int)ion_mean_pos[run_id*16])+"_canv.png");

  TCanvas *scint_cluster_pos_canv = new TCanvas("scint_cluster_pos_canv","scint_cluster_pos_canv",1800,1000);
  scint_cluster_pos_canv->Divide(2,2);
  for (int i=0; i<4; i++){
    scint_cluster_pos_canv->cd(i+1);
    scint_cluster_pos_canv->cd(i+1)->SetLogy();
    scint_cluster_pos[i]->Draw(); }  
  scint_cluster_pos_canv->SaveAs(TString::Format("plots/run%d/scint_cluster_pos",(int)ion_mean_pos[run_id*16])+"_canv.png");

  TCanvas *scint_cluster_width_canv = new TCanvas("scint_cluster_width_canv","scint_cluster_width_canv",1800,1000);
  scint_cluster_width_canv->Divide(2,2);
  for (int i=0; i<4; i++){
    scint_cluster_width_canv->cd(i+1);
    scint_cluster_width[i]->Draw(); }  
  scint_cluster_width_canv->SaveAs(TString::Format("plots/run%d/scint_cluster_width",(int)ion_mean_pos[run_id*16])+"_canv.png");

  TCanvas *stkmaxq_vs_scint_cluster_q_canv = new TCanvas("stkmaxq_vs_scint_cluster_q_canv","stkmaxq_vs_scint_cluster_q_canv",1800,1000);
  stkmaxq_vs_scint_cluster_q_canv->Divide(2,2);
  for (int i=0; i<4; i++){
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
  printf("[1] n_single_particle_trigger_events: %d\n",n_single_particle_trigger_events);
  printf("[2] n_good_upsteam_events: %d\n",n_good_upsteam_events);
  printf("[3] n_good_events: %d\n",n_good_events);
  printf("[4] n_fragmentation_events: %d\n",n_fragmentation_events);
  printf("[1]/[0]: %.4f\n",n_single_particle_trigger_events/(double)neventstoprocess);
  printf("[2]/[1]: %.4f\n",n_good_upsteam_events/(double)n_single_particle_trigger_events);
  printf("[3]/[2]: %.4f\n",n_good_events/(double)n_good_upsteam_events);
  printf("[4]/[2]: %.4f\n",n_fragmentation_events/(double)n_good_upsteam_events);
  //printf(": %d\n",);
  printf("\n");

  //------------------------------------------ END of the program -------------------------------//
  printf("\n\nAll done!!!\n\n");
  return;
}
