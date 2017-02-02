#include <vector>
#include <iostream>
#include <cmath>
#include "TSystem.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TObject.h"
#include "TStyle.h"
#include "TProfile.h"
#include <fstream>
#include <sstream>
#include <string>
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
//#include "DmpSoftware/trunk/Event/include/DmpStkEventMetadata.h"
#include "DmpSoftware/trunk/Event/include/DmpStkLadderAdc.hh"
#include "DmpSoftware/trunk/Event/include/DmpStkLadderCalibration.hh"
#include "DmpSoftware/trunk/Event/include/DmpStkSiCluster.h"
#include "DmpSoftware/trunk/Event/include/DmpStkStep.hh"
#include "DmpSoftware/trunk/Event/Ams/include/Scintillator.h"
#include "DmpSoftware/trunk/Event/Ams/include/ScintillatorCluster.h"


#define NScint 4

//GitHub test/1/

void STK_scint(char filename[200], char runnumber[200]){

  gSystem->Load("DmpSoftware/trunk/Event/libDmpEvent.so");

  char prefix[100]="Results/out_";
  char suffix[100]=".root";
  
  TFile* f = new TFile(filename,"READ");
  TFile* out = new TFile(strcat(strcat(prefix,runnumber),suffix), "recreate");
  
  TTree* t =(TTree*)f->Get("CollectionTree");
  cout << endl;
  f->ls();
  cout << endl;

  //Register SCINTILLATOR collections
  TClonesArray* scintillators = new TClonesArray("ScintillatorCluster");  
  t->SetBranchAddress("ScintillatorClusterCollection",&scintillators);
    
 int NEvents  = t->GetEntries();

 
 
 cout << "Number of events: " << NEvents <<endl;
 cout << endl;

 cout <<"Scintillators entries: " << scintillators->GetSize()<< endl;
 

 
 //////////////////Histos

 //scint test 
  TH1D *hScint[NScint];
 for(int i=0; i<NScint; i++){
   hScint[i] = new TH1D(Form("hScint_%d",i), Form("Energy scint %d",i),20,0,200);
   hScint[i]->GetXaxis()->SetTitle("ADC");
 }//close for1
 
 //scint tot test 
  TH1D *hScintClust[NScint];
 for(int i=0; i<NScint; i++){
   hScintClust[i] = new TH1D(Form("hScintClust_%d",i), Form("Energy Total scint %d",i),50,400,10000);
   hScintClust[i]->GetXaxis()->SetTitle("ADC");
 }//close for

 //scint ID
 TH1D *hScintID = new TH1D("scint_ID","Scint ID",8,-0.5,7.5);
 hScintID->GetXaxis()->SetTitle("ADC");


 // EVENT LOOP //

 cout << "Processing percentage: " << "00" << "% \n";

 for(int ientry=0; ientry < t->GetEntries(); ientry++){
   t->GetEntry(ientry);
      
   if ( (ientry%(NEvents/10+1))==int(NEvents/10) ) cout << "Processing percentage: " << 10.*int((ientry+1)/int(NEvents/10)) << "% \n";


   //CLUSTER LOOP//
	  

   
   for(int iscint=0; iscint < NScint; iscint++){
         ScintillatorCluster* scintillator = (ScintillatorCluster*)scintillators->ConstructedAt(iscint);
         for(int j=0; j<scintillator->GetSize(); j++){
	   
	   // cout << scintillator->GetSignal(j) << endl;
	   // cout << scintillator->GetPosition(j) << endl;
	   // cout << scintillator->GetScintillatorID() << endl;
	     hScintClust[iscint]->Fill(scintillator->GetSignalTotal());
	     hScint[iscint]->Fill(scintillator->GetSignal(j));
	     hScintID->Fill(scintillator->GetScintillatorID());
   	 }
      }
    
    if (ientry == int(NEvents-1)) cout << "Processing percentage: " << "100" << "% \n";
  
 }//close event loop
 
  TCanvas *c3 = new TCanvas("c3", "scint test ", 1000,1000);
  c3->Divide(2,2);
 for(int i=0; i<NScint;i++){
   c3->cd(i+1);
   //gPad->SetLogy();
   hScint[i]->Draw();
   hScint[i]->Write();
 }//close for
 c3->SaveAs("Results/scint.pdf"); 
 
  TCanvas *c4 = new TCanvas("c4", "scint total test ", 1000,1000);
  c4->Divide(2,2);
 for(int i=0; i<NScint;i++){
   c4->cd(i+1);
   //gPad->SetLogy();
   hScintClust[i]->Draw();
   hScintClust[i]->Write();
 }//close for
 c4->SaveAs("Results/scintTot.pdf");

 TCanvas *c5 = new TCanvas("c5", "Scint ID", 1000, 1000);
 hScintID->Draw();
 hScintID->Write();
 c5->SaveAs("Results/scintID.pdf");
 
}//close macro
