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
#include "TLine.h"
#include <fstream>
#include <sstream>
#include <string>


//#include "DmpSoftware/trunk/Event/include/DmpStkEventMetadata.h"
#include "DmpSoftware/trunk/Event/include/DmpStkLadderAdc.hh"
#include "DmpSoftware/trunk/Event/include/DmpStkLadderCalibration.hh"
#include "DmpSoftware/trunk/Event/include/DmpStkSiCluster.h"
#include "DmpSoftware/trunk/Event/include/DmpStkStep.hh"
#include "DmpSoftware/trunk/Event/Ams/include/Scintillator.h"

#define NLadders 1
#define NScint 2
#define NbadChan 11

bool isBad(int arr[], int size, int value);

void STK_2d(char filename[200], int runnumber){
  gSystem->Load("DmpSoftware/trunk/Event/libDmpEvent.so");

  
 TString term_cmd = "mkdir -pv "+TString::Format("plots/run%d",(int)runnumber);
 system(term_cmd.Data());
 
 TString term_cmd2 = "mkdir -pv "+TString::Format("Results/run%d",(int)runnumber);
 system(term_cmd2.Data());

 TString plotdirectory = TString::Format("plots/run%d",(int)runnumber);
 TString run_string = TString::Format("%04d",(int)runnumber);
 
TString outfile = TString::Format("Results/run%d",(int)runnumber)+"/out_2d_l126_"+run_string+".root";

  
  TFile* f = new TFile(filename,"READ");
  TFile* out = new TFile(outfile, "recreate");

  
  TTree* t =(TTree*)f->Get("CollectionTree");
  cout << endl;
  f->ls();
  cout << endl;
  
  //Register STK collections//
  TClonesArray* stkclusters = new TClonesArray("DmpStkSiCluster");
  t->SetBranchAddress("StkClusterCollection",&stkclusters);

    
 int NEvents  = t->GetEntries();
 cout << "Number of events: " << NEvents << endl;
 cout << endl;
 
 //Available ladders ={0,1,6,7,24,25,30,31,72,73,78,79,96,97,102,103,120,121,126,127}
 //Ordered   ladders ={31,30,25,24,72,73,103,102,97,96,127,126,121,120,78,79,0,1,6,7}
 // ladder names={"1Y","1X","2Y","2X","3Y","3X","LPFM116B","LGFM072T","LGFM018B","LGFM084T","LGFM031T","LPFM999B","LGFM071T","LPFM069B","4Y","4X,"5Y","5X","6Y","6X"}

 // const char *ladderName[NLadders]={"1Y","1X","2Y","2X","3Y","3X","LPFM116B","LGFM072T","LGFM018B","LGFM084T","LGFM031T","LPFM905T","LGFM071T","LPFM069B","4Y","4X","5Y","5X","6Y","6X"};
 //int ladderID[NLadders]={31,30,25,24,72,73,103,102,97,96,127,126,121,120,78,79,0,1,6,7};

 int ladderID[NLadders]={126};
 int badChan[NbadChan]={130,134,172,189,190,191,192,193,194,224,225};
 int ladder;

 int first;
 double GetEnergy;
 double getcentroid[4];
 double MaxEne;
 int iscluster;
 double clustercog;
 int Nstrip;
 

 //////////////////Histos

 //2d test
 TH2D *h2test[NLadders];
 for(int iladder=0; iladder<NLadders;iladder++){
   h2test[iladder] = new TH2D(Form("h2dTest_ladder_%d",ladderID[iladder]), Form("2d test ladder %d ",ladderID[iladder]),384,-0.5 ,383.5,250,0,500);
   h2test[iladder]->GetXaxis()->SetTitle("strip number");
   h2test[iladder]->GetYaxis()->SetTitle("ADC");
 }//close i for
 
 //2d test
 TH2D *h2testdoppiobond = new TH2D(Form("h2dTest(hPitch)_ladder_%d",126), Form("2d test (half pitch) ladder %d ",126),32, 160.5,192.5,250,0,500);
   h2testdoppiobond->GetXaxis()->SetTitle("strip number");
   h2testdoppiobond->GetYaxis()->SetTitle("ADC");

 //2d test
 TH2D *h2test8si = new TH2D(Form("h2dTest(8Si)_ladder_%d",126), Form("2d test (8 Si) ladder %d ",126),64, 192.5,256.5,250,0,500);
   h2test8si->GetXaxis()->SetTitle("strip number");
   h2test8si->GetYaxis()->SetTitle("ADC"); 


 // EVENT LOOP //

 cout << "Processing percentage: " << "00" << "% \n";

 for(int ientry=0; ientry < t->GetEntries(); ientry++){
   t->GetEntry(ientry);

   if ( (ientry%(NEvents/10+1))==int(NEvents/10) ) cout << "Processing percentage: " << 10.*int((ientry+1)/int(NEvents/10)) << "% \n";


   //CLUSTER LOOP//
   
   for(int i =0; i<4; i++)getcentroid[i]=0;
   //reset centroids
   
   MaxEne = 0;
   iscluster=0;
   //initialize energy
   
   for(int i=0; i<stkclusters->GetLast()+1; i++){
     DmpStkSiCluster* cluster = (DmpStkSiCluster*) stkclusters-> ConstructedAt(i);
   
     GetEnergy = cluster->getEnergy();
     ladder = cluster->getLadderHardware();
     clustercog = cluster->getLadderCentroid();
     Nstrip = cluster->getNstrip();
     int FirstStrip = cluster->getFirstStrip();
     int LastStrip = cluster->getLastStrip();
     
     for(int iladder = 0; iladder < NLadders; iladder++){
       if(ladder == ladderID[iladder]){
	 double MaxStripEnergy=0;
	 double SecondStripEnergy=0;
	 int PositionMax=0;
	 int PositionOth=0;
	 
	 if(Nstrip>0){
	   for(int j = 0; j < Nstrip; j++){
	     if(cluster->GetSignal(j)>= MaxStripEnergy){  
	       SecondStripEnergy=MaxStripEnergy;
	       MaxStripEnergy=cluster->GetSignal(j);
	       PositionMax=j;
	     }
	     else if(cluster->GetSignal(j) > SecondStripEnergy){
	       SecondStripEnergy=cluster->GetSignal(j);
	       PositionOth=j;
	     }
	   }
	   
	 int ClusterStrip= FirstStrip+PositionMax;
	 if(isBad(badChan,NbadChan,ClusterStrip)==0){  
	   if(ladderID[iladder]==126){
	     if(ClusterStrip>160 && ClusterStrip<193){
	       h2testdoppiobond->Fill(ClusterStrip,GetEnergy);
	     } else if(ClusterStrip>192 && ClusterStrip<257){
	       h2test8si->Fill(ClusterStrip,GetEnergy);
	     } else if(ClusterStrip<161 || ClusterStrip>192){
	       h2test[iladder]->Fill(ClusterStrip,GetEnergy);
	     }
	   }
	 }
       }//close ladder if
       }
     }//close for
   }//close cluster loop


    if (ientry == int(NEvents-1)) cout << "Processing percentage: " << "100" << "% \n";
  
 }//close event loop
 
 TCanvas *c3 = new TCanvas("c3", "2d test ", 900,1000);
 c3->Divide(1,1);
 for(int i=0; i< NLadders;i++){
   c3->cd(i+1);
   gPad->SetLogz();
   h2test[i]->Draw("colz");
   h2test[i]->Write();
   TLine *l = new TLine(64,0,64,500);
   l->Draw("SAME");
   TLine *l1 = new TLine(128,0,128,500);
   l1->Draw("SAME");
   TLine *l2 = new TLine(192,0,192,500);
   l2->Draw("SAME");
   TLine *l3 = new TLine(256,0,256,500);
   l3->Draw("SAME");
   TLine *l4 = new TLine(320,0,320,500);
   l4->Draw("SAME");
 }//close for
 c3->SaveAs(plotdirectory+"/2d_test_126_"+run_string+".png");

 TCanvas *c4 = new TCanvas("c4", "2d test (half pitch)", 900,1000);
 gPad->SetLogz();
 h2testdoppiobond->Draw("colz");
 h2testdoppiobond->Write();
 c4->SaveAs(plotdirectory+"/2d_test_126_hp_"+run_string+".png");
   
 TCanvas *c5 = new TCanvas("c5", "2d test (8 Si)", 900,1000);
 gPad->SetLogz();
 h2test8si->Draw("colz");
 h2test8si->Write();
 c5->SaveAs(plotdirectory+"/2d_test_126_8si_"+run_string+".png");

}//close macro

bool isBad(int arr[], int size, int value){
	
	bool isbad=0;
	int i=0;
		
		while(i<size && isbad==0){
			if(arr[i]==value)isbad=1;
			i=i+1;
		}
		
   return isbad;
}

