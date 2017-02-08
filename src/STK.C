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
#include "DmpStkLadderAdc.hh"
#include "DmpStkLadderCalibration.hh"
#include "DmpStkSiCluster.h"
#include "DmpStkStep.hh"
#include "Scintillator.h"
#include "ScintillatorCluster.h"


#define NLadders 20

//GitHub test 2//



  
void STK_test(char filename[200], int runnumber){

  gSystem->Load("libDmpEvent.so");

 TString term_cmd = "mkdir -pv "+TString::Format("plots/run%d",(int)runnumber);
 system(term_cmd.Data());
 
 TString term_cmd2 = "mkdir -pv "+TString::Format("Results/run%d",(int)runnumber);
 system(term_cmd2.Data());
 
 TString run_string = TString::Format("%04d",(int)runnumber);
 TString plotdirectory = TString::Format("plots/run%d",(int)runnumber);
 TString outfile = TString::Format("Results/run%d",(int)runnumber)+"/out_"+run_string+".root";

  
  
  TFile* f = new TFile(filename,"READ");
  TFile* out = new TFile(outfile, "recreate");
  
  TTree* t =(TTree*)f->Get("CollectionTree");
  cout << endl;
  f->ls();
  cout << endl;
  
  //Register STK collections//
  TClonesArray* stkclusters = new TClonesArray("DmpStkSiCluster");
  t->SetBranchAddress("StkClusterCollection",&stkclusters);

  //Register SCINTILLATOR collections
  TClonesArray* scintillators = new TClonesArray("ScintillatorCluster");  
  t->SetBranchAddress("ScintillatorClusterCollection",&scintillators);
    
 int NEvents  = t->GetEntries();
 
 cout << "Number of events: " << NEvents <<endl;
 cout << endl;
 
 //Available ladders ={0,1,6,7,24,25,30,31,72,73,78,79,96,97,102,103,120,121,126,127}
 //Ordered   ladders ={31,30,25,24,72,73,103,102,97,96,127,126,121,120,78,79,0,1,6,7}
 // ladder names={"1Y","1X","2Y","2X","3Y","3X","LPFM116B","LGFM072T","LGFM018B","LGFM084T","LGFM031T","LPFM999B","LGFM071T","LPFM069B","4Y","4X,"5Y","5X","6Y","6X"}

 const char *ladderName[NLadders]={"1Y","1X","2Y","2X","3Y","3X","LPFM116B","LGFM072T","LGFM001B","LGFM084T","LGFM031T","LPFM999B","LGFM071T","LPFM069B","4Y","4X","5Y","5X","6Y","6X"};
 int ladderID[NLadders]={31,30,25,24,72,73,103,102,97,96,127,126,121,120,78,79,0,1,6,7};
  
 int ladder;

 int Nclust[NLadders]={};

 int first;
 double GetEnergy;
 double getcentroid[4];
 double MaxEne;
 int iscluster;
 double clustercog;
 int Nstrip;
 
 //////////////////Histos
 TH1D *hEnergyCluster[NLadders];
 for(int iladder =0; iladder < NLadders; iladder++){
hEnergyCluster[iladder] = new TH1D(Form("hEnergyCluster_iladder_%d",ladderID[iladder]),ladderName[iladder],1000,0,200);
   hEnergyCluster[iladder]->GetXaxis()->SetTitle("ADC");
 }//close for
 
 TH1D *hEnergyClusterMin[NLadders];
 for(int iladder =0; iladder < NLadders; iladder++){
hEnergyClusterMin[iladder] = new TH1D(Form("hEnergyCluster__min_iladder_%d",ladderID[iladder]),ladderName[iladder],1000,0,200);
   hEnergyClusterMin[iladder]->GetXaxis()->SetTitle("ADC");
 }//close for

 TH1D *hEnergyClusterMax[NLadders];
 for(int iladder =0; iladder < NLadders; iladder++){
hEnergyClusterMax[iladder] = new TH1D(Form("hEnergyCluster__max_iladder_%d",ladderID[iladder]),ladderName[iladder],100,-0.5,20.5);
   hEnergyClusterMax[iladder]->GetXaxis()->SetTitle("ADC");
 }//close for
 
 
  TH1D *hEnergyOneStrip[NLadders];
 for(int iladder =0; iladder < NLadders; iladder++){
   hEnergyOneStrip[iladder] = new TH1D(Form("One Strip hEnergyCluster_iladder_%d",ladderID[iladder]), Form("One Strip Energy clusters iladder %d",ladderID[iladder]),100,0,400);
   hEnergyOneStrip[iladder]->GetXaxis()->SetTitle("ADC");
 }//close for

   TH1D *hEnergyTwoStrip[NLadders];
 for(int iladder =0; iladder < NLadders; iladder++){
   hEnergyTwoStrip[iladder] = new TH1D(Form("Two Strip hEnergyCluster_iladder_%d",ladderID[iladder]), Form("Two Strips Energy clusters iladder %d",ladderID[iladder]),100,0,400);
   hEnergyTwoStrip[iladder]->GetXaxis()->SetTitle("ADC");
 }//close for

  TH1D *hEnergyThreeStrip[NLadders];
 for(int iladder =0; iladder < NLadders; iladder++){
   hEnergyThreeStrip[iladder] = new TH1D(Form("Three Strip hEnergyCluster_iladder_%d",ladderID[iladder]), Form("Three Strips Energy clusters iladder %d",ladderID[iladder]),100,0,400);
   hEnergyThreeStrip[iladder]->GetXaxis()->SetTitle("ADC");
 }//close for

  TH1D *hEnergyMainStrip[NLadders];
 for(int iladder =0; iladder < NLadders; iladder++){
   hEnergyMainStrip[iladder] = new TH1D(Form("Main Strip_iladder_%d",ladderID[iladder]), Form("Main Strip Energy iladder %d",ladderID[iladder]),50,0,250);
   hEnergyMainStrip[iladder]->GetXaxis()->SetTitle("ADC");
 }//close for

  TH1D *hEnergyOtherStrip[NLadders];
 for(int iladder =0; iladder < NLadders; iladder++){
   hEnergyOtherStrip[iladder] = new TH1D(Form("Other Strip_iladder_%d",ladderID[iladder]), Form("Other Strip Energy iladder %d",ladderID[iladder]),20,0,100);
   hEnergyOtherStrip[iladder]->GetXaxis()->SetTitle("ADC");
 }//close for

 //eta function
 TH1D *hEta[NLadders];
 for(int iladder =0; iladder <NLadders;iladder++){
   hEta[iladder] = new TH1D(Form("eta_ladder_%d",ladderID[iladder]), Form("Eta %d", ladderID[iladder]), 100,-1,1);
   hEta[iladder]->GetXaxis()->SetTitle("eta");
 }

 TH2D *hEta2D[NLadders];
 for(int iladder=0; iladder<NLadders;iladder++){
   hEta2D[iladder] = new TH2D(Form("Eta 2D_ladder_%d",ladderID[iladder]), Form("Eta 2d test ladder %d ",ladderID[iladder]), 100,-1,1,100,0,250);
   hEta2D[iladder]->GetXaxis()->SetTitle("Eta");
   hEta2D[iladder]->GetYaxis()->SetTitle("ADC");
 }//close i for
 
 TH2D *hZcorr[3];
 for(int i=0; i<3;i++){
   hZcorr[i] = new TH2D("Z_corr", "Z corr", 20,-0.5,20.5,20,-0.5,20.5);
   hZcorr[i]->GetXaxis()->SetTitle("Z");
   hZcorr[i]->GetYaxis()->SetTitle("Z");
 }//close i for
 
  TH2D *hZcorrLADD[8];
 for(int i=0; i<8;i++){
   hZcorrLADD[i] = new TH2D("Z_corr LADD", "Z corr LADD", 20,-0.5,20.5,20,-0.5,20.5);
   hZcorrLADD[i]->GetXaxis()->SetTitle("Z");
   hZcorrLADD[i]->GetYaxis()->SetTitle("Z");
 }//close i for
 
 //beam profile
 TH1D *occupancy[NLadders];
 for(int iladder =0; iladder <NLadders;iladder++){
   occupancy[iladder] = new TH1D(Form("occupancy_ladder_%d",ladderID[iladder]), Form("Occupancy %d", ladderID[iladder]), 93,0,93);
   occupancy[iladder]->GetXaxis()->SetTitle("cluster cog [mm]");
   
 }//close iladder for

 TH1D *hNstrips[NLadders];
 for(int iladder =0; iladder <NLadders;iladder++){
   hNstrips[iladder]= new TH1D(Form("hNstrips_ladder_%d",ladderID[iladder]), Form("# Strips in each cluster ladder %d",ladderID[iladder]),10,0.5,10.5);
   hNstrips[iladder]->GetXaxis()->SetTitle("# strips");
 }//close iladder for

 TH1D *hNclusters[NLadders];
 for(int iladder =0; iladder <NLadders;iladder++){
   hNclusters[iladder]= new TH1D(Form("hNclusters_ladder_%d",ladderID[iladder]), Form("# clusters per event in ladder %d",ladderID[iladder]),51,-0.5,50.5);
   hNclusters[iladder]->GetXaxis()->SetTitle("# clusters");
 }//close iladder for

 
 //2d test
 TH2D *h2test[NLadders];
 for(int iladder=0; iladder<NLadders;iladder++){
   h2test[iladder] = new TH2D(Form("h2dTest_ladder_%d",ladderID[iladder]), Form("2d test ladder %d ",ladderID[iladder]), 93, 0,93,250,0,500);
   h2test[iladder]->GetXaxis()->SetTitle("centroid [mm]");
   h2test[iladder]->GetYaxis()->SetTitle("ADC");
 }//close i for

 // EVENT LOOP //

 cout << "Processing percentage: " << "00" << "% \n";

 for(int ientry=0; ientry < t->GetEntries(); ientry++){
   t->GetEntry(ientry);

   for(int i=0; i<NLadders; i++)Nclust[i]=0;
      
   if ( (ientry%(NEvents/10+1))==int(NEvents/10) ) cout << "Processing percentage: " << 10.*int((ientry+1)/int(NEvents/10)) << "% \n";


   //CLUSTER LOOP//
   
   for(int i =0; i<4; i++)getcentroid[i]=0;
   //reset centroids


   double MaxCluster[20];
   double MinCluster[20];

    for(int i=0; i<20; i++){
     MaxCluster[i]=0;
   }

   for(int i=0; i<20; i++){
     MinCluster[i]=999999;
   }
   
   // select greatest cluster per ladder 
   for(int i=0; i<stkclusters->GetLast()+1; i++){
     DmpStkSiCluster* cluster = (DmpStkSiCluster*) stkclusters-> ConstructedAt(i);
    
     ladder = cluster->getLadderHardware();
     GetEnergy = cluster->getEnergy();
     for(int iladder = 0; iladder < NLadders; iladder++){
       if(ladder == ladderID[iladder] && GetEnergy > MaxCluster[iladder]){
   	 MaxCluster[iladder]=GetEnergy;
       }
       if(ladder == ladderID[iladder] && GetEnergy < MinCluster[iladder] && GetEnergy > 0){
	 MinCluster[iladder]=GetEnergy;
       }
     }
   }// end maxcluster


   //Fill Min e Max Cluster Histos
   for(int iladder = 0; iladder < NLadders; iladder++){
     // cout << MaxCluster[iladder] << endl;
     hEnergyClusterMin[iladder]->Fill(MinCluster[iladder]);
   }
   for(int iladder = 0; iladder < NLadders; iladder++){
     // cout << MaxCluster[iladder] << endl;
     hEnergyClusterMax[iladder]->Fill(sqrt(MaxCluster[iladder]/60.0));
   }
   /////

 
   //Z coorelations//to do: truncated mean instead of mean
   double cluster_mean_up=0;
   //Upstream Telescope
   for(int i=0; i<6; i++){cluster_mean_up+=MaxCluster[i];}
   cluster_mean_up=sqrt(cluster_mean_up/6/60.0);
   // cout << cluster_mean_up << endl;
   
   double cluster_mean_dut=0;
   //Upstream Telescope
   for(int i=6; i<14; i++){cluster_mean_dut+=MaxCluster[i];}
   cluster_mean_dut=sqrt(cluster_mean_dut/6/60.0);
   // cout << cluster_mean_dut << endl;
   
   double cluster_mean_dwn=0;
   //Downstream Telescope
   for(int i=14; i<20; i++){cluster_mean_dwn+=MaxCluster[i];}
   cluster_mean_dwn=sqrt(cluster_mean_dwn/6/60.0);
   //cout << cluster_mean_dwn << endl;

   //corr UP,DUT,DWN
   hZcorr[0]->Fill(cluster_mean_up,cluster_mean_dwn);
   hZcorr[1]->Fill(cluster_mean_up,cluster_mean_dut);
   hZcorr[2]->Fill(cluster_mean_dut,cluster_mean_dwn);

   
   //corr UP, Ladders
   for(int i=0; i<8; i++){
     hZcorrLADD[i]->Fill(sqrt(MaxCluster[i+6]/60),cluster_mean_up);
			 }

   
     
   for(int i=0; i<stkclusters->GetLast()+1; i++){
     DmpStkSiCluster* cluster = (DmpStkSiCluster*) stkclusters-> ConstructedAt(i);

     ladder = cluster->getLadderHardware();
     clustercog = cluster->getLadderCentroid();
     Nstrip = cluster->getNstrip();
     int FirstStrip = cluster->getFirstStrip();
     int LastStrip = cluster->getLastStrip();

     for(int iladder = 0; iladder < NLadders; iladder++){
       if(ladder == ladderID[iladder]){
	 Nclust[iladder]=Nclust[iladder]+1;
       }
     }
   }

    for(int iladder = 0; iladder < NLadders; iladder++){
     hNclusters[iladder]->Fill(Nclust[iladder]);
   }

   // //number of cluster selection
   // bool go=1;
   // int i=0;
   //   while(i < NLadders && go==1 ){
   //     if(Nclust[i]!=3) go=0;
   //     i=i+1;
   //   }

   //   if(go==0) continue;
     
   //   for(int iladder = 0; iladder < NLadders; iladder++){
   //     hNclusters[iladder]->Fill(Nclust[iladder]);
   //   }

     /////
	  
   for(int i=0; i<stkclusters->GetLast()+1; i++){
     DmpStkSiCluster* cluster = (DmpStkSiCluster*) stkclusters-> ConstructedAt(i);

     GetEnergy = cluster->getEnergy();
     
     ladder = cluster->getLadderHardware();
     clustercog = cluster->getLadderCentroid();
     Nstrip = cluster->getNstrip();
     int FirstStrip = cluster->getFirstStrip();
     int LastStrip = cluster->getLastStrip();

     if(GetEnergy>0){
     
     for(int iladder = 0; iladder < NLadders; iladder++){
       if(ladder == ladderID[iladder]){
	 
	 double MaxStripEnergy=0;
	 double SecondStripEnergy=0;
	 int PositionMax=0;
	 int PositionOth=0;
	 int ClusterStrip= FirstStrip+PositionMax;
	 double eta;
	 double totalEnergy;


	 //Eta (2 strip cluster only)
	 if(Nstrip==2){
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

	   hEnergyMainStrip[iladder]->Fill(MaxStripEnergy);
	   hEnergyOtherStrip[iladder]->Fill(SecondStripEnergy);

	   totalEnergy=(MaxStripEnergy+SecondStripEnergy);

	   if(PositionMax<PositionOth){
	     eta=(MaxStripEnergy-SecondStripEnergy)/(MaxStripEnergy+SecondStripEnergy);
	   } else {
	     eta=(SecondStripEnergy-MaxStripEnergy)/(MaxStripEnergy+SecondStripEnergy);
	   }

	  hEta[iladder]->Fill(eta);
	  hEta2D[iladder]->Fill(eta,totalEnergy);
	 }
	 //end Eta loop

	 hEnergyCluster[iladder]->Fill(GetEnergy);
	 occupancy[iladder]->Fill(clustercog);
	 hNstrips[iladder]->Fill(Nstrip);
	     if(Nstrip == 1){
	       hEnergyOneStrip[iladder]->Fill(GetEnergy);
	     } else if(Nstrip == 2){
	       hEnergyTwoStrip[iladder]->Fill(GetEnergy);
	     } else if(Nstrip == 3){
	       hEnergyThreeStrip[iladder]->Fill(GetEnergy);
	     }
	 h2test[iladder]->Fill(clustercog,GetEnergy);
       }//close ladder if 
     }//close ladder for
     }
   }//close cluster loop

    if (ientry == int(NEvents-1)) cout << "Processing percentage: " << "100" << "% \n";
  
 }//close event loop
 
 // TCanvas *c = new TCanvas("c", "Cluster Energy", 1920,1080);
 // // gStyle->SetOptStat(0);
 // c->Divide(5,4);
 // for(int i=0; i< NLadders;i++){
 //   c->cd(i+1);
 //   hEnergyCluster[i]->SetLineColor(kRed);
 //   hEnergyCluster[i]->SetLineWidth(2);
 //   hEnergyCluster[i]->Draw();
   
 //   // hEnergyOneStrip[i]->SetLineColor(kBlue);
 //   // hEnergyOneStrip[i]->Draw("SAME");
   
 //   // hEnergyTwoStrip[i]->SetLineColor(kGreen);
 //   // hEnergyTwoStrip[i]->Draw("SAME");


 //   // hEnergyThreeStrip[i]->SetLineColor(kBlack);
 //   // hEnergyThreeStrip[i]->Draw("SAME");

 //   // gPad->BuildLegend();

 //  hEnergyCluster[i]->Write();
 //  // hEnergyOneStrip[i]->Write();
 //  // hEnergyTwoStrip[i]->Write();
 //  // hEnergyThreeStrip[i]->Write();
 // }//close for
 // c->SaveAs(plotdirectory+"/henclust_"+run_string+".png");
 
 // TCanvas *ca = new TCanvas("ca", "Cluster Energy Min", 1920,1080);
 // // gStyle->SetOptStat(0);
 // ca->Divide(5,4);
 // for(int i=0; i< NLadders;i++){
 //   ca->cd(i+1);
 //   hEnergyClusterMin[i]->SetLineColor(kRed);
 //   hEnergyClusterMin[i]->SetLineWidth(2);
 //   hEnergyClusterMin[i]->Draw();
 
 //  hEnergyClusterMin[i]->Write();

 // }//close for
 // ca->SaveAs(plotdirectory+"/henclust_min_"+run_string+".png");  

 // TCanvas *cb = new TCanvas("cb", "Cluster Energy Max", 1920,1080);
 // // gStyle->SetOptStat(0);
 // cb->Divide(5,4);
 // for(int i=0; i< NLadders;i++){
 //   cb->cd(i+1);
 //   hEnergyClusterMax[i]->SetLineColor(kRed);
 //   hEnergyClusterMax[i]->SetLineWidth(2);
 //   hEnergyClusterMax[i]->Draw();
 
 //  hEnergyClusterMax[i]->Write();

 // }//close for
 // cb->SaveAs(plotdirectory+"/henclust_max_"+run_string+".png");


 // TCanvas *c1 = new TCanvas("c1", "# strips ", 1920,1080);
 // c1->Divide(5,4);
 // for(int i=0; i< NLadders;i++){
 //   c1->cd(i+1);
 //   hNstrips[i]->Draw();
 //   hNstrips[i]->Write();
 // }
 // c1->SaveAs(plotdirectory+"/n_strips_"+run_string+".png");

 // TCanvas *c2 = new TCanvas("c2", "Occupancy ", 1920,1080);
 // c2->Divide(5,4);
 // for(int i=0; i< NLadders;i++){
 //   c2->cd(i+1);
 //   occupancy[i]->Draw();
 //   occupancy[i]->Write();
 // }//close for
 // c2->SaveAs(plotdirectory+"/occupancy_"+run_string+".png");
  
 //  TCanvas *c3 = new TCanvas("c3", "2d test ", 1920,1080);
 // c3->Divide(5,4);
 // for(int i=0; i< NLadders;i++){
 //   c3->cd(i+1);
 //   gPad->SetLogz();
 //   h2test[i]->Draw("colz");
 //   if(ladderID[i]==126){
 //   h2test[i]->GetXaxis()->SetRange(0,78);
 //   }
 //   h2test[i]->Write();
 // }//close for
 // c3->SaveAs(plotdirectory+"/2dTest_"+run_string+".png");

 
 // TCanvas *c5 = new TCanvas("c5", "Main Strip Energy",1920,1080);
 // c5->Divide(5,4);
 // for(int i=0; i<NLadders;i++){
 //   c5->cd(i+1);
 //   hEnergyMainStrip[i]->Draw();
 //   hEnergyMainStrip[i]->Write();
 // }//close for 
 // c5->SaveAs(plotdirectory+"/mainstrip_"+run_string+".png");
  
 // TCanvas *c6 = new TCanvas("c6", "Other Strip Energy", 1920,1080);
 // c6->Divide(5,4);
 // for(int i=0; i<NLadders;i++){
 //   c6->cd(i+1);
 //   hEnergyOtherStrip[i]->Draw();
 //   hEnergyOtherStrip[i]->Write();
 // }//close for
 // c6->SaveAs(plotdirectory+"/otherstrip_"+run_string+".png");

 // TCanvas *c7 = new TCanvas("c7", "Eta function", 1920,1080);
 // c7->Divide(5,4);
 // for(int i=0; i<NLadders;i++){
 //   c7->cd(i+1);
 //   hEta[i]->Draw();
 //   hEta[i]->Write();
 // }//close for
 // c7->SaveAs(plotdirectory+"/eta_"+run_string+".png");

 //  TCanvas *c8 = new TCanvas("c8", "Eta 2D Test", 1920,1080);
 // c8->Divide(5,4);
 // for(int i=0; i< NLadders;i++){
 //   c8->cd(i+1);
 //   gPad->SetLogz();
 //   hEta2D[i]->Draw("colz");
 //   hEta2D[i]->Write();
 // }//close for
 // c8->SaveAs(plotdirectory+"/eta2d_"+run_string+".png");

 // TCanvas *c8b= new TCanvas("c8b", "Eta 2D Test Profile", 1920,1080);
 // c8b->Divide(5,4);
 // for(int i=0; i< NLadders;i++){
 //   c8b->cd(i+1);
 //   gPad->SetLogz();
 //   TProfile *prof = hEta2D[i]->ProfileX();
 //   prof->Draw(); 
 //   prof->Write();
 // }//close for
 // c8b->SaveAs(plotdirectory+"/eta_profile_"+run_string+".png");

 
 // TCanvas *c9 = new TCanvas("c9", "# clusters ",1920 ,1080);
 // c9->Divide(5,4);
 // for(int i=0; i< NLadders;i++){
 //   c9->cd(i+1);
 //   hNclusters[i]->Draw();
 //   hNclusters[i]->Write();
 // }
 // c9->SaveAs(plotdirectory+"/n_clusters_"+run_string+".png");

 TCanvas *c10 = new TCanvas("c10", "Z corr ",1920 ,1080);
 c10->Divide(3,1);
 for(int i=0; i< 3;i++){
   c10->cd(i+1);
   gPad->SetLogz();
   hZcorr[i]->Draw("colz");
   hZcorr[i]->Write();
 }
 c10->SaveAs(plotdirectory+"/z_corr_"+run_string+".png");

  TCanvas *c11 = new TCanvas("c11", "Z corr LADD",1920 ,1080);
 c11->Divide(4,2);
 for(int i=0; i< 8;i++){
   c11->cd(i+1);
   gPad->SetLogz();
   hZcorrLADD[i]->Draw("colz");
   hZcorrLADD[i]->Write();
 }
 c11->SaveAs(plotdirectory+"/z_corr_LADD_"+run_string+".png");
}//close macro
