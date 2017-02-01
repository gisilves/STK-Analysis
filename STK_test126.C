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

void STK_test126(char filename[200], char runnumber[200]){
  gSystem->Load("DmpSoftware/trunk/Event/libDmpEvent.so");

  char prefix[100]="Results/out_126_";
  char suffix[100]=".root";
  
  TFile* f = new TFile(filename,"READ");
  TFile* out = new TFile(strcat(strcat(prefix,runnumber),suffix), "recreate");
  
  TTree* t =(TTree*)f->Get("CollectionTree");
  cout << endl;
  f->ls();
  cout << endl;
  
  //Register STK collections//
  TClonesArray* stkclusters = new TClonesArray("DmpStkSiCluster");
  t->SetBranchAddress("StkClusterCollection",&stkclusters);

  //Register SCINTILLATOR collections
  TClonesArray* scintillators = new TClonesArray("Scintillator");  
  t->SetBranchAddress("ScintillatorAdcCollection",&scintillators);
    
 int NEvents  = t->GetEntries();
 cout << "Number of events: " << NEvents << endl;
 cout << endl;
 

 const char *ladderName[NLadders]={"LPFM905T"};
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
 TH1D *hEnergyCluster = new TH1D("hEnergyCluster_iladder 126","Energy Clusters 126",50,0,200);
   hEnergyCluster->GetXaxis()->SetTitle("ADC");
 TH1D *hEnergyClusterhp = new TH1D("hEnergyCluster_iladder 126hp","Energy Clusters 126 Half Pitch",50,0,200);
   hEnergyClusterhp->GetXaxis()->SetTitle("ADC");
 TH1D *hEnergyCluster8si = new TH1D("hEnergyCluster_iladder 1268si","Energy Clusters  126 8 Si",50,0,200);
   hEnergyCluster8si->GetXaxis()->SetTitle("ADC");
   
 TH1D *hEnergyOneStrip= new TH1D("One Strip hEnergyCluster_iladder_126","One Strip Energy clusters iladder 126",50,0,200);
   hEnergyOneStrip->GetXaxis()->SetTitle("ADC");
 TH1D *hEnergyOneStriphp= new TH1D("One Strip hEnergyCluster_iladder_126 hp","One Strip Energy clusters iladder 126 hp",50,0,200);
   hEnergyOneStriphp->GetXaxis()->SetTitle("ADC");
 TH1D *hEnergyOneStrip8si= new TH1D("One Strip hEnergyCluster_iladder_126 8si","One Strip Energy clusters iladder 126 8si",50,0,200);
   hEnergyOneStrip8si->GetXaxis()->SetTitle("ADC");
   
 TH1D *hEnergyTwoStrip= new TH1D("Two Strip hEnergyCluster_iladder_126","Two Strip Energy clusters iladder 126",50,0,200);
   hEnergyTwoStrip->GetXaxis()->SetTitle("ADC");
 TH1D *hEnergyTwoStriphp= new TH1D("Two Strip hEnergyCluster_iladder_126hp","Two Strip Energy clusters iladder 126 Half Pitch",50,0,200);
   hEnergyTwoStriphp->GetXaxis()->SetTitle("ADC");
 TH1D *hEnergyTwoStrip8si= new TH1D("Two Strip hEnergyCluster_iladder_126 8si","Two Strip Energy clusters iladder 126 8si",50,0,200);
   hEnergyTwoStrip8si->GetXaxis()->SetTitle("ADC");
   
 TH1D *hEnergyThreeStrip= new TH1D("Three Strip hEnergyCluster_iladder_126","Three Strip Energy clusters iladder 126",50,0,200);
   hEnergyTwoStrip->GetXaxis()->SetTitle("ADC");
 TH1D *hEnergyThreeStriphp= new TH1D("Three Strip hEnergyCluster_iladder_126hp","Three Strip Energy clusters iladder 126 Half Pitch",50,0,200);
   hEnergyThreeStriphp->GetXaxis()->SetTitle("ADC");
 TH1D *hEnergyThreeStrip8si= new TH1D("Three Strip hEnergyCluster_iladder_126 8si","Three Strip Energy clusters iladder 126 8si",50,0,200);
   hEnergyThreeStrip8si->GetXaxis()->SetTitle("ADC");

   
  TH1D *hEnergyMainStrip = new TH1D("Main Strip_iladder_126","Main Strip Energy iladder 126",50,0,250);
   hEnergyMainStrip->GetXaxis()->SetTitle("ADC");
  TH1D *hEnergyMainStriphp = new TH1D("Main Strip_iladder_126hp","Main Strip Energy iladder 126 Half Pitch",50,0,250);
   hEnergyMainStriphp->GetXaxis()->SetTitle("ADC");
  TH1D *hEnergyMainStrip8si = new TH1D("Main Strip_iladder_1268si","Main Strip Energy iladder 126 8 Si",50,0,250);
   hEnergyMainStrip8si->GetXaxis()->SetTitle("ADC");

  TH1D *hEnergyOtherStrip = new TH1D("Other Strip_iladder_126","Other Strip Energy iladder 126",50,0,250);
   hEnergyOtherStrip->GetXaxis()->SetTitle("ADC");
  TH1D *hEnergyOtherStriphp = new TH1D("Other Strip_iladder_126hp","Other Strip Energy iladder 126 Half Pitch",50,0,250);
   hEnergyOtherStriphp->GetXaxis()->SetTitle("ADC");
  TH1D *hEnergyOtherStrip8si = new TH1D("Other Strip_iladder_1268si","Other Strip Energy iladder 126 8 Si",50,0,250);
   hEnergyOtherStrip8si->GetXaxis()->SetTitle("ADC");   

 //eta function
  TH1D *hEta = new TH1D("126 eta","Eta ladder 126", 50,0,1);
    hEta->GetXaxis()->SetTitle("eta");
  TH1D *hEtahp = new TH1D("126 etahp","Eta ladder 126 Half Pitch", 50,0,1);
    hEtahp->GetXaxis()->SetTitle("eta");
  TH1D *hEta8si = new TH1D("126 eta8si","Eta ladder 126 8 Si", 50,0,1);
    hEta8si->GetXaxis()->SetTitle("eta"); 

 TH2D *hEta2D = new TH2D("Eta 2D_ladder_126", "Eta 2d test ladder 126", 50,0,1,100,0,250);
   hEta2D->GetXaxis()->SetTitle("Eta");
   hEta2D->GetYaxis()->SetTitle("ADC");
 TH2D *hEta2Dhp = new TH2D("Eta 2D_ladder_126hp", "Eta 2d test ladder 126 Half Pitch", 50,0,1,100,0,250);
   hEta2Dhp->GetXaxis()->SetTitle("Eta");
   hEta2Dhp->GetYaxis()->SetTitle("ADC");
 TH2D *hEta2D8si = new TH2D("Eta 2D_ladder_1268si", "Eta 2d test ladder 126 8 Si", 50,0,1,100,0,250);
   hEta2D8si->GetXaxis()->SetTitle("Eta");
   hEta2D8si->GetYaxis()->SetTitle("ADC");

   
 //beam profile
 TH1D *occupancy = new TH1D("occupancy_ladder_126","Occupancy 126", 93,0,93);
   occupancy->GetXaxis()->SetTitle("cluster cog [mm]");
 TH1D *occupancyhp = new TH1D("occupancy_ladder_126hp","Occupancy 126 Half Pitch", 93,0,93);
   occupancyhp->GetXaxis()->SetTitle("cluster cog [mm]");
 TH1D *occupancy8si = new TH1D("occupancy_ladder_1268si","Occupancy 126 8 Si", 93,0,93);
   occupancy8si->GetXaxis()->SetTitle("cluster cog [mm]");

   
 // number of strips per cluster
 TH1D *hNstrips = new TH1D("hNstrips_ladder_126","# Strips in each cluster ladder 126",10,0.5,10.5);
   hNstrips->GetXaxis()->SetTitle("# strips");
 TH1D *hNstripshp = new TH1D("hNstrips_ladder_126hp","# Strips in each cluster ladder 126 Half Pitch",10,0.5,10.5);
   hNstripshp->GetXaxis()->SetTitle("# strips");
 TH1D *hNstrips8si = new TH1D("hNstrips_ladder_1268si","# Strips in each cluster ladder 126 8 Si",10,0.5,10.5);
   hNstrips8si->GetXaxis()->SetTitle("# strips");


 // //2d test
 TH2D *h2test = new TH2D("h2dTest_ladder_126","2d test ladder 126", 93, 0,93,250,0,500);
   h2test->GetXaxis()->SetTitle("centroid [mm]");
   h2test->GetYaxis()->SetTitle("ADC");
 TH2D *h2testhp = new TH2D("h2dTest_ladder_126hp","2d test ladder 126 Half Pitch", 93, 0,93,250,0,500);
   h2testhp->GetXaxis()->SetTitle("centroid [mm]");
   h2testhp->GetYaxis()->SetTitle("ADC");
 TH2D *h2test8si = new TH2D("h2dTest_ladder_1268si","2d test ladder 126 8 Si", 93, 0,93,250,0,500);
   h2test8si->GetXaxis()->SetTitle("centroid [mm]");
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
     
     //find the cluster with highest energy
     if(GetEnergy>MaxEne){
       MaxEne= GetEnergy;
       iscluster =1;
     }//close if
     else iscluster=0;
     

     ladder = cluster->getLadderHardware();
     clustercog = cluster->getLadderCentroid();
     Nstrip = cluster->getNstrip();
     int FirstStrip = cluster->getFirstStrip();
     int LastStrip = cluster->getLastStrip();
  
       if(ladder == 126){
	 double MaxStripEnergy=0;
	 double SecondStripEnergy=0;
	 int PositionMax=0;
	 int PositionOth=0;
	 int ClusterStrip= FirstStrip+PositionMax;
	 double eta;
	 double totalEnergy;
	 if(isBad(badChan,NbadChan,ClusterStrip)==0){  
        
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
	   if(ClusterStrip>160 && ClusterStrip<193){
	     hEnergyMainStriphp->Fill(MaxStripEnergy);
	     } else if(ClusterStrip>192 && ClusterStrip<257){
	     hEnergyMainStrip8si->Fill(MaxStripEnergy);
	     } else if(ClusterStrip<161 || ClusterStrip>256){
	     hEnergyMainStrip->Fill(MaxStripEnergy);
	   }

	   if(ClusterStrip>160 && ClusterStrip<193){
	     hEnergyOtherStriphp->Fill(MaxStripEnergy);
	     } else if(ClusterStrip>192 && ClusterStrip<257){
	     hEnergyOtherStrip8si->Fill(MaxStripEnergy);
	     } else if(ClusterStrip<161 || ClusterStrip>256){
	     hEnergyOtherStrip->Fill(MaxStripEnergy);
	   }

	   totalEnergy=(MaxStripEnergy+SecondStripEnergy);

	   if(PositionMax<PositionOth){
	     eta=(SecondStripEnergy)/(MaxStripEnergy+SecondStripEnergy);
	   } else {
	     eta=(MaxStripEnergy)/(MaxStripEnergy+SecondStripEnergy);
	   }

	   
	   if(ClusterStrip>160 && ClusterStrip<193){
	       hEtahp->Fill(eta);
	     } else if(ClusterStrip>192 && ClusterStrip<257){
	       hEta8si->Fill(eta);
	     } else if(ClusterStrip<161 || ClusterStrip>256){
	       hEta->Fill(eta);
	   }
	   
	  if(ClusterStrip>160 && ClusterStrip<193){
	    hEta2Dhp->Fill(eta,totalEnergy);
	     } else if(ClusterStrip>192 && ClusterStrip<257){
	    hEta2D8si->Fill(eta,totalEnergy);
	     } else if(ClusterStrip<161 || ClusterStrip>256){
	    hEta2D->Fill(eta,totalEnergy);
	  }
	 
	 }
	 //end Eta loop


	 //fill energy
	  if(ClusterStrip>160 && ClusterStrip<193){
	       hEnergyClusterhp->Fill(GetEnergy);
	     } else if(ClusterStrip>192 && ClusterStrip<257){
	       hEnergyCluster8si->Fill(GetEnergy);
	     } else if(ClusterStrip<161 || ClusterStrip>256){
		   hEnergyCluster->Fill(GetEnergy);
		 }
	  //end fill energy

	 if(ClusterStrip>160 && ClusterStrip<193){
			occupancyhp->Fill(clustercog);
	     } else if(ClusterStrip>192 && ClusterStrip<257){
			occupancy8si->Fill(clustercog);
	 } else if(ClusterStrip<161 || ClusterStrip>256){
			occupancy->Fill(clustercog);
	 }

	 //number of strip per cluster
	  if(ClusterStrip>160 && ClusterStrip<193){
			hNstripshp->Fill(Nstrip);
	     } else if(ClusterStrip>192 && ClusterStrip<257){
			hNstrips8si->Fill(Nstrip);
	     } else if(ClusterStrip<161 || ClusterStrip>256){
			hNstrips->Fill(Nstrip);
		}
	 //number of strip end

	 

	   if(ClusterStrip>160 && ClusterStrip<193){
	     if(Nstrip == 1){
	       hEnergyOneStriphp->Fill(GetEnergy);
	     } else if(Nstrip == 2){
	       hEnergyTwoStriphp->Fill(GetEnergy);
	     } else if(Nstrip == 3){
	       hEnergyThreeStriphp->Fill(GetEnergy);
	     }
	   } else if(ClusterStrip>192 && ClusterStrip<257){
	      if(Nstrip == 1){
	       hEnergyOneStrip8si->Fill(GetEnergy);
	     } else if(Nstrip == 2){
	       hEnergyTwoStrip8si->Fill(GetEnergy);
	     } else if(Nstrip == 3){
	       hEnergyThreeStrip8si->Fill(GetEnergy);
	     }
	   } else if(ClusterStrip<161 || ClusterStrip>256){
	   if(Nstrip == 1){
	       hEnergyOneStrip->Fill(GetEnergy);
	     } else if(Nstrip == 2){
	       hEnergyTwoStrip->Fill(GetEnergy);
	     } else if(Nstrip == 3){
	       hEnergyThreeStrip->Fill(GetEnergy);
	     }
		}


	   if(ClusterStrip>160 && ClusterStrip<193){
	     h2testhp->Fill(clustercog,GetEnergy);
	     } else if(ClusterStrip>192 && ClusterStrip<257){
	     h2test8si->Fill(clustercog,GetEnergy);
	   } else if(ClusterStrip<161 || ClusterStrip>256){
	   h2test->Fill(clustercog,GetEnergy);
		}
	   }
       }//close ladder if 

     
   }//close cluster loop

 
    if (ientry == int(NEvents-1)) cout << "Processing percentage: " << "100" << "% \n";
  
 }//close event loop
 
  // TCanvas *c = new TCanvas("c", "Cluster Energy", 900,1000);
  // gStyle->SetOptStat(0);
  //   hEnergyCluster->SetLineColor(kRed);
  //   hEnergyCluster->SetLineWidth(2);
  //   hEnergyCluster->Draw();
  //   hEnergyOneStrip->SetLineColor(kBlue);
  //   hEnergyOneStrip->Draw("SAME");
  //   hEnergyTwoStrip->SetLineColor(kGreen);
  //   hEnergyTwoStrip->Draw("SAME");
  //   hEnergyThreeStrip->SetLineColor(kBlack);
  //   hEnergyThreeStrip->Draw("SAME");
  //   gPad->BuildLegend();
  //  hEnergyCluster->Write();
  //  hEnergyOneStrip->Write();
  //  hEnergyTwoStrip->Write();
  //  hEnergyThreeStrip->Write();

  // TCanvas *c1 = new TCanvas("c1", "Cluster Energy Half Pitch", 900,1000);
  // gStyle->SetOptStat(0);
  //   hEnergyClusterhp->SetLineColor(kRed);
  //   hEnergyClusterhp->SetLineWidth(2);
  //   hEnergyClusterhp->Draw();
  //   hEnergyOneStriphp->SetLineColor(kBlue);
  //   hEnergyOneStriphp->Draw("SAME");
  //   hEnergyTwoStriphp->SetLineColor(kGreen);
  //   hEnergyTwoStriphp->Draw("SAME");
  //   hEnergyThreeStriphp->SetLineColor(kBlack);
  //   hEnergyThreeStriphp->Draw("SAME");
  //   gPad->BuildLegend();
  //  hEnergyClusterhp->Write();
  //  hEnergyOneStriphp->Write();
  //  hEnergyTwoStriphp->Write();
  //  hEnergyThreeStriphp->Write();

  //  TCanvas *c2 = new TCanvas("c2", "Cluster Energy 8 Si", 900,1000);
  // gStyle->SetOptStat(0);
  //   hEnergyCluster8si->SetLineColor(kRed);
  //   hEnergyCluster8si->SetLineWidth(2);
  //   hEnergyCluster8si->Draw();
  //   hEnergyOneStrip8si->SetLineColor(kBlue);
  //   hEnergyOneStrip8si->Draw("SAME");
  //   hEnergyTwoStrip8si->SetLineColor(kGreen);
  //   hEnergyTwoStrip8si->Draw("SAME");
  //   hEnergyThreeStrip8si->SetLineColor(kBlack);
  //   hEnergyThreeStrip8si->Draw("SAME");
  //   gPad->BuildLegend();
  //   hEnergyCluster8si->Write();
  //   hEnergyOneStrip8si->Write();
  //   hEnergyTwoStrip8si->Write();
  //   hEnergyThreeStrip8si->Write(); 

//    TCanvas *c24 = new TCanvas("c24", "Cluster Energy per zone", 900,1000);
//    	gStyle->SetOptStat(0);
//     hEnergyCluster->SetLineColor(kRed);
//     hEnergyCluster->Draw();
//     hEnergyClusterhp->SetLineColor(kBlue);
//     hEnergyClusterhp->Draw("SAME");
//     hEnergyCluster8si->SetLineColor(kGreen);
//     hEnergyCluster8si->Draw("SAME");
//     gPad->BuildLegend();
//     hEnergyCluster->Write();
//     hEnergyClusterhp->Write();
//     hEnergyCluster8si->Write();


// TCanvas *c25 = new TCanvas("c25", "# strip per zone", 900,1000);
// gStyle->SetOptStat(0);
// hNstrips->SetLineColor(kRed);
// hNstrips->Draw();
// hNstripshp->SetLineColor(kBlue);
// hNstripshp->Draw("SAME");
// hNstrips8si->SetLineColor(kGreen);
// hNstrips8si->Draw("SAME");
// gPad->BuildLegend();
// hNstrips->Write();
// hNstripshp->Write();
// hNstrips8si->Write();


//  gStyle->SetOptStat(1111111);
    
//  TCanvas *c3 = new TCanvas("c3", "# strips ", 900,1000);
//    hNstrips->Draw();
//    hNstrips->Write();
//  TCanvas *c4 = new TCanvas("c4", "# strips Half Pitch", 900,1000);
//    hNstripshp->Draw();
//    hNstripshp->Write();
//   TCanvas *c5 = new TCanvas("c5", "# strips 8 Si", 900,1000);
//    hNstrips8si->Draw();
//    hNstrips8si->Write(); 

 // TCanvas *c6 = new TCanvas("c6", "Eta function", 900,1000);
 // hEta->Draw();
 // hEta->Write();
 // TCanvas *c7 = new TCanvas("c7", "Eta function Half Pitch", 900,1000);
 // hEtahp->Draw();
 // hEtahp->Write();
 // TCanvas *c8 = new TCanvas("c8", "Eta function 8 Si", 900,1000);
 // hEta8si->Draw();
 // hEta8si->Write();

 
  // TCanvas *c9 = new TCanvas("c9", "Eta 2D Test", 900,1000);
  // gPad->SetLogz();
  // hEta2D->Draw("colz");
  // hEta2D->Write();
  // TCanvas *c9b = new TCanvas("c9b", "Eta 2D Test Profile", 900,1000);
  // gPad->SetLogz();
  // TProfile *prof = hEta2D->ProfileX();
  // prof->Draw();
  // TCanvas *c10 = new TCanvas("c10", "Eta 2D Test Half Pitch", 900,1000);
  // gPad->SetLogz();
  // hEta2Dhp->Draw("colz");
  // hEta2Dhp->Write();
TCanvas *c10b = new TCanvas("c10b", "Eta 2D Test Half Pitch Profile", 900,1000);
gPad->SetLogz();
TProfile *profhp = hEta2Dhp->ProfileX();
profhp->Draw();
  // TCanvas *c11 = new TCanvas("c11", "Eta 2D Test 8 Si", 900,1000);
  // gPad->SetLogz();
  // hEta2D8si->Draw("colz");
  // hEta2D8si->Write();
  // TCanvas *c11b = new TCanvas("c11b", "Eta 2D Test 8 Si Profile", 900,1000);
  // gPad->SetLogz();
  // TProfile *prof8si = hEta2D8si->ProfileX();
  // prof8si->Draw();
  
  
 // TCanvas *c12 = new TCanvas("c12", "Occupancy ", 900,1000);
 //   occupancy->Draw();
 //   occupancy->Write();
 // TCanvas *c13 = new TCanvas("c13", "Occupancy ", 900,1000);
 //   occupancyhp->Draw();
 //   occupancyhp->Write();
 // TCanvas *c14 = new TCanvas("c14", "Occupancy ", 900,1000);
 //   occupancy8si->Draw();
 //   occupancy8si->Write();

   
  // TCanvas *c15 = new TCanvas("c15", "2d test ", 900,1000);
  //  gPad->SetLogz();
  //  h2test->Draw("colz");
  //  h2test->Write();
  TCanvas *c16 = new TCanvas("c16", "2d test ", 900,1000);
   gPad->SetLogz();
   h2testhp->Draw("colz");
   h2testhp->Write();
  // TCanvas *c17 = new TCanvas("c17", "2d test ", 900,1000);
  //  gPad->SetLogz();
  //  h2test8si->Draw("colz");
  //  h2test8si->Write();

   
 // TCanvas *c18 = new TCanvas("c18", "Main Strip Energy", 900,1000);
 //  hEnergyMainStrip->Draw();
 //  hEnergyMainStrip->Write();
 // TCanvas *c19 = new TCanvas("c19", "Main Strip Energy", 900,1000);
 //  hEnergyMainStriphp->Draw();
 //  hEnergyMainStriphp->Write();
 // TCanvas *c20 = new TCanvas("c20", "Main Strip Energy", 900,1000);
 //  hEnergyMainStrip8si->Draw();
 //  hEnergyMainStrip8si->Write();
  

 // TCanvas *c21 = new TCanvas("c21", "Other Strip Energy", 900,1000);
 //  hEnergyOtherStrip->Draw();
 //  hEnergyOtherStrip->Write();
 // TCanvas *c22 = new TCanvas("c22", "Other Strip Energy", 900,1000);
 //  hEnergyOtherStriphp->Draw();
 //  hEnergyOtherStriphp->Write();
 //  TCanvas *c23 = new TCanvas("c23", "Other Strip Energy", 900,1000);
 //  hEnergyOtherStrip8si->Draw();
 //  hEnergyOtherStrip8si->Write();
 
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

