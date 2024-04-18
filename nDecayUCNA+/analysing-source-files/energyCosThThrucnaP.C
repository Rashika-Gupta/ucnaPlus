/*
Date : March 20, 2023
Details : converting root to txt file
*/

#include<fstream>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include "TRandom3.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TLeaf.h"
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TChain.h>

using namespace std;

void energyCosThThrucnaP(){
    Char_t temp[200];
    std::string filename ; 

    ofstream file, file1, file2,file01, file02;
    
    std::cout<<"opening root file "<<endl;
// declare variables to stoer thevalues of branches
    Double_t dEeSilicon1, dEeSilicon2, Te0, timeHit1, timeHit2, pz0_e, pZe, pOutFoilz,pInFoil1z,dEeFoil1,dEeFoil2;
    Double_t nBins = 160;
    Double_t keMax = 800;
    Double_t binWidth = keMax/nBins;
//saving the evnets that wouold also be emiited in opposite direction. 
    Int_t TotalNoHits;
    Double_t counter = 0;
    TChain chain("Tout");
 //   Double_t totEnergy=[];
 //   Double_t totcosTh=[];
 //   Double_t nTot=[];
//////v///
 //   chain.Add("/mnt/gpfs3_amd/scratch/rgu245/Now/ucnaPlus/nDecayUCNA+/root-files/ucna+/event-gen-mar-27-90mil/6.*.root");
    chain.Add("/mnt/gpfs3_amd/scratch/rgu245/Now/ucnaPlus/nDecayUCNA+/root-files/ucna+/asym-set-0/6.*.root");
   
    chain.SetBranchAddress("Te0", &Te0);
    chain.SetBranchAddress("dEeSilicon1", &dEeSilicon1);
    chain.SetBranchAddress("dEeSilicon2", &dEeSilicon2);
    
    chain.SetBranchAddress("dEeFoil1", &dEeFoil1);
    chain.SetBranchAddress("dEeFoil2", &dEeFoil2);
    
    chain.SetBranchAddress("timeHit1", &timeHit1);
    chain.SetBranchAddress("timeHit2", &timeHit2);
    chain.SetBranchAddress("pz0_e", &pZe);
    
    Double_t E_thr = 10;
    file01.open("/mnt/gpfs3_amd/scratch/rgu245/Now/ucnaPlus/nDecayUCNA+/type01_energycosth_foil_a=0.txt"); 
    file02.open("/mnt/gpfs3_amd/scratch/rgu245/Now/ucnaPlus/nDecayUCNA+/type02_energycosth_foil_a=0.txt");
    file1.open("/mnt/gpfs3_amd/scratch/rgu245/Now/ucnaPlus/nDecayUCNA+/dead_events_foil_a=0.txt");
    file2.open("/mnt/gpfs3_amd/scratch/rgu245/Now/ucnaPlus/nDecayUCNA+/dead_events_belowThreshold_foil_a=0.txt");
  
    std::cout <<"EThr : "<<E_thr;
    for (Long64_t i = 0; i < 80000000; i++) {
        chain.GetEntry(i);
        Double_t E1 = dEeSilicon1;
        Double_t E2 = dEeSilicon2;
        Double_t KE = Te0;
        Double_t eFoil1 = dEeFoil1 ;
        Double_t eFoil2 = dEeFoil2 ;  

        if(E1 < E_thr){
            Double_t x  = 0;
            E1 = x;
        }
        
        if(E2 < E_thr){
      
            Double_t y  = 0;
            E2 = y;
        }
        Double_t eTotReconst = E1 + E2;

      /*filling type 0 above threshold */
        /* deposits energy in dtector 1*/
        if((E1 !=0 && E2 == 0) ){
            file01<<KE<<" "<<pZe<<" "<<E1<<" "<<dEeSilicon2<<" "<<eFoil1<<" "<<eFoil2<<endl;  
            }
        if((E2 !=0 && E1 == 0) ){
            file02<<KE<<" "<<pZe<<" "<<E2<<" "<<dEeSilicon1<<" "<<eFoil1<<" "<<eFoil2<<endl;  
          
        }
// dead events 
        if(E1 == 0 && E2 == 0){
	       file1<<KE<<" "<<pZe<<" "<<eFoil1<<" "<<eFoil2<<endl;
           if(dEeSilicon1 != 0 || dEeSilicon2!=0){
            file2<<KE<<" "<<dEeSilicon1<<" "<<dEeSilicon2<<" "<<pZe<<endl;
           }
        }   
    }
     file1.close();
     file2.close();
    file02.close();
   
    file01.close();
   //fout->Write();
   //fout->Close();
}