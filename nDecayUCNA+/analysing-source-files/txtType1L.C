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

void txtType1L(){
    Char_t temp[200];
    std::string filename ; 

    ofstream file, file1, file2;
     file1.open("/mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/nDecay+_100Mill_01_oppDirection.txt");
     file2.open("/mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/nDecay+_100Mill_02_oppDirection.txt");		
     file.open("/mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/nDecay_100Mill_foilBothSide.txt");
 sprintf(temp, "/mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/nDecay_100Mill_foilBothSide.root");

// declare variables to stoer thevalues of branches
    Double_t dEeSilicon1, dEeSilicon2, Te0, timeHit1, timeHit2, pz0_e, pZe, pOutFoilz,pInFoil1z;
    TFile *fout = new TFile(temp,"RECREATE");  //creating output file
    
    TH1D *KE = new TH1D("KE" ,"KE", 8000., 0., 800.); //creating histogrms
    TH1D *EdepType01 = new TH1D("EdepType01" ,"EdepType01", 8000., 0., 800.); //creating histogrms
    TH1D *EdepType11 = new TH1D("EdepType11" ,"EdepType11", 8000., 0., 800.);
    TH1D *EdepType02 = new TH1D("EdepType02" ,"EdepType02", 8000., 0., 800.); //creating histogrms
    TH1D *EdepType12 = new TH1D("EdepType12" ,"EdepType12", 8000., 0., 800.);
    TH1D *EdepType0  = new TH1D("EdepType0" ,"EdepType0", 8000., 0., 800.); //creating histogrms
    TH1D *EdepType1  = new TH1D("EdepType1" ,"EdepType1", 8000., 0., 800.);
//saving the evnets that wouold also be emiited in opposite direction. 
    Int_t TotalNoHits;
    Double_t counter = 0;
    TChain chain("Tout");
    chain.Add("/mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/UCNA+/100Mill/6.*.root");
 //   
 //   chain.Add("/mnt/gpfs3_amd/scratch/rgu245/Now/foil-mag-fieldexpansion/mylar+z/root/root-75Mil-05/5.*.root");
 //   chain.Add("/mnt/gpfs3_amd/scratch/rgu245/Now/foil-mag-fieldexpansion/mylar+z/root/root-100Million-04/5.*.root");
    chain.SetBranchAddress("Te0", &Te0);
    chain.SetBranchAddress("dEeSilicon1", &dEeSilicon1);
    chain.SetBranchAddress("dEeSilicon2", &dEeSilicon2);
    chain.SetBranchAddress("timeHit1", &timeHit1);
    chain.SetBranchAddress("timeHit2", &timeHit2);
    chain.SetBranchAddress("pz0_e", &pZe);
    chain.SetBranchAddress("TotalNoHits", &TotalNoHits); 
    chain.SetBranchAddress("pInFoil1z", &pInFoil1z);
    chain.SetBranchAddress("pOutFoilz", &pOutFoilz);
    //std::cout<<" Number of Entries for sixFsixF at 5: "<<chain.GetEntries()<<endl;
   // const Long64_t maxEntries = 100000000;
    //const Long64_t numEntries = std::min(maxEntries, chain.GetEntries());
 //   std::cout<<"Entering for LOpp nEntries : "<<numEntries<<endl;
    for (Long64_t i = 0; i < chain.GetEntries(); i++) {
        chain.GetEntry(i);
        KE->Fill(Te0);

       if(dEeSilicon1 != 0 || dEeSilicon2 != 0){
            EdepType0->Fill(dEeSilicon1+dEeSilicon2);
        }
// deposits energy on detector 2 type 0 
        if(dEeSilicon1 == 0 && dEeSilicon2 != 0){
		    if(pZe > 0){
		        file2<<TotalNoHits<<" "<<Te0<<" "<<dEeSilicon1<<" "<<dEeSilicon2<<" "<<pZe<<endl;

		    }
            EdepType02->Fill(dEeSilicon2);
        }
// deposits energy on detector 1 type 0 
        if(dEeSilicon2 == 0 && dEeSilicon1 != 0){
           
		 EdepType01->Fill(dEeSilicon1);
            if(pZe < 0){
                file1<<TotalNoHits<<" "<<Te0<<" "<<dEeSilicon1<<" "<<dEeSilicon2<<" "<<pZe<<endl;
            }
        }
        if(dEeSilicon1 != 0 && dEeSilicon2 != 0){
		if(dEeSilicon1 > 30){	
                if(timeHit1< timeHit2){
                    EdepType11->Fill(dEeSilicon1 + dEeSilicon2);
                }
                if(timeHit1 > timeHit2){
                    EdepType12->Fill(dEeSilicon1 + dEeSilicon2);
                }}
                file<<TotalNoHits<<" "<<Te0<<" "<<dEeSilicon1<<" "<<dEeSilicon2<<" "<<timeHit1<<" "<<timeHit2<<" "<<pInFoil1z<<" "<<pOutFoilz<<" "<<pZe<<endl;
             
        }    
    }
     file.close();
   // file0.close();
   fout->Write();
   fout->Close();
}
