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

void energyThrucnaPE1(){
    Char_t temp[200];
    std::string filename ; 

    ofstream file, file1, file2,file0;
    // file1.open("/mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/nDecay+_100Mill_01_oppDirection.txt");
    // file2.open("/mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/nDecay+_100Mill_02_oppDirection.txt");		
      file.open("/mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucnaP_Ethr1_Type1.txt");
     file0.open("/mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucnaP_Ethr1_Type0.txt");
     fileCheck.open("/mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucnaP_Ethr1_check.txt");
    
 sprintf(temp, "/mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucnap_Ethr1.root");


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

    TH1D *EdepType01Thr0 = new TH1D("EdepType01Thr0" ,"EdepType01Thr0", 8000., 0., 800.); //creating histogrms
    TH1D *EdepType11Thr0 = new TH1D("EdepType11Thr0" ,"EdepType11Thr0", 8000., 0., 800.);
    TH1D *EdepType02Thr0 = new TH1D("EdepType02Thr0" ,"EdepType02Thr0", 8000., 0., 800.); //creating histogrms
    TH1D *EdepType12Thr0 = new TH1D("EdepType12Thr0" ,"EdepType12Thr0", 8000., 0., 800.);
    TH1D *EdepType0Thr0  = new TH1D("EdepType0Thr0" ,"EdepType0Thr0", 8000., 0., 800.); //creating histogrms
    TH1D *EdepType1Thr0  = new TH1D("EdepType1Thr0" ,"EdepType1Thr0", 8000., 0., 800.);
//saving the evnets that wouold also be emiited in opposite direction. 
    Int_t TotalNoHits;
    Double_t counter = 0;
    TChain chain("Tout");
    chain.Add("/mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/UCNA+/100-million/6.*.root");
 //   chain.SetBranchAddress("Te0", &Te0);
    chain.SetBranchAddress("dEeSilicon1", &dEeSilicon1);
    chain.SetBranchAddress("dEeSilicon2", &dEeSilicon2);
    chain.SetBranchAddress("timeHit1", &timeHit1);
    chain.SetBranchAddress("timeHit2", &timeHit2);
    chain.SetBranchAddress("pz0_e", &pZe);
    chain.SetBranchAddress("TotalNoHits", &TotalNoHits); 
    chain.SetBranchAddress("pInFoil1z", &pInFoil1z);
    chain.SetBranchAddress("pOutFoilz", &pOutFoilz);
    Double_t E_thr = 10;
    std::cout <<"EThr : "<<E_thr;
    for (Long64_t i = 0; i <  10; i++) {
        chain.GetEntry(i);
        KE->Fill(Te0);
        Double_t E1 = dEeSilicon1;
        Double_t E2 = dEeSilicon2;
        
        std::cout<<" E1 : "<<dEeSilicon1 <<" E2 : "<<dEeSilicon2<<endl;
        if(E1 < E_thr){
            Double_t x  = 0;
            E1 = x;
        }
        
        if(E2 < E_thr){
            Double_t y  = 0;
            E2 = y;
        }
        if((E1 !=0 && E2 == 0) ||(E2 !=0 && E1 == 0)){
            EdepType0->Fill(E1+E2);
            file0<<Te0<<" "<<E1<<" "<<E2<<" "<<pZe<<endl;
        }
// deposits energy on detector 2 type 0 
        if(E1 == 0 && E2 != 0){
		    EdepType02->Fill(E2);

        }
// deposits energy on detector 1 type 0 
        if(E2 == 0 && E1 != 0){
           
		 EdepType01->Fill(E1);
            
        }
        if((dEeSilicon1 !=0 && dEeSilicon2 == 0) ||(dEeSilicon2 !=0 && dEeSilicon1 == 0)){
            EdepType0Thr0->Fill(dEeSilicon1+dEeSilicon2);
            
        }
        if((dEeSilicon1 != E1) || (dEeSilicon2 != E2)){
            fileCheck<<Te0<<" "<<dEeSilicon1<<" "<<dEeSilicon1<<" "<<E1<<" "<<E2<<" "<<pZe<<endl;
        }
// deposits energy on detector 2 type 0 
        if(dEeSilicon1 == 0 && dEeSilicon2 != 0){
		    EdepType02->Fill(dEeSilicon1 + dEeSilicon2);

        }
// deposits energy on detector 1 type 0 
        if(dEeSilicon2 == 0 && dEeSilicon1 != 0){
           
		 EdepType01->Fill(dEeSilicon1 + dEeSilicon2);
            
        }
        if(E1 != 0 && E2 != 0){
		   EdepType1->Fill(E1 + E2);
        
                if(timeHit1< timeHit2){
                    EdepType11->Fill(E1 + E2);
                }
                if(timeHit1 > timeHit2){
                    EdepType12->Fill(E1 + E2);
                }
                file<<Te0<<" "<<E1<<" "<<E2<<" "<<timeHit1<<" "<<timeHit2<<" "<<pZe<<endl;
             
        }
        if(dEeSilicon1 != 0 && dEeSilicon2 != 0){
		   EdepType1Thr0->Fill(dEeSilicon1 + dEeSilicon2);
        
                if(timeHit1< timeHit2){
                    EdepType11Thr0->Fill(dEeSilicon1+dEeSilicon2);
                }
                if(timeHit1 > timeHit2){
                    EdepType12Thr0->Fill(dEeSilicon1+dEeSilicon2);
                }
        }    
    }
     file.close();
     fileCheck.close();
    file0.close();
   fout->Write();
   fout->Close();
}
