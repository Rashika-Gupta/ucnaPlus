#include <fstream>
#include <iostream>
#include "TBranch.h"
#include "TChain.h"
#include "TFile.h"

void mergefiles() {
    // Declare and initialize variables
    Double_t Te0, dEeSilicon1, dEeSilicon2,timeHit1,timeHit2, pZe;

    // Create a TChain
    TChain chain("Tout");
    chain.Add("/mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna/root-100-mil/1.*.root");
   
    chain.SetBranchAddress("dEeSilicon1", &dEeSilicon1);
    chain.SetBranchAddress("dEeSilicon2", &dEeSilicon2);
    chain.SetBranchAddress("timeHit1", &timeHit1);
    chain.SetBranchAddress("timeHit2", &timeHit2);
    chain.SetBranchAddress("pz0_e", &pZe);
    chain.SetBranchAddress("Te0", &Te0);
   
    // Create a TFile for output
    TFile *outputFile = new TFile("/mnt/gpfs3_amd/scratch/rgu245/Now/nDecayUCNA+/ucna.root", "RECREATE");
    TTree *outputTree = new TTree("outputTree", "Output Tree");

    // Create a TBranch for the accumulated energy
    Double_t e1, e2, ke, t1, t2, pz;
    TBranch *energy1 = outputTree->Branch("energyDet1", &e1, "energyDet1/D");
    TBranch *energy2 = outputTree->Branch("energyDet2", &e2, "energyDet2/D");
    TBranch *time1 = outputTree->Branch("timeHitDet1", &t1, "timeHitDet1");
    TBranch *time2 = outputTree->Branch("timeHitDet2", &t2, "timeHitDet2/D");
    TBranch *angle = outputTree->Branch("cosangle", &pz, "cosangle");
    TBranch *initKE= outputTree->Branch("KE", &ke, "KE");

       // Loop over events in the TChain
    for (Long64_t i = 0; i < chain.GetEntries(); i++) {
        chain.GetEntry(i);
        // Fill the TBranch with the appropriate variable (e.g., dEeSilicon1 + dEeSilicon2)
        
        e1 = dEeSilicon1 ;
        e2 =  dEeSilicon2; 
        t1 = timeHit1;
        t2 = timeHit2;
        ke = Te0;
        pz = pZe;
        energy1->Fill();
        energy2->Fill();
        time1->Fill();
        time2->Fill();
        angle->Fill();
        initKE->Fill();
}
      //  std::cout<<"energyDeposited " << energyDeposited << std::endl;

    // Write the TChain to the output file
    outputTree->Write();

    // Close the output file
    outputFile->Close();
  //  delete outputFile;
}
