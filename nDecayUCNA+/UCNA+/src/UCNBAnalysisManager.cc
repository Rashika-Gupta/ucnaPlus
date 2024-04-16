
#include "UCNBAnalysisManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "UCNBStackingAction.hh"
#include "Randomize.hh"
#include "G4Poisson.hh"
//-------------------------------------------------------------------------
UCNBAnalysisManager* UCNBAnalysisManager::fManager = 0;
//------------------------------------------------------------------------
UCNBAnalysisManager* UCNBAnalysisManager::getInstance()
{
  //G4cout<<"getting manager "<<G4endl;
  if(!fManager) {
    fManager = new UCNBAnalysisManager();
  }
  return fManager;
}
void UCNBAnalysisManager::dispose()
{
  G4cout<<"tryinh to delete from AM" <<G4endl;
  delete fManager;
  fManager = 0;
}
//.................................................................

UCNBAnalysisManager::UCNBAnalysisManager()
{
  verbose = 5;
  fileROOT = 0;
}

//................................................................

UCNBAnalysisManager::~UCNBAnalysisManager()
{
}
//////////////////////////////////////////////////////////////////////////////////////////

void UCNBAnalysisManager::bookROOT()
{
  G4cout << "AM == *booking root**" << ROOTfilename << G4endl;
  fileROOT = new TFile( ROOTfilename, "RECREATE", "UCNB Simulation");
  Tout = new TTree("Tout","Tout");
  Tout->Branch("Te0",         &Te0,         "Te0/D");
  Tout->Branch("Tp0",         &Tp0,         "Tp0/D"); //Initial Energies at an event vertex.
  Tout->Branch("Tn0",         &Tn0,         "Tn0/D");
  Tout->Branch("Tv0",         &Tv0,         "Tv0/D");
  //---------------------------------------------------------
  //---------------------------------------------------------
  Tout->Branch("x0",          &x0,          "x0/D");
  Tout->Branch("y0",          &y0,          "y0/D"); //Initial Positions
  Tout->Branch("z0",          &z0,          "z0/D");
  //-----------------------------------------------------------
  Tout->Branch("px0_e",       &pXe,          "px0_e/D");
  Tout->Branch("py0_e",       &pYe,          "py0_e/D");//Initial Momentum direction of Electron
  Tout->Branch("pz0_e",       &pZe,          "pz0_e/D");
  //----------------------------------------------------------
  Tout->Branch("px0_p",       &pX,          "px0_p/D");
  Tout->Branch("py0_p",       &pY,          "py0_p/D"); //Initial Momentum direction of Proton
  Tout->Branch("pz0_p",       &pZ,          "pz0_p/D");
  //-----------------------------------------------------------
  Tout->Branch("thetae0",     &thetae0,     "thetae0/D"); //Electron initial angle relative to z-direction.
  Tout->Branch("thetap0",     &thetap0,     "thetap0/D"); //Proton initial angle relative to z-direction.
  Tout->Branch("px_in_det1", &px_in_det1, "px_in_det1/D");
  Tout->Branch("px_in_det2", &px_in_det2, "px_in_det2/D");

  Tout->Branch("py_in_det1", &py_in_det1, "py_in_det1/D");
  Tout->Branch("py_in_det2", &py_in_det2, "py_in_det2/D");

  Tout->Branch("pz_in_det1", &pz_in_det1, "pz_in_det1/D");
  Tout->Branch("pz_in_det2", &pz_in_det2, "pz_in_det2/D");

  Tout->Branch("px_out_det1", &px_out_det1, "px_out_det1/D");
  Tout->Branch("px_out_det2", &px_out_det2, "px_out_det2/D");

  Tout->Branch("py_out_det1", &py_out_det1, "py_out_det1/D");
  Tout->Branch("py_out_det2", &py_out_det2, "py_out_det2/D");
  Tout->Branch("pz_out_det1", &pz_out_det1, "pz_out_det1/D");
  Tout->Branch("pz_out_det2", &pz_out_det2, "pz_out_det2/D");

  Tout->Branch("pInFoil1x", &pInFoil1x, "pInFoil1x/D");
  Tout->Branch("pInFoil1y", &pInFoil1y, "pInFoil1y/D");
  Tout->Branch("pInFoil1z", &pInFoil1z, "pInFoil1z/D");

  Tout->Branch("pOutFoilx", &pOutFoilx, "pOutFoilx/D");
  Tout->Branch("pOutFoily", &pOutFoily, "pOutFoily/D");
  Tout->Branch("pOutFoilz", &pOutFoilz, "pOutFoilz/D");

  Tout->Branch("xe1",    &xeSilicon1,    "xe1/D");
  Tout->Branch("ye1",    &yeSilicon1,    "ye1/D");
  Tout->Branch("ze1",    &zeSilicon1,    "ze1/D");

  Tout->Branch("xe2",    &xeSilicon2,    "xe2/D");
	Tout->Branch("ye2",    &yeSilicon2,    "ye2/D");
	Tout->Branch("ze2",    &zeSilicon2,    "ze2/D");

  Tout->Branch("xe1D",    &xeDead1,    "xe1D/D");
	Tout->Branch("ye1D",    &yeDead1,    "ye1D/D");
	Tout->Branch("ze1D",    &zeDead1,    "ze1D/D");

  Tout->Branch("xe2D",    &xeDead2,    "xe2D/D");
	Tout->Branch("ye2D",    &yeDead2,    "ye2D/D");
	Tout->Branch("ze2D",    &zeDead2,    "ze2D/D");
  
//-----------------------------------------------------------------------------
  Tout->Branch("dEeOther",    &dEeOther,    "dEeOther/D");
  Tout->Branch("dEeSilicon1", &dEeSilicon1, "dEeSilicon1/D");
  Tout->Branch("dEeSilicon2", &dEeSilicon2, "dEeSilicon2/D");
//  Tout->Branch("dEeDead1",    &dEeDead1,    "dEeDead1/D");
//  Tout->Branch("dEeDead2",    &dEeDead2,    "dEeDead2/D");
  Tout->Branch("dEeFoil1",    &dEeFoil1,    "dEeFoil1/D");
  Tout->Branch("dEeFoil2",    &dEeFoil2,    "dEeFoil2/D");
  Tout->Branch("dEeBeTube1",    &dEeBeTube1,    "dEeBeTube1/D");
  Tout->Branch("dEeBeTube2",    &dEeBeTube2,    "dEeBeTube2/D");

//----------------------------------------------------------------------------------
  //Tout->Branch("dEpOther",    &dEpOther,    "dEpOther/D");
  //Tout->Branch("dEpSilicon1", &dEpSilicon1, "dEpSilicon1/D");
  //Tout->Branch("dEpSilicon2", &dEpSilicon2, "dEpSilicon2/D");
  //Tout->Branch("dEpDead1",    &dEpDead1,    "dEpDead1/D");
  //Tout->Branch("dEpDead2",    &dEpDead2,    "dEpDead2/D");
//---------------------------------------------------------------------------------
        Tout->Branch("dEeSilicon1Gauss", &dEeSilicon1Gauss, "dEeSilicon1Gauss/D");
        Tout->Branch("dEeSilicon2Gauss", &dEeSilicon2Gauss, "dEeSilicon2Gauss/D");
//---------------------------------------------------------------------------------
  Tout->Branch("eTOF",        &eTOF,        "eTOF/D");
  Tout->Branch("pTOF",        &pTOF,        "pTOF/D");

  Tout->Branch("etof",        &etof,        "etof/D");
  Tout->Branch("ptof",        &ptof,        "ptof/D");

  Tout->Branch("Etof",        &Etof,        "Etof/D");
  Tout->Branch("Ptof",        &Ptof,        "Ptof/D");
//---------------------------------------------------------------
  Tout->Branch("timeHit1",    &timeHit1,    "timeHit1/D");
  Tout->Branch("timeHit2",    &timeHit2,    "timeHit2/D");

//  Tout->Branch("timeDead1",   &timeDead1,   "timeDead1/D");
//  Tout->Branch("timeDead2",   &timeDead2,   "timeDead2/D");
  
  Tout->Branch("eEdepTime1",  &eEdepTime1,  "eEdepTime1[500]/D");
  Tout->Branch("eEdepTime2",  &eEdepTime2,  "eEdepTime2[500]/D");

//  Tout->Branch("eEdepDeadTime1",  &eEdepDeadTime1,  "eEdepDeadTime1[500]/D");
//  Tout->Branch("eEdepDeadTime2",  &eEdepDeadTime2,  "eEdepDeadTime2[500]/D");
  
//-----------------------------------------------------------------------------
  Tout->Branch("timeHit1P",    &timeHit1P,    "timeHit1P/D");
  Tout->Branch("timeHit2P",    &timeHit2P,    "timeHit2P/D");

//  Tout->Branch("timeDead1P",   &timeDead1P,   "timeDead1P/D");
//  Tout->Branch("timeDead2P",   &timeDead2P,   "timeDead2P/D");
  
  Tout->Branch("pEdepTime1P",  &pEdepTime1,  "pEdepTime1[500]/D");
  Tout->Branch("pEdepTime2P",  &pEdepTime2,  "pEdepTime2[500]/D");

//  Tout->Branch("pEdepDeadTime1",  &pEdepDeadTime1,  "pEdepDeadTime1[500]/D");
//  Tout->Branch("pEdepDeadTime2",  &pEdepDeadTime2,  "pEdepDeadTime2[500]/D");
  
//-------------------------------------------------------------------------------
   Tout->Branch("dESi1Tr", &dESi1Tr, "dESi1Tr[200]/D");
  Tout->Branch("dESi2Tr", &dESi2Tr, "dESi2Tr[200]/D");
  Tout->Branch("BtimeHitSi1",  &BtimeHitSi1,  "BtimeHitSi1[200]/D"); // Deprecated "multi-hit/backscatter
  Tout->Branch("BtimeHitSi2",  &BtimeHitSi2,  "BtimeHitSi2[200]/D");  //
  Tout->Branch("BdESi1Tr", &BdESi1Tr, "BdESi1Tr[200]/D"); //Deprecated!
  Tout->Branch("BdESi2Tr", &BdESi2Tr, "BdESi2Tr[200]/D"); //

  Tout->Branch("dESi1Hit", dESi1Hit, "dESi1Hit[5]/D"); //
  Tout->Branch("dESi2Hit", dESi2Hit, "dESi2Hit[5]/D"); 
  Tout ->Branch("TotalNoHits", &TotalNoHits, "TotalNoHits/I");
 // Tout->Branch("dESi1Hit", &dESi1Hit, "dESi1Hit[5]/D"); //
 // Tout->Branch("dESi2Hit", &dESi2Hit, "dESi2Hit[5]/D"); // formerly kept track of timing from 2 hits/detector
  Tout->Branch("dESi1HitTime", &dESi1HitTime, "dESi1HitTime[5]/D");//
  Tout->Branch("dESi2HitTime", &dESi2HitTime, "dESi2HitTime[5]/D");// Still deprecated
//-------------------------------------------------------------------------------

  Tout->Branch("Xpho",   &bremsX,   "Xpho/D");
  Tout->Branch("Ypho",   &bremsY,   "Ypho/D");
  Tout->Branch("Zpho",   &bremsZ,   "Zpho/D");
  Tout->Branch("pXpho0", &pXpho0,   "pXpho0/D");
  Tout->Branch("pYpho0", &pYpho0,   "pYpho0/D");
  Tout->Branch("pZpho0", &pZpho0,   "pZpho0/D");
  
     Tout->Branch("iKill",       &iKill,       "iKill/I");
     Tout->Branch("numGamma",  &numGamma,  "numGamma/I");
     Tout->Branch("numOther",  &numOther,  "numOther/I");//commented
     Tout->Branch("PID",  &PID,  "PID/I");


 // Tout->Branch("dESi1Hit", &dESi1Hit,   "dESi1Hit[5]/D");
 // Tout->Branch("dESi2Hit", &dESi2Hit,   "dESi2Hit[5]/D");
 // Tout->Branch("dESi1HitTime", &dESi1HitTime,   "dESi1HitTime[5]/D");
  //Tout->Branch("dESi2HitTime", &dESi2HitTime,   "dESi2HitTime[5]/D");
  
//  Tout->Branch("dEDead1Hit", &dEDead1Hit,   "dEDead1Hit[5]/D");
//  Tout->Branch("dEDead2Hit", &dEDead2Hit,   "dEDead2Hit[5]/D");
//  Tout->Branch("dEDead1HitTime", &dEDead1HitTime,   "dEDead1HitTime[5]/D");
//  Tout->Branch("dEDead2HitTime", &dEDead2HitTime,   "dEDead2HitTime[5]/D");
//  
//  Tout->Branch("dEDead1HitP", &dEDead1HitP,   "dEDead1HitP[5]/D");
//  Tout->Branch("dEDead2HitP", &dEDead2HitP,   "dEDead2HitP[5]/D");
//  Tout->Branch("dEDead1HitTimeP", &dEDead1HitTimeP,   "dEDead1HitTimeP[5]/D");
//  Tout->Branch("dEDead2HitTimeP", &dEDead2HitTimeP,   "dEDead2HitTimeP[5]/D");
  
  Tout->Branch("dESi1HitP", &dESi1HitP,   "dESi1HitP[5]/D");
  Tout->Branch("dESi2HitP", &dESi2HitP,   "dESi2HitP[5]/D");
  Tout->Branch("dESi1HitTimeP", &dESi1HitTimeP,   "dESi1HitTimeP[5]/D");
  Tout->Branch("dESi2HitTimeP", &dESi2HitTimeP,   "dESi2HitTimeP[5]/D");

}
///////////////////////////////////////////////////////////////////////////

void UCNBAnalysisManager::BeginOfRun()
{
  evNo=-1;
  bookROOT();
  G4cout<<"AM Beginning of run ." << G4endl;
}
//  G4cout << ">>>>>>>>>> UCNBAnalysisManager: ROOT Ntuple is booked ..." << G4endl;

//---------------------------------------------------------------------------
void UCNBAnalysisManager::EndOfRun()
{
  G4cout << "AM- END OF RUN " << G4endl;
  saveROOT();
}
//----------------------------------------------------------------------------

void UCNBAnalysisManager::BeginOfEvent()
{
  evNo++;
  if(evNo%1==0) G4cout << "EVENT NUMBER ------- " << evNo << G4endl;
  G4int Det1Hits = 0.;
  G4int Det2Hits = 0.;
  dETotal = 0.;
  dummydummy=0;
  dummydummy2=0;
  G4cout << evNo << "Event ----------------------------------------" << G4endl;
  HitNo=0.;
 G4double HitNo1=0.;
  HitNo2=0.;
  
  HitNo1d=0.;
  HitNo2d=0.;
  
  HitNo11=0.;
  HitNo22=0.;
  
  HitNo1s=0.;
  HitNo2s=0.; 
  TotalNoHits = 0.;
  //px_in_det1 = 0.;
 // G4cout<<" AM    TotalNoHits "<< TotalNoHits <<G4endl;
  for(Int_t pp=0; pp<5; pp++){
  //  G4cout<<" pp value is :"<<pp<<G4endl;
  dESi1Hit[pp]=0.0;
  dESi2Hit[pp]=0.0;
  dESi1HitTime[pp]=0.0;
  //G4cout<<"AM     dESi1HitTime [pp]"<<dESi1HitTime[pp]<<G4endl;
  dESi2HitTime[pp]=0.0;
 // dEDead1Hit[pp]=0.0;
 // dEDead2Hit[pp]=0.0;
 // dEDead1HitTime[pp]=0.0;
 // dEDead2HitTime[pp]=0.0;
 // 
 // dEDead1HitP[pp]=0.0;
 // dEDead2HitP[pp]=0.0;
 // dEDead1HitTimeP[pp]=0.0;
 // dEDead2HitTimeP[pp]=0.0;  
  dESi1HitP[pp]=0.0;
  dESi2HitP[pp]=0.0;
  dESi1HitTimeP[pp]=0.0;
  dESi2HitTimeP[pp]=0.0;
  }
   //G4cout<<"AM     dESi1HitTime [pp]"<<dESi1HitTime[pp]<<G4endl;
  for(Int_t qik=0; qik<15; qik++)
  dEeOther = 0.;
  dEeSilicon1 = 0.;
  dEeSilicon2 = 0.;
//  dEeDead1 = 0.;
//  dEeDead2 = 0.;
  dEeFoil1 = 0.;
  dEeFoil2 = 0.;

  timeFirst1=0;
  timeLast1=0;
  timeFirst2=0;
  timeLast2=0;
  PID=0;
  dEeSilicon1Gauss = 0.;
  dEeSilicon2Gauss = 0.;

  dEpOther = 0.;
  dEpSilicon1 = 0.;
  dEpSilicon2 = 0.;
  dEpDead1 = 0.;
  dEpDead2 = 0.;

  globalTimeHit1 = -5.0*1e-9;
 // G4cout <<"globalTimeHit1 ln 283 AM : "<<globalTimeHit1<<G4endl;
  globalTimeHit2 = -5.0*1e-9;

//  globalTimeDead1 = -5.0*1e-9;
//  globalTimeDead2 = -5.0*1e-9;

  globalTimeHit1P = -5.0*1e-9;
  globalTimeHit2P = -5.0*1e-9;

//  globalTimeDead1P = -5.0*1e-9;
//  globalTimeDead2P = -5.0*1e-9;


  pInFoil1x = 0;
  pInFoil1y = 0;
  pInFoil1z = 0;
  pOutFoilx = 0;
  pOutFoily = 0;
  pOutFoilz = 0;
  eTOF = 0.;
  pTOF = 0.;
  etof=0.;
  ptof=0.;

  siHit=0;

//---------------------------------------------------------------------
  iKill = -1;
  bremsEdep = 0;
  bremsID = 0;
  bremsIter = 0;
  eSecIter = 0;
  gammaIter = 0;
  otherIter = 0;
  postStepCount = 0;
  numGamma=0;
  numOther=0;

  for (G4int kk=0; kk<200; kk++) {
    EdepTimeBin1[kk] = 0.;
    EdepTimeBin2[kk] = 0.;
  //  EdepDeadTimeBin1[kk] = 0.;
  //  EdepDeadTimeBin2[kk] = 0.;

    EdepTimeBin1P[kk] = 0.;
    EdepTimeBin2P[kk] = 0.;
  //  EdepDeadTimeBin1P[kk] = 0.;
  //  EdepDeadTimeBin2P[kk] = 0.;

    timeHitSi1[kk]=0;
    timeHitSi2[kk]=0;

    dESi1Tr[kk]=0;
    dESi2Tr[kk]=0;

    BtimeHitSi1[kk]=0;
    BtimeHitSi2[kk]=0;

    BdESi1Tr[kk]=0;
    BdESi2Tr[kk]=0;
  }
}
//-----------------------------------------------------------------------------------------------------------
void UCNBAnalysisManager::EndOfEvent()
{
//  G4cout<<"Edead 1 : "<<dEeDead1<<G4endl;
  G4cout <<"AM == end of event "<< G4endl;
  if (globalTimeHit1 > 0. && globalTimeHit2 > 0.) {
    etof = globalTimeHit1 - globalTimeHit2;
  }
  else {
    etof = 0.;
  }

 if (globalTimeDead1P > 0. && globalTimeDead2P > 0.) {
    ptof = globalTimeDead1P - globalTimeDead2P;
  }
  else {
    ptof = 0.;
  }

 /////////////////New proton Time of flight detected by the detectors/////////////////
  if (globalTimeHit1P > 0. && globalTimeHit2P > 0.) {
    Ptof = globalTimeHit1P - globalTimeHit2P;
  }
  else {
    Ptof = 0.;
  }
// if (globalTimeDead1 > 0. && globalTimeDead2 > 0.) {
//    Etof = globalTimeDead1 - globalTimeDead2;
//  }
//  else {
//    Etof = 0.;
//  }

TotalNoHits=TotalNoHits;

//------------------Electron time---------------------------------------------------------------
  timeHit1  = globalTimeHit1;
 // G4cout<<"Time hit det 1 AM ln 374   :"<<timeHit1<<G4endl;
  timeHit2  = globalTimeHit2;
  timeDead1 = globalTimeDead1;
  timeDead2 = globalTimeDead2;

  //---------------Proton time-----------------------------------------------------------------------
  timeHit1P  = globalTimeHit1P;
  timeHit2P  = globalTimeHit2P;
  timeDead1P = globalTimeDead1P;
  timeDead2P = globalTimeDead2P;

  for (int kk=0; kk<500; kk++) {
    eEdepTime1[kk]     = EdepTimeBin1[kk];
    eEdepTime2[kk]     = EdepTimeBin2[kk];
  //  eEdepDeadTime1[kk] = EdepDeadTimeBin1[kk];
  //  eEdepDeadTime2[kk] = EdepDeadTimeBin2[kk];

    pEdepTime1[kk]     = EdepTimeBin1P[kk];
    pEdepTime2[kk]     = EdepTimeBin2P[kk];
    pEdepDeadTime1[kk] = EdepDeadTimeBin1P[kk];
    pEdepDeadTime2[kk] = EdepDeadTimeBin2P[kk];
  }


  if (dEeSilicon1 > 0.) dEeSilicon1Gauss = G4RandGauss::shoot(dEeSilicon1, 0.85);
  if (dEeSilicon2 > 0.) dEeSilicon2Gauss = G4RandGauss::shoot(dEeSilicon2, 0.85);

  if (dEeSilicon1Gauss < 0.) dEeSilicon1Gauss = 0.;
  if (dEeSilicon2Gauss < 0.) dEeSilicon2Gauss = 0.;
  Tout->Fill();
}  


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void UCNBAnalysisManager::p1incident(G4double px, G4double py, G4double pz)
{ px_in_det1 = px;
  py_in_det1 = py; 
  pz_in_det1 = pz;
}
void UCNBAnalysisManager::p2incident(G4double px, G4double py, G4double pz)
{ px_in_det2 = px;
  py_in_det2 = py; 
  pz_in_det2 = pz;
}
void UCNBAnalysisManager::p1out(G4double px, G4double py, G4double pz)
{ px_out_det1 = px;
  py_out_det1 = py;
  pz_out_det1 = pz;
}
void UCNBAnalysisManager::pInFoil(G4double px, G4double py, G4double pz)
{ pInFoil1x = px;
  pInFoil1y = py; 
  pInFoil1z = pz;
}
void UCNBAnalysisManager::pOutFoil(G4double px, G4double py, G4double pz)
{ pOutFoilx = px;
  pOutFoily = py; 
  pOutFoilz = pz;
}
void UCNBAnalysisManager::p2out(G4double px, G4double py, G4double pz)
{ px_out_det2 = px;
  py_out_det2 = py;
  pz_out_det2 = pz;
}
void UCNBAnalysisManager::AddUpTotalEnergyDeposit(G4double x)
{
//  G4cout<<"Add up total energy deposit "<<G4endl;
  dETotal = dETotal + x;
}
void UCNBAnalysisManager::AddUpElectronOtherEnergyDeposition(G4double x)
{
  dEeOther = dEeOther + x;
}
void UCNBAnalysisManager::AddUpProtonOtherEnergyDeposition(G4double x)
{
  dEpOther = dEpOther + x;
}

void UCNBAnalysisManager::AddUpElectronSilicon1EnergyDeposition(G4double x)
{
  dEeSilicon1 = dEeSilicon1 + x;
}
void UCNBAnalysisManager::AddUpElectronSilicon2EnergyDeposition(G4double x)
{
  dEeSilicon2 = dEeSilicon2 + x;
}
void UCNBAnalysisManager::AddUpElectronFoil1EnergyDeposition(G4double x)
{
  dEeFoil1 = dEeFoil1 + x;
}
void UCNBAnalysisManager::AddUpElectronFoil2EnergyDeposition(G4double x)
{
  dEeFoil2 = dEeFoil2 + x;
}
void UCNBAnalysisManager::AddUpElectronBeTube1EnergyDeposition(G4double x)
{
  dEeBeTube1 = dEeBeTube1 + x;
}
void UCNBAnalysisManager::AddUpElectronBeTube2EnergyDeposition(G4double x)
{
  dEeBeTube2 = dEeBeTube2 + x;
}
void UCNBAnalysisManager::AddUpProtonSilicon1EnergyDeposition(G4double x)
{
  dEpSilicon1 = dEpSilicon1 + x;  
}
void UCNBAnalysisManager::AddUpProtonSilicon2EnergyDeposition(G4double x)
{
  dEpSilicon2 = dEpSilicon2 + x;
}
//void UCNBAnalysisManager::AddUpElectronDeadLayer1EnergyDeposition(G4double x)
//{
//  dEeDead1 = dEeDead1 + x;
//}
//void UCNBAnalysisManager::AddUpElectronDeadLayer2EnergyDeposition(G4double x)
//{
//  dEeDead2 = dEeDead2 + x;
//}
void UCNBAnalysisManager::AddUpProtonDeadLayer1EnergyDeposition(G4double x)
{
  dEpDead1 = dEpDead1 + x;
}
void UCNBAnalysisManager::AddUpProtonDeadLayer2EnergyDeposition(G4double x)
{
  dEpDead2 = dEpDead2 + x;
}
void UCNBAnalysisManager::recordBremPos(G4double x, G4double y,G4double z, G4double px, G4double py, G4double pz)
{
  bremsX = x;
  bremsY = y;
  bremsZ = z;
  pXpho0 = px;
  pYpho0 = py;
  pZpho0 = pz;
}
void UCNBAnalysisManager::storePostPhoeE(G4double Ee)
{
  postStepEe[postStepCount]=Ee;
  postStepCount++;
}
void UCNBAnalysisManager::AddUpProtonDriftTime(G4double x)
{
  pTOF = pTOF + x;
}

void UCNBAnalysisManager::AddUpBremsEnergyDeposit(G4double x)
{
  bremsEdep = bremsEdep + x;
}

void UCNBAnalysisManager::AddUpElectronDriftTime(G4double x)
{
 eTOF = eTOF + x;
}
void UCNBAnalysisManager::saveEventVertex(G4double xe,G4double ye,G4double ze, G4double Te, G4double Tp,
  G4double pxi, G4double pyi, G4double pzi, G4double the, G4double thp, G4double Tn, G4double Tv)
{
  x0 = xe; // added on dec 9, 22
  y0 = ye; // added on dec 9, 22
  z0 = ze; // added on dec 9, 22
  Te0 = Te;
  Tp0 = Tp;
  pXe = pxi;
  pYe = pyi;
  pZe = pzi;
  thetae0 = the;
  thetap0 = thp;
  Tn0 = Tn;
  Tv0 = Tv;
}

void UCNBAnalysisManager::saveSourceVertex(G4double Te, G4double Tp,
  G4double xvtx, G4double yvtx, G4double zvtx, G4double pXi, G4double pYi, G4double pZi, G4int PIDi2)
{
  Te0 = Te;
  Tpho0 = Tp;
  x0 = xvtx;
  y0 = yvtx;
  z0 = zvtx;
  pX = pXi;
  pY = pYi;
  pZ = pZi;
  PID= PIDi2;
}
//---------------------Electron Position-------------------------------------------------------
void UCNBAnalysisManager::recordSilicon1ePosition(G4double x, G4double y, G4double z)
{
  xeSilicon1= x;
  yeSilicon1 = y;
  zeSilicon1 = z;  
}

void UCNBAnalysisManager::recordSilicon2ePosition(G4double x, G4double y,G4double z)
{
  xeSilicon2 = x;
  yeSilicon2 = y;
  zeSilicon2 = z;
}
void UCNBAnalysisManager::recordDead1ePosition(G4double x, G4double y, G4double z)
{
  xeDead1 = x;
  yeDead1 = y;
  zeDead1 = z;  
}
void UCNBAnalysisManager::recordDead2ePosition(G4double x, G4double y,G4double z)
{
  xeDead2 = x;
  yeDead2 = y;
  zeDead2 = z;
}

//-------------------------Proton Position--------------------------------------------------
void UCNBAnalysisManager::recordSilicon1pPosition(G4double x, G4double y, G4double z)
{
  xpSilicon1 = x;
  ypSilicon1 = y;
  zpSilicon1 = z;
}
void UCNBAnalysisManager::recordSilicon2pPosition(G4double x, G4double y, G4double z)
{
  xpSilicon2 = x;
  ypSilicon2 = y;
  zpSilicon2 = z;
}
void UCNBAnalysisManager::recordDead1pPosition(G4double x, G4double y,G4double z)
{
  xpDead1 = x;
  ypDead1 = y;
  zpDead1 = z;
}
void UCNBAnalysisManager::recordDead2pPosition(G4double x, G4double y, G4double z)
{
  xpDead2 = x;
  ypDead2 = y;
  zpDead2 = z;
}
//----------------------------------------------------------------------------------------
void UCNBAnalysisManager::killEventFlag(G4int i)
{
  iKill = i;
}
void UCNBAnalysisManager::saveROOT()
{
  fileROOT->Write();
  fileROOT->Close();
}

