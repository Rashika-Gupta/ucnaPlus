#include "UCNBPrimaryGeneratorAction.hh"
#include "UCNBDetectorConstruction.hh"
#include "UCNBAnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4Proton.hh"
#include "G4VEmModel.hh"
 #include "globals.hh"

UCNBPrimaryGeneratorAction::UCNBPrimaryGeneratorAction(UCNBDetectorConstruction* myDC)
  :UCNBDetector(myDC)
{
}

UCNBPrimaryGeneratorAction::~UCNBPrimaryGeneratorAction()
{
  delete particleGun;
}

void UCNBPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
 // G4cout<<" Generating primaries    :PGA"<<G4endl;

  // Electron-Proton Event Generator

  G4double MN = 939.565379*1e6; // eV
  G4double MP = 938.272046*1e6; // eV
  G4double ME = 510.998928*1e3; // eV

  G4double MNk = 939565.379*keV; // keV
  G4double MPk = 938272.046*keV; // keV
  G4double MEk = 510.998928*keV; // keV

  G4double elambda = -1.2701;
  G4double alpha = 1./137.036;

  G4double POL = 1.0;
//--------------------------------------------------------
  G4double MN1 = MN/ME;
  G4double MP1 = MP/ME;
  G4double ME1 = 1.;
 
  //--------------------------------------------------------

  G4double A = -2.*elambda*(1.+elambda)/(1.+3.*elambda*elambda);
  G4double B = -2.*elambda*(1.-elambda)/(1.+3.*elambda*elambda);
  G4double ALIT = (1.-elambda*elambda)/(1.+3.*elambda*elambda);
//  G4double el1 = elambda;
  G4double BLIT = 0.0;
  //--------------------------------------------------------
  G4double gg = 1.08;
  G4double ff = 0.00;
  G4double EE, PE,cos_theta, PHIE, thetaE;

  while (gg>ff) {
    EE = ENERGY2() + 1.;
    PE = sqrt(EE*EE-1.);
 //    THETAE = M_PI*G4UniformRand();
    G4double cos_theta = std::cos(M_PI * G4UniformRand()); // Sample cos(theta) uniformly from 0 to 1, then take the cosine to ensure the range is from -1 to 1
    G4double thetaE = acos(cos_theta);
    PHIE = 2.*M_PI*G4UniformRand();
    ff = aprob2(A,elambda,thetaE,EE);
    gg = 1.08*G4UniformRand();
  }
  
  G4double EU = cos(PHIE)*sin(thetaE);
  G4double EV = sin(PHIE)*sin(thetaE);
  G4double EW = cos_theta;
// 
//  // Electron variables
  G4double Te0 = (EE-1.0) * ME/1000*keV;
  G4double px_hat_e = EU;
  G4double py_hat_e = EV;
  G4double pz_hat_e = EW;
  G4double thetaElectron = acos(pz_hat_e);
 
  // Sample vertex position
   G4double x_vertex, y_vertex, z_vertex;
  x_vertex = 0.0*m;//x_test*m;
  y_vertex = 0.0*m;//y_test*m;
  z_vertex = (-1.5 + G4UniformRand()*3.0)*m;
  
  
  // Generate electron and proton
  G4int n_particle = 1;
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  particleGun = new G4ParticleGun(n_particle);
  G4ParticleDefinition* particle1 = particleTable->FindParticle("e-");
  particleGun->SetParticleDefinition(particle1);
  particleGun->SetParticleEnergy(Te0);
  particleGun->SetParticleMomentumDirection(G4ThreeVector(px_hat_e,py_hat_e,pz_hat_e));

  particleGun->SetParticlePosition(G4ThreeVector(x_vertex,y_vertex,z_vertex));
  particleGun->GeneratePrimaryVertex(anEvent);

}

double UCNBPrimaryGeneratorAction::ENERGY2(){
  while( 1 ) {
    G4double E0 = 2.5295196;
    G4double b = 0.0;
    G4double f = 0.00;
    G4double y = 1.80;

    G4double E, FERMI;
    while (f<y) {
      E=(E0-1.0)*G4UniformRand();
      y=1.80*G4UniformRand();
      FERMI = 1.;
      f=FERMI*sqrt(E*E+2*E)*(E0-(E+1))*(E0-(E+1))*(E+1)*(1+b*1/(E+1));
    }
    return E;

  }
}

double UCNBPrimaryGeneratorAction::aprob2(double A, double lambda,
					  double theta, double E) {
  while( 1 ) {
    double beta = sqrt(E*E - 1.)/E;
    double temp = sin(theta)*(1.0+A*beta*cos(theta));
    return temp;
  }
}
