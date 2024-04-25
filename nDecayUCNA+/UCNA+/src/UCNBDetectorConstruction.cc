#include "UCNBDetectorConstruction.hh"
#include "UCNBPrimaryGeneratorAction.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4UserLimits.hh"
#include "G4SystemOfUnits.hh"
#include "G4FieldManager.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4TransportationManager.hh"
#include "G4PropagatorInField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4EqMagElectricField.hh"
#include "G4UniformMagField.hh"


#include "G4SystemOfUnits.hh"
#include "G4Polycone.hh"
UCNBDetectorConstruction::UCNBDetectorConstruction()
  : solidWorld(0), logicalWorld(0), physicalWorld(0), stepLimitMyl(0), stepLimitDead(0), stepLimitCar(0)
{
}

UCNBDetectorConstruction::~UCNBDetectorConstruction()
{
  delete stepLimitMyl;
  delete stepLimitDead;
  delete stepLimitCar;
}

G4VPhysicalVolume* UCNBDetectorConstruction::Construct()
{

  // Material definitions
  G4double a, z;
  G4double density, temperature, pressure;
  G4int nel, ncomponents, natoms;

  // Air at STP
  density = 1.293*mg/cm3;

 G4double torr = atmosphere/760.;
G4double nAtoms;
  G4Element* N = new G4Element("Nitrogen", "N", z=7., a=14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8., a=16.00*g/mole);
  G4Material* Air = new G4Material("Air", density, nel=2);
  Air->AddElement(N, 70*perCent);
  Air->AddElement(O, 30*perCent);
  // "Vacuum": 1e-6 Torr (low density air taken proportional to pressure)
  G4double densityAirSTP = density;
  pressure = 1e-99; // Torr
  density  = (pressure/760.) * densityAirSTP;
  //G4double density1= 1e-25*g/cm3;
  G4Material* Vacuum = new G4Material("Vacuum", density, nel=2);
  Vacuum->AddElement(N, 70*perCent);
  Vacuum->AddElement(O, 30*perCent);

  // Silicon: properties from wikipedia
  G4Material *Silicon = new G4Material("Silicon", z=14., a=28.086*g/mole, density=2.329*g/cm3);
  // Copper: properties from wikipedia
  G4Material *Copper = new G4Material("Copper", z=29., a=63.546*g/mole, density=8.94*g/cm3);
  // Stainless steel
  G4int nSS = 6;
  G4double fractionmass;
  density = 8.06*g/cm3;
  G4Material* SS = new G4Material("SS", density, nSS);
  G4Element* C = new G4Element("Carbon", "C", z=6., a=12.011*g/mole);
  G4Element* Si = new G4Element("Silicon", "Si", z=14., z=28.086*g/mole);
  G4Element* Cr = new G4Element("Chromium", "Cr", z=24., a=51.996*g/mole);
  G4Element* Mn = new G4Element("Manganese", "Mn", z=25., a=54.938*g/mole);
  G4Element* Fe = new G4Element("Iron", "Fe", z=26., a=55.845*g/mole);
  G4Element* Ni = new G4Element("Nickel", "Ni", z=28., a=58.693*g/mole);
  SS->AddElement(C,  fractionmass=0.001);
  SS->AddElement(Si, fractionmass=0.007);
  SS->AddElement(Cr, fractionmass=0.18);
  SS->AddElement(Mn, fractionmass=0.01);
  SS->AddElement(Fe, fractionmass=0.712);
  SS->AddElement(Ni, fractionmass=0.09);

  G4Material* Beryllium = new G4Material("Beryllium", z =4.,a = 9.01*g/mole, density = 1.848*g/cm3);
  G4Material* Aluminum =  new G4Material("Aluminum", z=13, a=26.9815*g/mole, density=2.70*g/cm3);
  G4Material* Tin = new G4Material("Tin", z=50, a=118.710*g/mole, density=7.365*g/cm3);
  G4Element* H = new G4Element("Hydrogen", "H", z=1., a=1.0008*g/mole);
  G4Material* Mylar = new G4Material("Mylar",1.370*g/cm3,3);
  Mylar->AddElement(C, fractionmass=0.62500);
  Mylar->AddElement(H, fractionmass=0.04167);
  Mylar->AddElement(O, fractionmass=0.33333);

  G4Material* sixFsixF = new G4Material("sixFsixF",1.480*g/cm3,5);
  G4Element* F = new G4Element("Fluorine", "F", z=9., a=18.9984*g/mole);
  sixFsixF->AddElement(H, fractionmass=0.27027);
  sixFsixF->AddElement(C, fractionmass=0.48648);
  sixFsixF->AddElement(N, fractionmass=0.02703);
  sixFsixF->AddElement(O, fractionmass=0.05405);
  sixFsixF->AddElement(F, fractionmass=0.16216);

  G4Material* PTFE = new G4Material("PTFE",2.20*g/cm3,2);
  PTFE->AddElement(C, fractionmass=0.24);
  PTFE->AddElement(F, fractionmass=0.76);

  G4Material* Alumina = new G4Material("Alumina",1.370*g/cm3,2);
  G4Element* Al = new G4Element("Aluminum","Al",z=13., a=26.9815*g/mole);
  Alumina->AddElement(Al,fractionmass=0.53);
  Alumina->AddElement(O, fractionmass=0.47);
 
 // G4Material* Be = new G4Material("Beryllium", z=4, a=9.0121*g/mole, density=1.85*g/cm3);
  G4Material* Germanium = new G4Material("Germanium", z=32., a=72.630*g/mole, density=5.323*g/cm3);
 
//***********************Material added from DetectorConstructionUtils.cc for the scintillator and photoguide editted on aug 16, 2020

// Scintillator, per Eljen EJ-204 datasheet
  G4Material*  Sci=new G4Material("Scintillator",1.032*g/cm3,2);
  Sci->AddElement(C,fractionmass=0.47609);
  Sci->AddElement(H,fractionmass=0.52390);
 //Wirechamber fill: N2 @95 torr
//  double P_MWPC = 100*torr;
//  double T_MWPC = 298*kelvin;

 // double P_N2 = P_MWPC - 5*torr;
 // G4Material* WCNitrogen = new G4Material("MWPC_N2",(28*mg)/(22.4*cm3)*P_N2/(760*torr)*(273.15*kelvin)/T_MWPC,1,kStateGas,T_MWPC,P_N2);
 // WCNitrogen->AddElement(G4Element::GetElement("N"),nAtoms=2);

//---------------------ooooOOOOOOOOOOooooooo---------------------------------------------

 
// Print all the material definitions
 G4cout << G4endl << "Material Definitions : " << G4endl << G4endl;
  //G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  // Generate & Add Material Properties Table ----------------------------------
  G4double stepSize=1.0;
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Definitions of Solids, Logical Volumes, Physical Volumes -----------------
  // World volume
  G4double world_x = 12.0*m;
  G4double world_y = 12.0*m;
  G4double world_z = 12.0*m;
  solidWorld    = new G4Box("World", world_x, world_y, world_z);
  logicalWorld  = new G4LogicalVolume(solidWorld, Vacuum, "World");
  physicalWorld = new G4PVPlacement(0, G4ThreeVector(), logicalWorld, "World", 0, false, 0);

  //G4double srcOffset_z = -13.0*cm;

  // Decay trap wall
 
  G4double innerRadiusDecayTrap = 0.065*m;
  G4double outerRadiusDecayTrap = 0.07*m;
  G4double halfLengthDecayTrap = (3.0/2.0)*m;
  G4double startAngleDecayTrap = 0.*deg;
  G4double spanningAngleDecayTrap = 360.*deg;

  G4double xDecayTrap = 0.0*m;
  G4double yDecayTrap = 0.0*m;
  G4double zDecayTrap = 0.0*m;

  DecayTrap = new G4Tubs("DecayTrap", innerRadiusDecayTrap, outerRadiusDecayTrap, halfLengthDecayTrap,
                         startAngleDecayTrap, spanningAngleDecayTrap);
  logicalDecayTrap = new G4LogicalVolume(DecayTrap, Copper, "DecayTrap");
  physicalDecayTrap = new G4PVPlacement(0, G4ThreeVector(xDecayTrap,yDecayTrap,zDecayTrap),
                                      logicalDecayTrap, "DecayTrap", logicalWorld, false, 0);
 

  // Decay trap interior
  G4double innerRadiusDecayTrapInt = 0.000*m;
  G4double outerRadiusDecayTrapInt = 0.065*m;
  G4double halfLengthDecayTrapInt = (3.0/2.0)*m;
  G4double startAngleDecayTrapInt = 0.*deg;
  G4double spanningAngleDecayTrapInt = 360.*deg;

  G4double xDecayTrapInt = 0.0*m;
  G4double yDecayTrapInt = 0.0*m;
  G4double zDecayTrapInt = 0.0*m;

  DecayTrapInt = new G4Tubs("DecayTrapInt", innerRadiusDecayTrapInt, outerRadiusDecayTrapInt, halfLengthDecayTrapInt,
                         startAngleDecayTrapInt, spanningAngleDecayTrapInt);
 logicalDecayTrapInt = new G4LogicalVolume(DecayTrapInt, Vacuum, "DecayTrapInt");
 physicalDecayTrapInt = new G4PVPlacement(0, G4ThreeVector(xDecayTrapInt,yDecayTrapInt,zDecayTrapInt),
                                     logicalDecayTrapInt, "DecayTrapInt", logicalWorld, false, 0);
 // Region between decay trap and Silicon detectors
  G4double innerRadiusDrift12 = 0.00*m;
  G4double outerRadiusDrift12 = 0.08*m;
 
  G4double startAngleDrift12 = 0.*deg;
  G4double spanningAngleDrift12 = 360.*deg;
  
  G4double halfLengthDrift12 = 0.35*m;
  G4double zDrift12 = 1.85*m; 
  G4double zDrift2 = -1.85*m; 
  
  G4double xDrift12 = 0.0*m;
  G4double yDrift12 = 0.0*m;
  Drift1 = new G4Tubs("Drift1", innerRadiusDrift12, outerRadiusDrift12, halfLengthDrift12,
                       startAngleDrift12, spanningAngleDrift12);
  logicalDrift1 = new G4LogicalVolume(Drift1, Vacuum, "Drift1");
  physicalDrift1 = new G4PVPlacement(0, G4ThreeVector(xDrift12,yDrift12,zDrift12),
                                    logicalDrift1, "Drift1", logicalWorld, false, 0);

//  G4double halfLengthDrift13 = 0.05*m;  for placing detector at different thickness
//  G4double zDrift13 = 1.55*m;
  Drift2 = new G4Tubs("Drift2", innerRadiusDrift12, outerRadiusDrift12, halfLengthDrift12,
		      startAngleDrift12, spanningAngleDrift12);
  logicalDrift2 = new G4LogicalVolume(Drift2, Vacuum, "Drift2");
  physicalDrift2 = new G4PVPlacement(0, G4ThreeVector(xDrift12,yDrift12,zDrift2),
				    logicalDrift2, "Drift2", logicalWorld, false, 0);
/*-------------------------- THIN FOIL ---------------------------------*/
/*need to place within the logical volume of the drift tu*/
  G4double dVary = 0.0*m;
  G4double dCoatingThick = 150*nm ;
  G4double dFoilThick = 130*nm ;
  G4double leftDrift = -0.35*m ;
  G4double rightDrift = +0.35*m ;
   // the leftmost point of the drift region wrt center of drift placed at 1.85*m
  G4double zEndCap = leftDrift + dVary;
  G4double zCoating = zEndCap + dCoatingThick/2;
  G4double zFoil  = dCoatingThick/2+zCoating+dFoilThick/2;
  //G4double zCoating =  dCoatingThick/2 + halfLengthDecayTrapInt ;
  //G4double zFoil =  dCoatingThick + halfLengthDecayTrapInt + dFoilThick/2  ;
 // G4double zFoil =  1.6*m  ;
  
  G4cout<<"zFoil    : "<<zFoil<<G4endl;
  BeTube1 = new G4Tubs("BeTube1", 0., outerRadiusDecayTrapInt, dCoatingThick/2., 0., 2*M_PI);
  logicalBeTube1 = new G4LogicalVolume(BeTube1, Beryllium, "BeTube1");
  physicalBeTube1 = new G4PVPlacement(0, G4ThreeVector(0,0,    zCoating), logicalBeTube1, "BeTube1", logicalDrift1, false, 0);
 
  foil1 = new G4Tubs("foil1", 0., outerRadiusDecayTrapInt, dFoilThick/2., 0., 2*M_PI);
  logicalfoil1 = new G4LogicalVolume(foil1, sixFsixF , "foil1");
  physicalfoil1 = new G4PVPlacement(0, G4ThreeVector(0,0,    zFoil), logicalfoil1, "foil1", logicalDrift1, false, 0);

  G4double zEndCap2 = rightDrift + dVary;
  G4double zCoating2 = zEndCap2 + dCoatingThick/2;
 // G4double zFoil2  = dCoatingThick/2+zCoating2+dFoilThick/2;
 // G4cout<<"ZFoil2  "<<zFoil2<<G4endl;
 // G4cout<<"zCoating 2 : "<<zCoating2<<G4endl;
  BeTube2 = new G4Tubs("BeTube2", 0., outerRadiusDecayTrapInt, dCoatingThick/2., 0., 2*M_PI);
  logicalBeTube2 = new G4LogicalVolume(BeTube2, Beryllium, "BeTube2");
  physicalBeTube2 = new G4PVPlacement(0, G4ThreeVector(0,0,    -1.*zCoating), logicalBeTube2, "BeTube2", logicalDrift2, false, 0);
 
  foil2 = new G4Tubs("foil2", 0., outerRadiusDecayTrapInt, dFoilThick/2., 0., 2*M_PI);
  logicalfoil2 = new G4LogicalVolume(foil2, sixFsixF , "foil2");
  physicalfoil2 = new G4PVPlacement(0, G4ThreeVector(0,0,    -1.*zFoil), logicalfoil2, "foil2", logicalDrift2, false, 0);

// Print out information about Drift2
G4cout << "Drift1 Position: " << physicalDrift1->GetTranslation() << G4endl;
G4cout << "Drift1 Dimensions: " << halfLengthDrift12 << " (half length)" << G4endl;
// Print out information about foil2
G4cout << "foil1 Position: " << physicalfoil1->GetTranslation() << G4endl;
G4cout << "foil2 Position: " << physicalfoil2->GetTranslation() << G4endl;


// removing the qindow source holder. 
/*Introducing source holder - windowTube */
// G4double windowRadius = 9.7*mm; //xuan
////G4double dSourceWindowThickness = 120*nm; //changing to 120nm 6F6F after 0.5 micron from a.young
//G4double dSourceWindowThickness = 0.5e-6*m; //changing to  0.5 so that can compare with same thickness material.
////G4double zWindow = halfLengthDecayTrapInt + dCoatingThick + dFoilThick +  dSourceWindowThickness/2 ;
//G4double zWindow = 0.0*m ;
//
//WindowTube = new G4Tubs("WindowTube", 0., windowRadius, dSourceWindowThickness/2,
//		      0., 360.*deg);
//logical_WindowTube = new G4LogicalVolume(WindowTube, Mylar, "WindowTube");
////physical_WindowTube = new G4PVPlacement(0, G4ThreeVector(0,0,0),
////				    logical_WindowTube, "WindowTube", logicalDecayTrapInt, false, 0);
//physical_WindowTube = new G4PVPlacement(0, G4ThreeVector(0,0,zWindow),
//				    logical_WindowTube, "WindowTube", logicalWorld, false, 0);
//
//G4cout<<"[DC]window tube constructed at : "<<zWindow/m<<G4endl;
 G4double si_dis=2.2;
  // Simple model for Silicon dead region
//  G4double tDead = 3e-6; // thickness of the dead layer was 80.0 e-9 for silicon was changed to 3 um from xuan's since the thickness of dead layer now is 3 um
//G4double tDead = 40.0e-9;  // Double thick dead layer
 //   G4double si_dis=1.5; //when there is no drift present the total length of the decay trap = 1.5 m  
/* G4double innerRadiusDead = 0.000*m;
  G4double outerRadiusDead = 0.08*m;
  G4double thicknessDead = (tDead)*m;
  G4double heightCylinderDead = thicknessDead/2.0;
  G4double startAngleDead = 0.*deg;
  G4double spanningAngleDead = 360.*deg;

 G4double xDead = 0.0*m;
 G4double yDead = 0.0*m;
 G4double zDead = ( si_dis+ tDead/2.0)*m;     


 Dead1 = new G4Tubs("Dead1", innerRadiusDead, outerRadiusDead, heightCylinderDead,
                     startAngleDead, spanningAngleDead);
  logicalDead1 = new G4LogicalVolume(Dead1, Sci, "Dead1"); //made Sci instead of Silicon - changed the material of the detector
 physicalDead1 = new G4PVPlacement(0, G4ThreeVector(xDead,yDead,zDead),
                                    logicalDead1, "Dead1", logicalWorld, false, 0);

 G4double zDead2 = ( si_dis + tDead/2.0)*m;   

  Dead2 = new G4Tubs("Dead2", innerRadiusDead, outerRadiusDead, heightCylinderDead,
                   startAngleDead, spanningAngleDead);
 logicalDead2 = new G4LogicalVolume(Dead2, Sci, "Dead2");
 physicalDead2 = new G4PVPlacement(0, G4ThreeVector(xDead,yDead,-1.*zDead2),
				    logicalDead2, "Dead2", logicalWorld, false, 0);
*/
  //G4double maxStepDL = stepSize*heightCylinderDead;
 //stepLimitDead = new G4UserLimits(maxStepDL);
 //logicalDead1->SetUserLimits(stepLimitDead);
  //logicalDead2->SetUserLimits(stepLimitDead);// *************what to do of it ??? *****

  // Simple model for Silicon active region
  // in the logical volume the silicon has been changed to scintillator material
  G4double tSi = 0.0035; //removing the deadlayer region = 0.0035 - tdead the thickness of the silicon detector was changed from 2 mm to 3.5 mm for scinitillator as in xuan's code but the variable is still the same.
 G4double innerRadiusSilicon = 0.0*m;
  G4double outerRadiusSilicon = 0.08*m;
//  G4double thicknessSilicon = (tSi - tDead)*m;
  G4double thicknessSilicon = tSi*m;
  G4double heightCylinder = thicknessSilicon/2.0;
  G4double startAngleSilicon = 0.*deg;
 G4double spanningAngleSilicon = 360.*deg;

 G4double xSilicon = 0.0*m;
 G4double ySilicon = 0.0*m;
// G4double zSilicon = (( ( si_dis+tDead) + ( si_dis+tSi) ) /2.0)*m;
 G4double zSilicon = (( ( si_dis) + ( si_dis+tSi) ) /2.0)*m;

 Silicon1 = new G4Tubs("Silicon1", innerRadiusSilicon, outerRadiusSilicon, heightCylinder,
                        startAngleSilicon, spanningAngleSilicon);
 logicalSilicon1 = new G4LogicalVolume(Silicon1, Sci, "Silicon1");
 physicalSilicon1 = new G4PVPlacement(0, G4ThreeVector(xSilicon,ySilicon,zSilicon),
                                     logicalSilicon1, "Silicon1", logicalWorld, false, 0);
// G4double zSilicon2 = (( ( si_dis+tDead) + ( si_dis+tSi) ) /2.0)*m;
 G4double zSilicon2 = (( ( si_dis) + ( si_dis+tSi) ) /2.0)*m;

 Silicon2 = new G4Tubs("Silicon2", innerRadiusSilicon, outerRadiusSilicon, heightCylinder,
			startAngleSilicon, spanningAngleSilicon);
 logicalSilicon2 = new G4LogicalVolume(Silicon2, Sci, "Silicon2");
 physicalSilicon2 = new G4PVPlacement(0, G4ThreeVector(xSilicon,ySilicon,-1.*zSilicon2),
                                   logicalSilicon2, "Silicon2", logicalWorld, false, 0);
G4cout<<"Detector constructed :     DC"<<G4endl;
/////////////////////////////////////////////////////////////////////////////////////////////
//************Detector properties from Xuan Scintillator Construction.hh ************///////

//G4double inch = 2.54*cm;
//G4double dScintRadius = 7.5*cm;
//G4double dBackingRadius = 10*cm;
//G4double dScintThick = 3.5*mm;
//G4double dDeadThick = 3.0*um;
//G4double dBackingThick = 1.*inch;
//G4double dLightGuideThick = 1.0*cm;

//G4double  dN2Volume_Z = dLightGuideThick + dBackingThick;
//G4double  dScintFace_PosZ = -dN2Volume_Z/2.;

//---- Create the shapes used in the scintillator object
//  // Overall container layer for the scintillator
  //  G4Tubs* N2VolTube = new G4Tubs("N2_vol_tube", 0., dBackingRadius, dN2Volume_Z/2., 0., 2*M_PI);
//scintOverall_shape = N2VolTube; 

// dead layer in scint
// G4Tubs* deadLayerTube = new G4Tubs("Dead_scint_tube", 0, dScintRadius, dDeadThick/2., 0., 2*M_PI);
// G4VisAttributes* visDScint= new G4VisAttributes(G4Colour(1.0,0.0,1.0,0.5));



//////////////////////

  // Electric and magnetic field definitions
  Field = new UCNBField();
  pEquation = new G4EqMagElectricField(Field);

  G4int nvar = 8;
  pStepper = new G4ClassicalRK4 (pEquation, nvar);

  pFieldMgr = G4TransportationManager::GetTransportationManager()->GetFieldManager();
  //pIntgrDriver = new G4MagInt_Driver(0.000001*mm,pStepper,pStepper->GetNumberOfVariables() );
  pIntgrDriver = new G4MagInt_Driver(1e-5*m,pStepper,pStepper->GetNumberOfVariables() );
  //pIntgrDriver = new G4MagInt_Driver(1.0*m,pStepper,pStepper->GetNumberOfVariables() );


  pChordFinder = new G4ChordFinder(pIntgrDriver);
  pFieldMgr->SetChordFinder( pChordFinder );
  pFieldMgr->GetChordFinder()->SetDeltaChord(1e-5*m);
  //pFieldMgr->GetChordFinder()->SetDeltaChord(1e-2*m);
  pFieldMgr->SetFieldChangesEnergy(true);
  pFieldMgr->SetDetectorField(Field);

  // Test: as recommended by Emil Frlez
  G4double myepsmin = 1.0e-5;
  pFieldMgr->SetMinimumEpsilonStep(myepsmin);
  pFieldMgr->SetMaximumEpsilonStep(myepsmin);
  pFieldMgr->SetDeltaOneStep(1.0e-4*mm);



  G4TransportationManager::GetTransportationManager()->GetPropagatorInField()->SetMaxLoopCount(INT_MAX);

  // Test: as recommended by Emil Frlez
  G4TransportationManager::GetTransportationManager()->GetPropagatorInField()->SetLargestAcceptableStep(0.1*mm); //

  // Visualization attributes

  logicalWorld->SetVisAttributes (G4VisAttributes::Invisible);
//  logicalDrift1->SetVisAttributes (G4VisAttributes::Invisible);
//  logicalDrift2->SetVisAttributes (G4VisAttributes::Invisible);
  logicalDecayTrapInt->SetVisAttributes (G4VisAttributes::Invisible);

  G4VisAttributes* simpleBoxVisAttGreen = new G4VisAttributes(G4Colour(0.0,1.0,0.0));
  simpleBoxVisAttGreen->SetVisibility(true);
  simpleBoxVisAttGreen->SetForceSolid(true);

 logicalSilicon1->SetVisAttributes(simpleBoxVisAttGreen);
 logicalSilicon2->SetVisAttributes(simpleBoxVisAttGreen);

  G4VisAttributes* simpleBoxVisAttRed = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
  simpleBoxVisAttRed->SetVisibility(true);
  simpleBoxVisAttRed->SetForceSolid(true);

// logicalDead1->SetVisAttributes(simpleBoxVisAttRed);
// logicalDead2->SetVisAttributes(simpleBoxVisAttRed);
logicalDecayTrapInt->SetVisAttributes(simpleBoxVisAttRed);
  // logicalDrift1->SetVisAttributes(simpleBoxVisAttRed);
  // logicalDrift2->SetVisAttributes(simpleBoxVisAttRed);

  G4VisAttributes* simpleBoxVisAttBlue = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  simpleBoxVisAttBlue->SetVisibility(true);
  simpleBoxVisAttBlue->SetForceSolid(true);

 logicalDecayTrap->SetVisAttributes(simpleBoxVisAttBlue);
 logicalDecayTrap->SetVisAttributes(G4VisAttributes::Invisible);
  // Return the geometry
  return physicalWorld;
}
