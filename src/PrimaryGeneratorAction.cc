#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "G4EmCalculator.hh"
#include "G4NistManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4GeneralParticleSource.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4UImanager.hh"
#include "G4Proton.hh"
#include "G4Alpha.hh"
#include "G4Deuteron.hh"
#include "globals.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include "G4Proton.hh"
#include "G4IonTable.hh"
#include  <fstream>
#include  <sstream>
#include <math.h>
#include "G4NavigationHistory.hh"
#include "G4TouchableHistory.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TransportationManager.hh"
#include "TFile.h"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
using namespace std;

PrimaryGeneratorAction* PrimaryGeneratorAction::theGenerator = NULL;

PrimaryGeneratorAction::PrimaryGeneratorAction(G4double energy, G4int ANumber, G4int Nprot):ENER(energy), A(ANumber), Nprotons(Nprot)
{ 

  theGenerator = this;
  theDetector = DetectorConstruction::GetInstance();

  particleSource  = new G4GeneralParticleSource();

  //Generic beam (no phase space)
  if(A==1) particle = G4Proton::Proton();//G4IonTable::GetIonTable()->GetIon(1,1,0); // proton
  else if(A==2) particle = G4Deuteron::Deuteron();
  else if(A==4) particle = G4Alpha::Alpha();
  else particle = G4IonTable::GetIonTable()->GetIon(int(A/2),A,0); // rest
  particleSource->SetParticleDefinition(particle);
  

  // Mono-energetic
  //eneDist->SetEnergyDisType("Mono");
  //eneDist->SetMonoEnergy(ENER*A*MeV);
  
  ////////////////////////////////////////////////////////////////////////////////////
  //Source 1
  //Position

  //Angular Specturm
  PencilBeamStdAng = 25;
  angDist = particleSource->GetCurrentSource()->GetAngDist();
  angDist->SetParticleMomentumDirection(G4ThreeVector(1,0,0));
  angDist->SetBeamSigmaInAngX(PencilBeamStdAng);
  angDist->SetBeamSigmaInAngY(PencilBeamStdAng);  
  angDist->SetMinTheta(0.0);
  angDist->SetMaxTheta(1.0);

  PencilBeamStdY = 6.393*mm;
  PencilBeamStdZ = 6.5524*mm;
  particleSource->SetCurrentSourceIntensity(0.75);
  posDist = particleSource->GetCurrentSource()->GetPosDist();
  posDist->SetPosRot1(G4ThreeVector(0,0,1));
  posDist->SetPosDisType("Beam");
  posDist->SetPosDisShape("Circle");
  posDist->SetBeamSigmaInX(PencilBeamStdY);
  posDist->SetBeamSigmaInY(PencilBeamStdZ);  


  
  //Energy Spectrum
  eneDist = particleSource->GetCurrentSource()->GetEneDist();
  eneDist = particleSource->GetCurrentSource()->GetEneDist();
  eneDist->SetEnergyDisType("Arb");
  eneDist->ArbEnergyHistoFile("macro/source/protonKEspectrumAtZ-0.1cm_71_163.9MeV.dat");
  eneDist->ArbInterpolate("Lin");

  //Particle type
  particleSource->SetParticleDefinition(particle);

  ////////////////////////////////////////////////////////////////////////////////////
  //Source 2
  //Position

  particleSource->AddaSource(0.25);
  PencilBeamStdY = 12.8434*mm;
  PencilBeamStdZ = 13.0028*mm;
  posDist = particleSource->GetCurrentSource()->GetPosDist();
  posDist->SetPosRot1(G4ThreeVector(0,0,1));
  posDist->SetPosDisType("Beam");
  posDist->SetPosDisShape("Circle");
  posDist->SetBeamSigmaInX(PencilBeamStdY);
  posDist->SetBeamSigmaInY(PencilBeamStdZ);  
  
  //Angular Specturm
  PencilBeamStdAng = 25;
  angDist = particleSource->GetCurrentSource()->GetAngDist();
  angDist->SetParticleMomentumDirection(G4ThreeVector(1,0,0));
  angDist->SetBeamSigmaInAngY(PencilBeamStdAng);
  angDist->SetBeamSigmaInAngX(PencilBeamStdAng);  
  angDist->SetMinTheta(0.0);
  angDist->SetMaxTheta(1.0);
  
  //Energy Spectrum
  eneDist = particleSource->GetCurrentSource()->GetEneDist();
  eneDist->SetEnergyDisType("Arb");
  eneDist->ArbEnergyHistoFile("macro/source/protonKEspectrumAtZ-0.1cm_71_163.9MeV.dat");
  eneDist->ArbInterpolate("Lin");

  //Particle type
  particleSource->SetParticleDefinition(particle);

  nProtonsGenerated = 0;
  fieldSizeY = 2*(theDetector->ScintHalfY);
  fieldSizeZ = 2*(theDetector->ScintHalfZ);  
  
  //Pencil beam parameters
  PencilBeamPosY = linspace(-fieldSizeY/2, fieldSizeY/2, NPBY); // mm
  PencilBeamPosZ = linspace(-fieldSizeZ/2, fieldSizeZ/2, NPBZ); // mm
  for(int i =0; i<NPBY; i++){
    cout<<PencilBeamPosY[i]<<endl;
  }
    
}
  
PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  theGenerator = NULL;
  delete particleGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  //Pencil beam
  ProtonsPerPB = int(Nprotons/(NPBY*NPBZ));
  if(ProtonsPerPB == 0) ProtonsPerPB = 1; // Not enough protons to feed all PB, we set it to 1  
  idPBGlobal   = G4UniformRand()*NPBY*NPBZ; //5050
  idPBY        = (idPBGlobal - (idPBGlobal%NPBY))/NPBY;
  idPBZ        = idPBGlobal%NPBY;
  x0           = -1*theDetector->PhantomHalfX - theDetector->midX -1*mm -10*mm;
  y0           = PencilBeamPosY[idPBY]; 
  z0           = PencilBeamPosZ[idPBZ];

  // Set Source Position
  particleSource->SetCurrentSourceto(0);
  particleSource->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(x0,y0,z0));
 
  particleSource->SetCurrentSourceto(1);
  particleSource->GetCurrentSource()->GetPosDist()->SetCentreCoords(G4ThreeVector(x0,y0,z0));
 
  //particleSource->SetVerbosity(10);
  //cout<<"  "<<endl;
  particleSource->GeneratePrimaryVertex(anEvent);
  Einit = particleSource->GetCurrentSource()->GetParticleEnergy();
  nProtonsGenerated++;
  if(nProtonsGenerated%20000==0) cout << nProtonsGenerated << endl;
  
}

vector<G4double> PrimaryGeneratorAction::linspace(double start, double stop, int NStep) {
  std::vector<double> array;

  if (NStep == 0) return array; 
  if (NStep == 1)
    {
      array.push_back(start);
      return array;
    }
  double delta = (stop - start) / (NStep);
  for(int i=0; i < NStep; ++i)  array.push_back(start + delta * i);
  array.push_back(stop);
  return array;
}



