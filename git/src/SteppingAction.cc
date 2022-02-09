#include <vector>
#include <iostream>
#include <iomanip>
#include "Analysis.hh"
#include "G4EventManager.hh"
#include "G4EnergyLossTables.hh"
#include "G4TrackVector.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4FastStep.hh"
#include "G4ParticleChange.hh"
#include "G4DynamicParticle.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4SteppingManager.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include "iomanip"
#include "G4UImanager.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "SteppingAction.hh"
#include "G4RunManager.hh"
class DetectorConstruction;
class EventAction;
SteppingAction *SteppingAction::theSteppingAction=NULL;
SteppingAction::~SteppingAction()
{theSteppingAction=NULL;}
SteppingAction::SteppingAction()
{
  theAnalysis = Analysis::GetInstance();
  theSteppingAction=this;
}

void SteppingAction::UserSteppingAction(const G4Step*)
{
  //if ( aStep->GetTrack()->GetNextVolume()->GetName()=="physWorld" ) aStep->GetTrack()->SetTrackStatus(fStopAndKill);
}


