//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file  SteppingAction.cc
/// \brief Implementation of the  SteppingAction class

#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4AnalysisManager.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "random"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

 SteppingAction:: SteppingAction( EventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fScoringVolume1(nullptr)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

 SteppingAction::~ SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void  SteppingAction::UserSteppingAction(const G4Step* step)
{

    if ( !fScoringVolume1 ) {
   const  DetectorConstruction* detectorConstruction
      = static_cast<const  DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume1 = detectorConstruction->GetFoliLog();


  }


  // get volume of the current step
  G4LogicalVolume* volume 
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();

  auto track = step->GetTrack();
  auto status = step->GetPostStepPoint()->GetStepStatus();



  // check if we are in scoring volume
  G4bool transmit = (status <= fGeomBoundary);  
  if (transmit==false) return;
  if (status == fWorldBoundary) return; 

  // //G4cout << "ok!" << G4endl;


  // auto preLogicalVolume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
  auto postLogicalVolume = step->GetPostStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();

  G4AnalysisManager* analysis = G4AnalysisManager::Instance();

  const G4ParticleDefinition* particle = track->GetParticleDefinition();
  G4String name     = particle->GetParticleName();
  //G4double meanLife = particle->GetPDGLifeTime();
  G4double ekin     = track->GetKineticEnergy();

  if(postLogicalVolume==fScoringVolume1){
    if(name=="e+"){
      //G4cout << name << G4endl;
      //analysis->FillH1(1,ekin);
    }
  }







}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

