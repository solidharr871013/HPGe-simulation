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
/// \file SD.cc
/// \brief Implementation of the class

#include "SD.hh"
#include "DetHit.hh"

#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include "G4EventManager.hh"
#include "G4SteppingManager.hh"

#include "random"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  SD::SD(G4String name)
: G4VSensitiveDetector(name), 
  fHitsCollection(nullptr), fHCID(-1)
{
  collectionName.insert("recorderColl");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

  SD::~SD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SD::Initialize(G4HCofThisEvent* hce)
{
  fHitsCollection 
    = new HitsCollection(SensitiveDetectorName,collectionName[0]);

  if (fHCID<0) { 
     fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(fHitsCollection); 
  }
  hce->AddHitsCollection(fHCID,fHitsCollection);

 // G4cout << "the hit collection ID of   is " << fHCID << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool SD::ProcessHits(G4Step* step, G4TouchableHistory*)
{

  auto track = step->GetTrack();

  //if(track->GetParticleDefinition()->GetParticleName()=="gamma"){

      //auto preStepPoint = step->GetPreStepPoint();
      auto postStepPoint = step->GetPostStepPoint();

      auto worldPos = postStepPoint->GetPosition();

      //G4double Edep = track->GetTotalEnergy();
      G4double Edep = step->GetTotalEnergyDeposit();
      auto hit = new DetHit;
      hit->SetScatteringPos(worldPos);
      hit->SetEdep(Edep);
      hit->SetTrackID(track->GetTrackID());
      //hit->SetTime(preStepPoint->GetGlobalTime());


      fHitsCollection->insert(hit);


 //}

  

  
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
