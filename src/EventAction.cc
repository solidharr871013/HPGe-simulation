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
/// \file  EventAction.cc
/// \brief Implementation of the  EventAction class

#include "EventAction.hh"
#include "RunAction.hh"
#include "G4AnalysisManager.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4Track.hh"
#include "G4RunManager.hh"
#include "G4THitsMap.hh"

#include "DetHit.hh"
#include "SD.hh"

#include "G4SDManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4GenericMessenger.hh"

#include "random"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace {

// Utility function which finds a hit collection with the given Id
// and print warnings if not found
G4VHitsCollection* GetHC(const G4Event* event, G4int collId) {
  auto hce = event->GetHCofThisEvent();
  if (!hce) {
      G4ExceptionDescription msg;
      msg << "No hits collection of this event found." << G4endl;
      G4Exception(" EventAction::EndOfEventAction()",
                  " Code001", JustWarning, msg);
      return nullptr;
  }

  auto hc = hce->GetHC(collId);
  if ( ! hc) {
    G4ExceptionDescription msg;
    msg << "Hits collection " << collId << " of this event not found." << G4endl;
    G4Exception(" EventAction::EndOfEventAction()",
                " Code001", JustWarning, msg);
  }
  return hc;
}
}

 EventAction:: EventAction( RunAction* runAction)
: G4UserEventAction(),
  fMessenger(nullptr),
  fRunAction(runAction),
  fEdep(0.),
  f3LayerCoincidence(0),
  fTotalNumber(0),
  fOrdered3(0),
  fOrdered2(0),
  f2LayerCoincidence(0),
  fRecorderID(-1),
  fCollID_detectorEnergy(-1)
{
    DefineCommands();
     G4RunManager::GetRunManager()->SetPrintProgress(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

 EventAction::~ EventAction()
{ delete fMessenger; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void  EventAction::BeginOfEventAction(const G4Event*)
{    
  fEdep = 0.;
  f3LayerCoincidence = 0;
  fTotalNumber = 0;
  fOrdered3 = 0;
  fOrdered2 = 0;
  f2LayerCoincidence = 0;

  if(fRecorderID == -1 ){
      auto sdManager = G4SDManager::GetSDMpointer();

      fRecorderID = sdManager->GetCollectionID("recorder/recorderColl");
  }


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void  EventAction::EndOfEventAction(const G4Event* event)
{   
    G4AnalysisManager* analysis = G4AnalysisManager::Instance();


    auto HCE = event->GetHCofThisEvent();
    if(!HCE) return;

    if(fCollID_detectorEnergy<0){
        auto SDMan = G4SDManager::GetSDMpointer();
        fCollID_detectorEnergy = SDMan->GetCollectionID("detectorEnergy/energyDeposition");
    }

///////////////////test/////////////////////

    auto spectrumSys = GetHC(event, fRecorderID);
    if(!spectrumSys) return;

    if(spectrumSys->GetSize()!=0){

        for(unsigned long itr = 0; itr != spectrumSys->GetSize(); ++itr){

            auto hit = static_cast<DetHit*>(spectrumSys->GetHit(itr));
            G4double Edep = hit->GetEdep();
            fEdep += Edep;
        }

        if(fEdep>30*eV){
            analysis->FillH1(0,fEdep);
        }

    }
    
///////////////////////////////
    // G4double eThreshold = 30*eV;

    // G4THitsMap<G4double>* evtMap =
    //                    static_cast<G4THitsMap<G4double>*>(HCE->GetHC(fCollID_detectorEnergy));

    // std::map<G4int,G4double*>::iterator itr;
    // for(itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); ++itr){

    //     //G4double CopyNo = static_cast<G4double>(itr->first);
    //     G4double eDep = *(itr->second);

    //     fEdep += eDep;
    // }

    // if(fEdep>eThreshold){
    //         analysis->FillH1(0,fEdep);
    //        // G4cout << "CopyNo is " << CopyNo << G4endl;
    //     }  

    /*************************
    auto hcRecorder = GetHC(event,fRecorderID);
    if(!hcRecorder) return;

    if(hcRecorder->GetSize()!=0){
        ++fTotalNumber;
    }

    G4int firstTID = -1;
    for(unsigned long itr=0; itr!=hcRecorder->GetSize(); ++itr){
        auto hit = static_cast< Layer1Hit*>(hcRecorder->GetHit(itr));
        G4int presentTID = hit->GetTrackID();
        if(firstTID!=presentTID){
            firstTID = presentTID;
            analysis->FillH1(0,hit->GetEdep());
        }
        else{}
    }

    fRunAction->AddTotalNumber(fTotalNumber);
***********************/

}

G4double  EventAction::AddFluction(G4double val){

    //calculate the coefficient of the resolution, based on 662keV
    G4double Resolution662 =0.08,
             coefficient = Resolution662*std::sqrt(0.662);

    G4double valVar = coefficient*std::sqrt(val)/(2*std::sqrt(2*std::log(2)));
    //G4double valVar = 0.08*val/(2*std::sqrt(2*std::log(2)));


    std::random_device rd;
    std::mt19937 generator(rd());

    std::normal_distribution<G4double> distribution(val,valVar);
    G4double var_new = distribution(generator);

    return var_new;
}


void  EventAction::Set2TotalELow(G4double val){f2TotalE_low = val;}
void  EventAction::Set2TotalEHigh(G4double val){f2TotalE_high = val;}
void  EventAction::Set2ScatterELow(G4double val){f2ScatterE_low = val;}
void  EventAction::Set2ScatterEHigh(G4double val){f2ScatterE_high = val;}

void  EventAction::DefineCommands(){
    fMessenger = new G4GenericMessenger(this,
                                        "/ Camera/event/",
                                        "Event Control");


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
