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
/// \file  RunAction.cc
/// \brief Implementation of the  RunAction class

#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "G4AnalysisManager.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

 RunAction:: RunAction()
: G4UserRunAction(),
  fEdep(0.),
  fEdep2(0.),
  f3LayerCoincidence(0),
  fTotalNumber(0),
  fOrdered3(0),
  fOrdered2(0),
  f2LayerCoincidence(0)
{ 

  // Register accumulable to the accumulable manager
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->RegisterAccumulable(fEdep);
  accumulableManager->RegisterAccumulable(fEdep2); 
  accumulableManager->RegisterAccumulable(f3LayerCoincidence);
  accumulableManager->RegisterAccumulable(fTotalNumber);
  accumulableManager->RegisterAccumulable(fOrdered3);
  accumulableManager->RegisterAccumulable(fOrdered2);
  accumulableManager->RegisterAccumulable(f2LayerCoincidence);

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  analysisManager->SetDefaultFileType("csv");

  //analysisManager->SetNtupleMerging(true);
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetFileName("target");

  analysisManager->CreateH1("energy","energySpectrum",2999,0,3000*keV);//ID=0
  //analysisManager->CreateH1("source-energy","energySpectrum2",2999,0,3000*keV);//ID=1

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

 RunAction::~ RunAction()
{delete G4AnalysisManager::Instance();}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void  RunAction::BeginOfRunAction(const G4Run*)
{ 
  // inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

  // reset accumulables to their initial values
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Reset();

  G4AnalysisManager* analysis = G4AnalysisManager::Instance();

  analysis->Reset();
  analysis->OpenFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void  RunAction::EndOfRunAction(const G4Run* run)
{
  G4int nofEvents = run->GetNumberOfEvent();
  if (nofEvents == 0) return;

  G4double edepTest = fEdep.GetValue();
  G4cout << "the test energy value is " << edepTest/keV << G4endl;

  // Merge accumulables 
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();

  // Compute dose = total energy deposit in a run and its variance
  //
  G4double edep  = fEdep.GetValue();
  G4double edep2 = fEdep2.GetValue();
  // cast int to double type
  G4double NO3LayerCoincidence = f3LayerCoincidence.GetValue();
  G4double NOTotalEvent = fTotalNumber.GetValue();
  G4double NoOrdered3 = fOrdered3.GetValue();
  G4double NoOrdered2 = fOrdered2.GetValue();
  G4double No2LayerCoincidence = f2LayerCoincidence.GetValue();

  G4double rms = edep2 - edep*edep/nofEvents;
  if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;  

  const  DetectorConstruction* detectorConstruction
   = static_cast<const  DetectorConstruction*>
     (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  //G4double mass = detectorConstruction->GetLayer2LogicalVolume()->GetMass();
  //G4double dose = edep/mass;
  //G4double rmsDose = rms/mass;

  G4double  Event3Ratio = NO3LayerCoincidence/NOTotalEvent,
            Event2Ratio = No2LayerCoincidence/NOTotalEvent;

  // Run conditions
  //  note: There is no primary generator action object for "master"
  //        run manager for multi-threaded mode.
  const  PrimaryGeneratorAction* generatorAction
   = static_cast<const  PrimaryGeneratorAction*>
     (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  G4String runCondition;
  // if (generatorAction)
  // {
  //   //const G4ParticleGun* particleGun = generatorAction->GetParticleGun();
  //   //runCondition += particleGun->GetParticleDefinition()->GetParticleName();
  //   //runCondition += " of ";
  //   //G4double particleEnergy = particleGun->GetParticleEnergy();
  //   //runCondition += G4BestUnit(particleEnergy,"Energy");
  // }

  //G4double thicknessLayer1 = detectorConstruction->GetLayer1Thickness();


  // Print
  //  
  if (IsMaster()) {
    G4cout
     << G4endl
     << "--------------------End of Global Run-----------------------";
  }
  else {
    G4cout
     << G4endl
     << "--------------------End of Local Run------------------------";
  }
  
  G4cout
     << G4endl
     //<< " The thickness of layer1 is " << G4BestUnit(thicknessLayer1, "Length") << G4endl
     << " The run consists of " << nofEvents << " "<< runCondition
     << G4endl
     //<< " Cumulated dose per run, in scoring volume : "
     //<< G4BestUnit(dose,"Dose") << " rms = " << G4BestUnit(rmsDose,"Dose")
     //<< G4endl
     << "The number of total event is " << NOTotalEvent << G4endl
     << "The number of 3LayerCoincidence event is " << NO3LayerCoincidence << G4endl
     << "The number of ordered 3LayerCoincident event is " << NoOrdered3 << G4endl
     << "The number of 2LayerCoincidence event is " << No2LayerCoincidence << G4endl
     << "The number of ordered 2LayerCoincident event is " << NoOrdered2 << G4endl
     << "The 3   event ratio is " << std::setw(3) <<  Event3Ratio << G4endl
     << "The 2   event ratio is " << std::setw(3) <<  Event2Ratio << G4endl
     << "The total   event ratio is " << std::setw(3)
     <<  Event3Ratio+ Event2Ratio
     << G4endl
     << "------------------------------------------------------------"
     << G4endl
     << G4endl;

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  analysisManager->Write();
  analysisManager->CloseFile();


/****************Start txt format output***************

  std::ofstream Single Efficiency;


    Single Efficiency.open(fOutput + "Single Efficiency.txt", std::ios_base::app);
    if (Single Efficiency.is_open()){

      Single Efficiency << "The thickness is " << G4BestUnit(thicknessLayer1, "Length") << "; "
                              << "the efficiency is" << std::setw(4) <<  EventRatio << ";"
                              << std::endl;
      }
    else Single Efficiency << "Unable to open file" << std::endl;

****************End txt format output*********************/

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void  RunAction::AddEdep(G4double edep)
{
  fEdep  += edep;
  fEdep2 += edep*edep;
}

void  RunAction::Add2LayerCoincidence(G4int number){
    f2LayerCoincidence += number;
}

void  RunAction::Add3LayerCoincidence(G4int number){
    f3LayerCoincidence += number;
}

void  RunAction::AddOrdered2(G4int number){
    fOrdered2 += number;
}

void  RunAction::AddOrdered3(G4int number){
    fOrdered3 += number;
}

void  RunAction::AddTotalNumber(G4int totalnumber){
    fTotalNumber += totalnumber;
}





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

