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
//

#include "globals.hh"

#include "ModularPhysicsList.hh"

#include "G4EmPenelopePhysics.hh"
#include "G4EmLivermorePhysics.hh"

#include "G4SystemOfUnits.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"

#include "G4PhotoElectricEffect.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4RayleighScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"

#include "G4PenelopeComptonModel.hh"
#include "G4PenelopeRayleighModel.hh"
#include "G4PenelopePhotoElectricModel.hh"
#include "G4PenelopeGammaConversionModel.hh"
#include "G4PenelopeIonisationModel.hh"
#include "G4PenelopeBremsstrahlungModel.hh"
#include "G4PenelopeAnnihilationModel.hh"

#include "G4LivermoreComptonModel.hh"
#include "G4LivermoreRayleighModel.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4LivermoreGammaConversionModel.hh"
#include "G4LivermoreIonisationModel.hh"
#include "G4LivermoreBremsstrahlungModel.hh"

#include "G4eMultipleScattering.hh"
#include "G4hMultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"
#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"

#include "G4hImpactIonisation.hh"

#include "G4ProductionCutsTable.hh"


//  PhysicsList:: PhysicsList():  G4VUserPhysicsList()
// {
//   defaultCutValue = 0.001*mm;
//    SetVerboseLevel(1);

//    std::cout << "==============================================================================="
// 	     << std::endl
// 	     << "Geant4   example - based on a simplified version of   simulation"
// 	     << std::endl
// 	     << "Further details can be found in:"
// 	     << std::endl
// 	     << "M.G. Pia et al., 'PIXE Simulation With Geant4', "
// 	     << "IEEE Trans. Nucl. Sci., vol. 56, no. 6, pp. 3614-3649, 2009"
// 	     << std::endl
// 	     << "N. Meidinger et al., 'Development of the focal plane PNCCD camera system for the X-ray space telescope  ', " 
// 	     << std::endl
// 	     <<"NIM A 624, 321-329, 2010"
// 	     << std::endl
// 	     << "==============================================================================="
// 	     << std::endl;

//    std::cout<< std::endl;
   
//    std::cout << "==============================================================================="
// 	     << std::endl
// 	     << " The use of G4LowEnergyIonisation, G4LowEnergyBremsstrahlung, "
// 	     << std::endl
// 	     << "G4LowEnergyPhotoElectric, G4LowEnergyCompton, G4LowEnergyGammaConversion"
// 	     << std::endl
// 	     << "in this example is intentional. These classes will be replaced by other classes"
// 	     << std::endl
// 	     << "appropriate to the problem domain in a forthcoming Geant4 version"
// 	     << std::endl
// 	     << "==============================================================================="
// 	     << std::endl;
// }


//  PhysicsList::~ PhysicsList()
// {}


// void  PhysicsList::ConstructParticle()
// {
//   ConstructBosons();
//   ConstructLeptons();
//   ConstructMesons();
//   ConstructBaryons();
// }


// void  PhysicsList::ConstructBosons()
// {
//   // pseudo-particles
//   //G4Geantino::GeantinoDefinition();
//   //G4ChargedGeantino::ChargedGeantinoDefinition();

//   // gamma
//   G4Gamma::GammaDefinition();
// }


// void  PhysicsList::ConstructLeptons()
// {
//   // leptons
//   //  e+/-
//   G4Electron::ElectronDefinition();
//   G4Positron::PositronDefinition();
//   // mu+/-
//   //G4MuonPlus::MuonPlusDefinition();
//   //G4MuonMinus::MuonMinusDefinition();
//   // nu_e
//   //G4NeutrinoE::NeutrinoEDefinition();
//   //G4AntiNeutrinoE::AntiNeutrinoEDefinition();
//   // nu_mu
//   //G4NeutrinoMu::NeutrinoMuDefinition();
//   //G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
// }


// void  PhysicsList::ConstructMesons()
// {
//   //  mesons
//   //    light mesons
//   //G4PionPlus::PionPlusDefinition();
//   //G4PionMinus::PionMinusDefinition();
//   //G4PionZero::PionZeroDefinition();
//   //G4Eta::EtaDefinition();
//   //G4EtaPrime::EtaPrimeDefinition();
//   //G4KaonPlus::KaonPlusDefinition();
//   //G4KaonMinus::KaonMinusDefinition();
//   //G4KaonZero::KaonZeroDefinition();
//   //G4AntiKaonZero::AntiKaonZeroDefinition();
//   //G4KaonZeroLong::KaonZeroLongDefinition();
//   //G4KaonZeroShort::KaonZeroShortDefinition();
// }


// void  PhysicsList::ConstructBaryons()
// {
//   //  barions
//   G4Proton::ProtonDefinition();
//   G4AntiProton::AntiProtonDefinition();

//   //G4Neutron::NeutronDefinition();
//   //G4AntiNeutron::AntiNeutronDefinition();
// }


// void  PhysicsList::ConstructProcess()
// {
//   AddTransportation();
//   ConstructEM();
//   ConstructGeneral();
//   //AddStepMax();
// }



// void  PhysicsList::ConstructEM()
// {
//   auto theParticleIterator=GetParticleIterator();
//   theParticleIterator->reset();
//   while( (*theParticleIterator)() ){
//     G4ParticleDefinition* particle = theParticleIterator->value();
//     G4ProcessManager* processManager = particle->GetProcessManager();
//     G4String particleName = particle->GetParticleName();
     
//     if (particleName == "gamma") {

//       // photon   

//       G4PhotoElectricEffect* photoelectric = new G4PhotoElectricEffect;
//       //photoelectric->SetEmModel(new G4PenelopePhotoElectricModel);
//       photoelectric->SetEmModel(new G4LivermorePhotoElectricModel);
//       //photoelectric->ActivateAuger(true);
//       //photoelectric->SetCutForLowEnSecPhotons(0.250 * keV);
//       //photoelectric->SetCutForLowEnSecElectrons(0.250 * keV);
//       G4ComptonScattering* compton = new G4ComptonScattering;
//       //compton->SetEmModel(new G4PenelopeComptonModel);
//       compton->SetEmModel(new G4LivermoreComptonModel);

//       G4GammaConversion* gammaConversion = new G4GammaConversion;
//       //gammaConversion->SetEmModel(new G4PenelopeGammaConversionModel);
//       gammaConversion->SetEmModel(new G4LivermoreGammaConversionModel);

//       G4RayleighScattering* rayleigh = new G4RayleighScattering;
//       //rayleigh->SetEmModel(new G4PenelopeRayleighModel);
//       rayleigh->SetEmModel(new G4LivermoreRayleighModel);


//       processManager -> AddDiscreteProcess(photoelectric);
//       processManager -> AddDiscreteProcess(compton);
//       processManager -> AddDiscreteProcess(gammaConversion);
//       processManager -> AddDiscreteProcess(rayleigh);

    
//     } else if (particleName == "e-") {

//       // electron

//       G4eMultipleScattering* eMultipleScattering = new G4eMultipleScattering();
      

//       G4eIonisation* eIonisation = new G4eIonisation();
//       //eIonisation->SetEmModel(new G4PenelopeIonisationModel);
//       eIonisation->SetEmModel(new G4LivermoreIonisationModel);

//       G4eBremsstrahlung* eBremsstrahlung = new G4eBremsstrahlung();
//       //eBremsstrahlung->SetEmModel(new G4PenelopeBremsstrahlungModel);
//       eBremsstrahlung->SetEmModel(new G4LivermoreBremsstrahlungModel);

//       processManager -> AddProcess(eMultipleScattering, -1, 1, 1);
//       processManager -> AddProcess(eIonisation, -1, 2, 2);
//       processManager -> AddProcess(eBremsstrahlung, -1, -1, 3);   

//     } else if (particleName == "e+") {
//       // positron     
//       G4eIonisation* eplusIonisation = new G4eIonisation();
//       eplusIonisation->SetEmModel(new G4PenelopeIonisationModel);
//       //eplusIonisation->SetEmModel(new G4LivermoreIonisationModel);

//       G4eBremsstrahlung* eplusBremsstrahlung = new G4eBremsstrahlung();
//       //eplusBremsstrahlung->SetEmModel(new G4PenelopeBremsstrahlungModel);
//       eplusBremsstrahlung->SetEmModel(new G4LivermoreBremsstrahlungModel);

//       G4eplusAnnihilation* eplusAnnihilation = new G4eplusAnnihilation();
//       eplusAnnihilation->SetEmModel(new G4PenelopeAnnihilationModel);

//       processManager->AddProcess(new G4eMultipleScattering, -1, 1, 1);
//       processManager->AddProcess(eplusIonisation,         -1, 2, 2);
//       processManager->AddProcess(eplusBremsstrahlung,     -1, 3, 3);
//       processManager->AddProcess(eplusAnnihilation,    0,-1, 4);

//       //} else if( particleName == "mu+" || 
//       //        particleName == "mu-"    ) {
//       //muon  
//       //processManager->AddProcess(new G4MuMultipleScattering, -1, 1, 1);
//       //processManager->AddProcess(new G4MuIonisation,         -1, 2, 2);
//       //processManager->AddProcess(new G4MuBremsstrahlung,     -1, 3, 3);
//       //processManager->AddProcess(new G4MuPairProduction,     -1, 4, 4);       
             
//     } else if( particleName == "proton" ||
//                particleName == "pi-" ||
//                particleName == "pi+"    ) {
//       //proton  
//       /*
//       G4hImpactIonisation* hIonisation = new G4hImpactIonisation();
//       hIonisation->SetPixeCrossSectionK("ecpssr");
//       hIonisation->SetPixeCrossSectionL("ecpssr");
//       hIonisation->SetPixeCrossSectionM("ecpssr");
//       hIonisation->SetPixeProjectileMinEnergy(1.* keV);
//       hIonisation->SetPixeProjectileMaxEnergy(200. * MeV);
//       hIonisation->SetCutForSecondaryPhotons(250. * eV);
//       hIonisation->SetCutForAugerElectrons(250. * eV);
//       */
//       G4hIonisation* hIonisation = new G4hIonisation();

//       G4hMultipleScattering* hMultipleScattering = new G4hMultipleScattering();

//       processManager -> AddProcess(hMultipleScattering, -1, 1, 1);   
//       processManager -> AddProcess(hIonisation, -1, 2, 2);
     
//     } else if( particleName == "alpha" || 
// 	       particleName == "He3" ||
// 	       particleName == "pi-" ||
//                particleName == "pi+" ||
// 	       particleName == "GenericIon" ) {

//       // pions, alpha, ions (should never occur in the current example) 
//       processManager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
//       processManager->AddProcess(new G4ionIonisation,       -1, 2, 2);
                      
//     } else if ((!particle->IsShortLived()) &&
// 	       (particle->GetPDGCharge() != 0.0) && 
// 	       (particle->GetParticleName() != "chargedgeantino")) {
//       //all others charged particles except geantino
//       processManager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
//       processManager->AddProcess(new G4hIonisation,         -1, 2, 2);
//     }
//   }
// }

// #include "G4Decay.hh"

// void  PhysicsList::ConstructGeneral()
// {
//   // Add Decay Process
//   G4Decay* theDecayProcess = new G4Decay();
//   auto theParticleIterator=GetParticleIterator();
//   theParticleIterator->reset();
//   while( (*theParticleIterator)() ){
//     G4ParticleDefinition* particle = theParticleIterator->value();
//     G4ProcessManager* processManager = particle->GetProcessManager();
//     if (theDecayProcess->IsApplicable(*particle)) { 
//       processManager ->AddProcess(theDecayProcess);
//       // set ordering for PostStepDoIt and AtRestDoIt
//       processManager ->SetProcessOrdering(theDecayProcess, idxPostStep);
//       processManager ->SetProcessOrdering(theDecayProcess, idxAtRest);
//     }
//   }
// }
  

// /*
// #include "G4StepLimiter.hh"
// #include "G4UserSpecialCuts.hh"

// void  PhysicsList::AddStepMax()
// {
//   // Step limitation seen as a process
//   G4StepLimiter* stepLimiter = new G4StepLimiter();
//   ////G4UserSpecialCuts* userCuts = new G4UserSpecialCuts();
  
//   theParticleIterator->reset();
//   while ((*theParticleIterator)()){
//       G4ParticleDefinition* particle = theParticleIterator->value();
//       G4ProcessManager* processManager = particle->GetProcessManager();

//       if (particle->GetPDGCharge() != 0.0)
//         {
// 	  processManager ->AddDiscreteProcess(stepLimiter);
// 	  ////processManager ->AddDiscreteProcess(userCuts);
//         }
//   }
// }
// */

// void  PhysicsList::SetCuts()
// {
//   //G4VUserPhysicsList::SetCutsWithDefault method sets 
//   //the default cut value for all particle types 
//   //
//   SetCutsWithDefault();

//   // Set the secondary production cut lower than 990. eV
//   // Very important for processes at low energies
 
//   G4double lowLimit = 250. * eV;
//   G4double highLimit = 100. * GeV;
//   G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowLimit, highLimit);
     
//   if (verboseLevel>0) DumpCutValuesTable();
// }

#include "ModularPhysicsList.hh"

#include "G4DecayPhysics.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"

 ModularPhysicsList::ModularPhysicsList()
: G4VModularPhysicsList(){
  SetVerboseLevel(1);

  // Default physics
  RegisterPhysics(new G4DecayPhysics());

  // EM physics
  RegisterPhysics(new G4EmPenelopePhysics());
  //RegisterPhysics(new G4EmLivermorePhysics());
  //RegisterPhysics(new G4EmStandardPhysics());

  // Radioactive decay
  RegisterPhysics(new G4RadioactiveDecayPhysics());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

 ModularPhysicsList::~ModularPhysicsList()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ModularPhysicsList::SetCuts()
{
  G4VUserPhysicsList::SetCuts();
    //SetCutValue(100*mm, "proton");
  SetCutValue(1*mm, "e-");
  SetCutValue(1*mm, "e+");
  SetCutValue(1*mm, "gamma");  
}
