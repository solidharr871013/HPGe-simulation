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
/// \file  DetectorConstruction.cc
/// \brief Implementation of the  DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
//#include "G4Cons.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4SubtractionSolid.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4VisAttributes.hh"
#include "G4GenericMessenger.hh"

#include "G4PhysicalConstants.hh"

#include "SD.hh"

#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"

#include "G4Region.hh"
#include "G4ProductionCuts.hh"

#include "G4GDMLParser.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

 DetectorConstruction:: DetectorConstruction()
: G4VUserDetectorConstruction(),
  fMessenger(nullptr)
{
    DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

 DetectorConstruction::~ DetectorConstruction()
{
    //delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume*  DetectorConstruction::Construct()
{
   //
   // Construct the material to be used
   //
   ConstructMaterial();

   G4Material* worldmat = G4Material::GetMaterial("G4_AIR");
   G4Material* copper = G4Material::GetMaterial("G4_Cu");
   G4Material* aluminum = G4Material::GetMaterial("G4_Al");
    G4Material* germanium = G4Material::GetMaterial("G4_Ge");
    G4Material* mylar = G4Material::GetMaterial("G4_MYLAR");
    G4Material* PMMA = G4Material::GetMaterial("PMMA");
    G4Material* vacuum = G4Material::GetMaterial("Vacuum");
    G4Material* gold = G4Material::GetMaterial("G4_Au");
    G4Material* manganese = G4Material::GetMaterial("G4_Mn");
    G4Material* lead = G4Material::GetMaterial("G4_Pb");
    G4Material* lithium = G4Material::GetMaterial("G4_Li");


   G4GDMLParser parser;

   parser.SetOverlapCheck(true);
   parser.Read("HPGe-source15cm.gdml");
   //parser.Read("HPGe-holdertest.gdml");
   
   auto worldPhys = parser.GetWorldVolume();

   auto worldLog = worldPhys->GetLogicalVolume();
        worldLog->SetMaterial(worldmat);
 
   fRecorderLog = parser.GetVolume("V-HPGeSensitive");
   fRecorderLog->SetMaterial(germanium);
   
   if(parser.GetVolume("V-innerDeadLayer") != nullptr){
     auto innerDeadLayerLog = parser.GetVolume("V-innerDeadLayer");
     innerDeadLayerLog->SetMaterial(germanium);
   }

   if(parser.GetVolume("V-outterDeadLayer") != nullptr){
     auto outerDeadLayerLog = parser.GetVolume("V-outterDeadLayer");
     outerDeadLayerLog->SetMaterial(germanium);
   }
    
   if(parser.GetVolume("V-copperFinger") != nullptr){
     auto copperFingerLog = parser.GetVolume("V-copperFinger");
     copperFingerLog->SetMaterial(copper);
   }
   
        

   auto topLayerAlLog = parser.GetVolume("V-topLayerAl");
        topLayerAlLog->SetMaterial(aluminum);

   auto topLayerMylar = parser.GetVolume("V-topLayerMylar");
        topLayerMylar->SetMaterial(mylar);

   auto CupLog = parser.GetVolume("V-Cup");
        CupLog->SetMaterial(aluminum);

   auto shellLog = parser.GetVolume("V-shell");
        shellLog->SetMaterial(aluminum);

   auto sampleHolderLog = parser.GetVolume("V-sampleHolder");
        sampleHolderLog->SetMaterial(PMMA);

   auto PlasticShieldLog = parser.GetVolume("V-plasticShielding");
        PlasticShieldLog->SetMaterial(PMMA);

   auto CopperShieldLog = parser.GetVolume("V-copperShielding");
        CopperShieldLog->SetMaterial(copper);

   auto LeadShieldLog = parser.GetVolume("V-LeadShielding");
        LeadShieldLog->SetMaterial(lead);


   /////////////////
   //place a foil
   /////////////////

   G4double holderHeight = 147.4*mm;// 97.4mm, 147.4mm, 197.4mm
   G4double foilThickness = 0.0508*mm, foilRadius = 6.19*mm; //thickness Au-Mn-0.0508mm, radius Au-6.19mm, Mn-6.335mm

   G4ThreeVector foilPosition = G4ThreeVector(0,0,(6+holderHeight)*mm+0.5*foilThickness);
   G4Tubs* foil = new G4Tubs("foil",0,foilRadius,0.5*foilThickness,0,360*deg);
   fFoilLog = new G4LogicalVolume(foil,gold,"foilLog");
   new G4PVPlacement(nullptr,foilPosition,fFoilLog,"foilPhys",worldLog,false,0);

    ////////////////////////////
    // if the foil is copper
    ////////////////////////////

//     G4Tubs* foilCu = new G4Tubs("foilCu",0,foilRadius,0.5*foilThickness,0,360*deg);
//     G4LogicalVolume* foilCuLog = new G4LogicalVolume(foilCu,copper,"foilCuLog");

//     G4double alCoverThickness = 1*mm;
//     G4Tubs* AlCover = new G4Tubs("AlCover", 0, foilRadius+0.1*mm, 0.5*alCoverThickness, 0, 360*deg);
//     fFoilLog = new G4LogicalVolume(AlCover,aluminum,"AlCoverLog");

//    //6mm is the thickness of the sum of the stage and the holder stage, the original line is based on HPGe. 
//    G4ThreeVector AlCover1Position = G4ThreeVector(0,0,(6+holderHeight)*mm+0.5*alCoverThickness),
//                 AlCover2Position = G4ThreeVector(0,0,(6+holderHeight)*mm+1.5*alCoverThickness+foilThickness),
//                 foilPositionCu = G4ThreeVector(0,0,(6+holderHeight)*mm+alCoverThickness+0.5*foilThickness);

//      new G4PVPlacement(nullptr,AlCover1Position,fFoilLog,"AlCover1Phys",worldLog,false,0,true);
//      new G4PVPlacement(nullptr,AlCover2Position,fFoilLog,"AlCover2Phys",worldLog,false,0,true);
//      new G4PVPlacement(nullptr,foilPositionCu,foilCuLog,"foilCuPhys",worldLog,false,0,true);

   //fFoilLog = AlCoverLog;

   //
   // always return world physical volume
   //

   return worldPhys;


}

void  DetectorConstruction::ConstructMaterial(){
    // Get NistManager
    G4NistManager* man = G4NistManager::Instance();

    man->FindOrBuildMaterial("G4_AIR");
    man->FindOrBuildMaterial("G4_Cu");
    man->FindOrBuildMaterial("G4_Al");
    man->FindOrBuildMaterial("G4_Ge");
    man->FindOrBuildMaterial("G4_Au");
    man->FindOrBuildMaterial("G4_Mn");
    man->FindOrBuildMaterial("G4_MYLAR");
    man->FindOrBuildMaterial("G4_Pb");
    man->FindOrBuildMaterial("G4_Li");

    G4Element* C = man->FindOrBuildElement("C",true);
    G4Element* H = man->FindOrBuildElement("H",true);
    G4Element* O = man->FindOrBuildElement("O",true);
    // G4Element* B = man->FindOrBuildElement("B",true);
    // G4Element* Li = man->FindOrBuildElement("Li",true);

    G4Material* PMMA = new G4Material("PMMA",1.18*g/cm3,3);
    PMMA->AddElement(C,5);
    PMMA->AddElement(H,8);
    PMMA->AddElement(O,2);


    G4double atomicNumber = 1.;
    G4double massOfMole = 1.008*g/mole;
    G4double density = 1.e-25*g/cm3;
    G4double temperature = 2.73*kelvin;
    G4double pressure = 3.e-18*pascal;
    G4Material* Vacuum
            = new G4Material("Vacuum", atomicNumber,
                             massOfMole, density, kStateGas,
                             temperature, pressure);
    /***********************/

}

void  DetectorConstruction::ConstructSDandField(){

    auto sdManger = G4SDManager::GetSDMpointer();
    G4String SDname;

    auto recorder = new SD(SDname = "/recorder");
    sdManger->AddNewDetector(recorder);
    fRecorderLog->SetSensitiveDetector(recorder);

    auto detectorEnergy = new G4MultiFunctionalDetector("detectorEnergy");
    sdManger->AddNewDetector(detectorEnergy);
    G4VPrimitiveScorer* primitive = new G4PSEnergyDeposit("energyDeposition");
    detectorEnergy->RegisterPrimitive(primitive);
    SetSensitiveDetector(fRecorderLog,detectorEnergy);


}


void  DetectorConstruction::DefineCommands(){

    fMessenger = new G4GenericMessenger(this,
                                        "/ Camera/detector/",
                                        "Detector control");
/********************
    //  command setting the thickness of Layer1
    auto& Layer1ThicknessCmd
      = fMessenger->DeclareMethodWithUnit("Layer1Thickness","mm",
                                  & DetectorConstruction::SetLayer1Thickness,
                                  "Set thickness of the Layer1.");
    Layer1ThicknessCmd.SetParameterName("Thickness", true);
    Layer1ThicknessCmd.SetRange("Thickness>=0. && Thickness<100.");
    Layer1ThicknessCmd.SetDefaultValue("8.");

    // command setting the thickness of layer2
    auto& Layer2ThicknessCmd
      = fMessenger->DeclareMethodWithUnit("Layer2Thickness","mm",
                                  & DetectorConstruction::SetLayer2Thickness,
                                  "Set thickness of the Layer2.");
    Layer2ThicknessCmd.SetParameterName("Thickness", true);
    Layer2ThicknessCmd.SetRange("Thickness>=0. && Thickness<100.");
    Layer2ThicknessCmd.SetDefaultValue("8.");

    // command setting the thickness of layer3
    auto& Layer3ThicknessCmd
      = fMessenger->DeclareMethodWithUnit("Layer3Thickness","mm",
                                  & DetectorConstruction::SetLayer3Thickness,
                                  "Set thickness of the Layer3.");
    Layer3ThicknessCmd.SetParameterName("Thickness", true);
    Layer3ThicknessCmd.SetRange("Thickness>=0. && Thickness<100.");
    Layer3ThicknessCmd.SetDefaultValue("8.");

    auto& Cell1SizeCmd
      = fMessenger->DeclareMethodWithUnit("Cell1Size","mm",
                                  & DetectorConstruction::SetCell1Size,
                                  "Set size of the Cell1.");
    Cell1SizeCmd.SetParameterName("cellsize", true);
    Cell1SizeCmd.SetRange("cellsize>=0. && cellsize<100.");
    Cell1SizeCmd.SetDefaultValue("2.");

    auto& Cell2SizeCmd
      = fMessenger->DeclareMethodWithUnit("Cell2Size","mm",
                                  & DetectorConstruction::SetCell2Size,
                                  "Set size of the Cell2.");
    Cell2SizeCmd.SetParameterName("cellsize", true);
    Cell2SizeCmd.SetRange("cellsize>=0. && cellsize<100.");
    Cell2SizeCmd.SetDefaultValue("2.");

    auto& Cell3SizeCmd
      = fMessenger->DeclareMethodWithUnit("Cell3Size","mm",
                                  & DetectorConstruction::SetCell3Size,
                                  "Set size of the Cell3.");
    Cell3SizeCmd.SetParameterName("cellsize", true);
    Cell3SizeCmd.SetRange("cellsize>=0. && cellsize<100.");
    Cell3SizeCmd.SetDefaultValue("2.");

    auto& Layer1ZPositionCmd
      = fMessenger->DeclareMethodWithUnit("Layer1Z","mm",
                                  & DetectorConstruction::SetLayer1ZPosition,
                                  "Set z position of the Layer1.");
    Layer1ZPositionCmd.SetParameterName("ZPosition", true);
    Layer1ZPositionCmd.SetRange("ZPosition>=-100. && ZPosition<100.");
    Layer1ZPositionCmd.SetDefaultValue("-18");

    auto& Layer2ZPositionCmd
      = fMessenger->DeclareMethodWithUnit("Layer2Z","mm",
                                  & DetectorConstruction::SetLayer2ZPosition,
                                  "Set z position of the Layer2.");
    Layer2ZPositionCmd.SetParameterName("ZPosition", true);
    Layer2ZPositionCmd.SetRange("ZPosition>=-50. && ZPosition<50.");
    Layer2ZPositionCmd.SetDefaultValue("0");

    auto& Layer3ZPositionCmd
      = fMessenger->DeclareMethodWithUnit("Layer3Z","mm",
                                  & DetectorConstruction::SetLayer3ZPosition,
                                  "Set z position of the Layer3.");
    Layer3ZPositionCmd.SetParameterName("ZPosition", true);
    Layer3ZPositionCmd.SetRange("ZPosition>=-100. && ZPosition<100.");
    Layer3ZPositionCmd.SetDefaultValue("18");
*********************/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
