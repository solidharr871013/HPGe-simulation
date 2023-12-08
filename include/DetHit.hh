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
/// \file DetHit.hh
/// \brief Definition of the DetHit class

#ifndef DetHit_h
#define DetHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class G4AttDef;
class G4AttValue;

/// Drift chamber hit
///
/// It records:
/// - the layer ID
/// - the particle time
/// - the particle local and global positions

class DetHit : public G4VHit
{
  public:
    DetHit();
    DetHit(const DetHit &right);
    virtual ~DetHit();

    const DetHit& operator=(const DetHit &right);
    G4bool operator==(const DetHit &right) const;
    
    inline void *operator new(size_t);
    inline void operator delete(void *aHit);
    
    virtual void Draw();
    virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
    virtual std::vector<G4AttValue>* CreateAttValues() const;
    virtual void Print();

    void SetEdep(G4double edep) { fEdep = edep; }
    G4double GetEdep() const { return fEdep; }

    void SetScatteringPos(G4ThreeVector xyz) { fScatteringPos = xyz; }
    G4ThreeVector GetScatteringPos() const { return fScatteringPos; }

    void SetTrackID(G4int PID){fTrackID = PID;}
    G4int GetTrackID() const {return  fTrackID;}
    
  private:
    G4double fEdep;
    G4ThreeVector fScatteringPos;
    G4int fTrackID;

};

using HitsCollection = G4THitsCollection<DetHit>;

extern G4ThreadLocal G4Allocator<DetHit>* HitAllocator;

inline void* DetHit::operator new(size_t)
{
  if (!HitAllocator) {
      HitAllocator = new G4Allocator<DetHit>;
  }
  return (void*)  HitAllocator->MallocSingle();
}

inline void DetHit::operator delete(void* aHit)
{
    HitAllocator->FreeSingle((DetHit*) aHit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
