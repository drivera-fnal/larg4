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
// $Id: G4EmExtraPhysics.hh 66704 2013-01-10 18:20:17Z gunter $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4EmExtraPhysics
//
// Author: 2002 J.P. Wellisch
//
// Modified:
// 10.11.2005 V.Ivanchenko edit to provide a standard
// 19.06.2006 V.Ivanchenko add mu-nuclear process
// 16.10.2012 A.Ribon: renamed G4EmExtraBertiniPhysics as G4EmExtraPhysics
// 04.04.2020 D. Rivera -- copy based on :
//    geant4.10.03.p03/source/physics_lists/constructors/gamma_lepto_nuclear/include/G4EmExtraPhysics.hh
//
//----------------------------------------------------------------------------
//

#ifndef MyG4EmExtraPhysics_h
#define MyG4EmExtraPhysics_h 1

#include "Geant4/G4VPhysicsConstructor.hh"
#include "Geant4/globals.hh"

//#include "Geant4/G4EmMessenger.hh"
#include "MyG4EmMessenger.hh"

//<--class G4BertiniElectroNuclearBuilder;
class MyG4BertiniElectroNuclearBuilder;
class G4SynchrotronRadiation;
class G4GammaConversionToMuons;
class G4AnnihiToMuPair;
class G4eeToHadrons;

class MyG4EmExtraPhysics : public G4VPhysicsConstructor
{
public:

  MyG4EmExtraPhysics(G4int ver = 1);

  // obsolete
  MyG4EmExtraPhysics(const G4String& name);

  virtual ~MyG4EmExtraPhysics();

  void ConstructParticle();
  void ConstructProcess();

  void Synch(G4bool val);
  void SynchAll(G4bool val);
  void GammaNuclear(G4bool val);
  void MuonNuclear(G4bool val);
  void GammaToMuMu(G4bool val);
  void PositronToMuMu(G4bool val);
  void PositronToHadrons(G4bool val);

private:

  static G4bool gnActivated;
  static G4bool munActivated;
  static G4bool synActivated;
  static G4bool synActivatedForAll;
  static G4bool gmumuActivated;
  static G4bool pmumuActivated;
  static G4bool phadActivated;

  static G4ThreadLocal MyG4BertiniElectroNuclearBuilder* theGNPhysics;
  static G4ThreadLocal G4SynchrotronRadiation* theSynchRad;
  static G4ThreadLocal G4GammaConversionToMuons* theGammaToMuMu;
  static G4ThreadLocal G4AnnihiToMuPair* thePosiToMuMu;
  static G4ThreadLocal G4eeToHadrons* thePosiToHadrons;

  //<--G4EmMessenger* theMessenger;
  MyG4EmMessenger* theMessenger;
  G4int verbose;
};

#endif





