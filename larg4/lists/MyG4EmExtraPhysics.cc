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
// $Id: G4EmExtraPhysics.cc 66704 2013-01-10 18:20:17Z gunter $
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
// 10.04.2014 A.Dotti: Add MT functionality for messenger
// 24.04.2014 A.Ribon: switched on muon-nuclear by default
// 04.04.2020 D. Rivera -- copy based on :
//    geant4.10.03.p03/source/physics_lists/constructors/gamma_lepto_nuclear/src/G4EmExtraPhysics.cc
//
//----------------------------------------------------------------------------
//

//#include "Geant4/G4EmExtraPhysics.hh"
#include "MyG4EmExtraPhysics.hh"

#include "Geant4/G4SystemOfUnits.hh"

#include "Geant4/G4ParticleDefinition.hh"
#include "Geant4/G4ParticleTable.hh"
#include "Geant4/G4Gamma.hh"
#include "Geant4/G4Electron.hh"
#include "Geant4/G4Positron.hh"
#include "Geant4/G4MuonPlus.hh"
#include "Geant4/G4MuonMinus.hh"

#include "Geant4/G4SynchrotronRadiation.hh"

//#include "Geant4/G4BertiniElectroNuclearBuilder.hh"
#include "MyG4BertiniElectroNuclearBuilder.hh"

#include "Geant4/G4MuonNuclearProcess.hh"
#include "Geant4/G4MuonVDNuclearModel.hh"

#include "Geant4/G4GammaConversionToMuons.hh"
#include "Geant4/G4AnnihiToMuPair.hh"
#include "Geant4/G4eeToHadrons.hh"

#include "Geant4/G4PhysicsListHelper.hh"
#include "Geant4/G4BuilderType.hh"
#include "Geant4/G4AutoDelete.hh"
 
// factory
#include "Geant4/G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(MyG4EmExtraPhysics);

G4bool MyG4EmExtraPhysics::gnActivated  = true;
G4bool MyG4EmExtraPhysics::munActivated = true;
G4bool MyG4EmExtraPhysics::synActivated = false;
G4bool MyG4EmExtraPhysics::synActivatedForAll = false;
G4bool MyG4EmExtraPhysics::gmumuActivated = false;
G4bool MyG4EmExtraPhysics::pmumuActivated = false;
G4bool MyG4EmExtraPhysics::phadActivated = false;

G4ThreadLocal MyG4BertiniElectroNuclearBuilder* MyG4EmExtraPhysics::theGNPhysics = nullptr;
G4ThreadLocal G4SynchrotronRadiation* MyG4EmExtraPhysics::theSynchRad = nullptr;
G4ThreadLocal G4GammaConversionToMuons* MyG4EmExtraPhysics::theGammaToMuMu = nullptr;
G4ThreadLocal G4AnnihiToMuPair* MyG4EmExtraPhysics::thePosiToMuMu = nullptr;
G4ThreadLocal G4eeToHadrons* MyG4EmExtraPhysics::thePosiToHadrons = nullptr;

MyG4EmExtraPhysics::MyG4EmExtraPhysics(G4int ver): 
  G4VPhysicsConstructor("G4GammaLeptoNuclearPhys"),
  verbose(ver)
{
  theMessenger = new MyG4EmMessenger(this);
  SetPhysicsType(bEmExtra);
  if(verbose > 1) G4cout << "### MyG4EmExtraPhysics" << G4endl;
}

MyG4EmExtraPhysics::MyG4EmExtraPhysics(const G4String&)
  : MyG4EmExtraPhysics(1)
{}

MyG4EmExtraPhysics::~MyG4EmExtraPhysics()
{
  delete theMessenger;
  theMessenger = nullptr;
}

void MyG4EmExtraPhysics::Synch(G4bool val)
{
  synActivated = val;
}

void MyG4EmExtraPhysics::SynchAll(G4bool val)
{
  synActivatedForAll = val;
  if(synActivatedForAll) { synActivated = true; }
}

void MyG4EmExtraPhysics::GammaNuclear(G4bool val)
{
  gnActivated = val;
}

void MyG4EmExtraPhysics::MuonNuclear(G4bool val)
{
  munActivated = val;
}

void MyG4EmExtraPhysics::GammaToMuMu(G4bool val)
{
  gmumuActivated = val;
}

void MyG4EmExtraPhysics::PositronToMuMu(G4bool val)
{
  pmumuActivated = val;
}

void MyG4EmExtraPhysics::PositronToHadrons(G4bool val)
{
  phadActivated = val;
}

void MyG4EmExtraPhysics::ConstructParticle()
{
  G4Gamma::Gamma();
  G4Electron::Electron();
  G4Positron::Positron();
  G4MuonPlus::MuonPlus();
  G4MuonMinus::MuonMinus();
}

void MyG4EmExtraPhysics::ConstructProcess()
{
  G4ParticleDefinition* gamma = G4Gamma::Gamma();
  G4ParticleDefinition* electron = G4Electron::Electron();
  G4ParticleDefinition* positron = G4Positron::Positron();
  G4ParticleDefinition* muonplus = G4MuonPlus::MuonPlus();
  G4ParticleDefinition* muonminus = G4MuonMinus::MuonMinus();

  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  if(gnActivated) {
    theGNPhysics = new MyG4BertiniElectroNuclearBuilder();
    theGNPhysics->Build();
    //G4AutoDelete::Register(theGNPhysics);
  }
  if(munActivated) {
    G4MuonNuclearProcess* muNucProcess = new G4MuonNuclearProcess();
    G4MuonVDNuclearModel* muNucModel = new G4MuonVDNuclearModel();
    muNucProcess->RegisterMe(muNucModel);
    ph->RegisterProcess( muNucProcess, muonplus);
    ph->RegisterProcess( muNucProcess, muonminus);
  }
  if(gmumuActivated) {
    theGammaToMuMu = new G4GammaConversionToMuons();
    ph->RegisterProcess(theGammaToMuMu, gamma);
  }  
  if(pmumuActivated) {
    thePosiToMuMu = new G4AnnihiToMuPair();
    ph->RegisterProcess(thePosiToMuMu, positron);
  }  
  if(phadActivated) {
    thePosiToHadrons = new G4eeToHadrons();
    ph->RegisterProcess(thePosiToHadrons, positron);
  }  
  if(synActivated) {
    theSynchRad = new G4SynchrotronRadiation();
    ph->RegisterProcess( theSynchRad, electron);
    ph->RegisterProcess( theSynchRad, positron);
    //G4AutoDelete::Register(theSynchRad);
    if(synActivatedForAll) {
      auto myParticleIterator=GetParticleIterator();
      myParticleIterator->reset();
      G4ParticleDefinition* particle = nullptr;

      while( (*myParticleIterator)() ) {
        particle = myParticleIterator->value();
        if( particle->GetPDGStable() && particle->GetPDGCharge() != 0.0) {
          if(verbose > 1) {
            G4cout << "### G4SynchrotronRadiation for "
                   << particle->GetParticleName() << G4endl;
          }
          ph->RegisterProcess( theSynchRad, particle);
        }
      }
    }
  }
}

