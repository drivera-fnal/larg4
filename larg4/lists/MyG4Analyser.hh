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
// $Id: MyG4Analyser.hh 66241 2012-12-13 18:34:42Z gunter $
//
// 20101010  M. Kelsey -- Migrate to integer A and Z
// 20191209  D. Rivera -- Copy based on:
//    geant4.10.03.p03/source/processes/hadronic/models/cascade/cascade/include/G4Analyser.hh

#ifndef MYG4ANALYSER_HH
#define MYG4ANALYSER_HH

#define WITH_NUCLEI

#include "Geant4/G4CollisionOutput.hh"
#include "Geant4/G4InuclParticle.hh"
#include "Geant4/G4InuclElementaryParticle.hh"
#include "Geant4/G4InuclNuclei.hh"
#include "Geant4/G4NuclWatcher.hh"
//#include "G4ExitonConfiguration.hh"

#include <vector>
#include <map>

const std::map<G4int,std::string> modelMap = {
  {0,"default"},
  {1,"bullet"},
  {2,"target"},
  {3,"G4ElementaryParticleCollider"},
  {4,"G4IntraNucleiCascader"},
  {5,"G4NonEquilibriumEvaporator"},
  {6,"G4EquilibriumEvaporator"},
  {7,"G4Fissioner"},
  {8,"G4BigBanger"},
  {9,"G4PreCompound"},
  {10,"G4CascadeCoalescence"}
};

class MyG4Analyser {

public:

  MyG4Analyser();
  void setInelCsec(G4double csec, G4bool withn);
  void setWatchers(const std::vector<G4NuclWatcher>& watchers);
  void try_watchers(G4int a, G4int z, G4bool if_nucl);
  void analyse(const G4CollisionOutput& output);
  void printResults();
  void printResultsSimple();
  void handleWatcherStatistics();
  void printResultsNtuple();

private: 

  G4int verboseLevel;
  G4double eventNumber;
  G4double averageMultiplicity;
  G4double averageProtonNumber;
  G4double averageNeutronNumber;
  G4double averagePionNumber;
  G4double averagePhotonNumber;
  G4double averageNucleonKinEnergy;
  G4double averageProtonKinEnergy;
  G4double averageNeutronKinEnergy;
  G4double averagePionKinEnergy;
  G4double averagePhotonKinEnergy;
  G4double averageExitationEnergy;
  G4double averageOutgoingNuclei;
  G4double fissy_prob;
  G4double averagePionPl;
  G4double averagePionMin;
  G4double averagePion0;
  G4double averageA;
  G4double averageZ;
  std::vector<G4NuclWatcher> ana_watchers;
  G4double inel_csec;
  G4bool withNuclei;

  std::map<G4int, std::pair<std::string,G4int> > modelCounterMap;
};        

#endif // MYG4ANALYSER_HH
