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
// 20191209  D. Rivera -- Created to store outgoing particle info at each reaction

#ifndef MYG4ANALYSER2_HH
#define MYG4ANALYSER2_HH

#define WITH_NUCLEI

#include "Geant4/G4CollisionOutput.hh"
#include "Geant4/G4InuclParticle.hh"
#include "Geant4/G4InuclElementaryParticle.hh"
#include "Geant4/G4InuclNuclei.hh"
#include "Geant4/G4NuclWatcher.hh"
//#include "G4ExitonConfiguration.hh"
#include "Geant4/G4String.hh"

#include <vector>
#include <map>

class MyG4Analyser2 {

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

public:

  MyG4Analyser2();
  void setOutputFile(G4String outFileName);
  void printLineBreak();
  void setInelCsec(G4double csec, G4bool withn);
  void setWatchers(const std::vector<G4NuclWatcher>& watchers);
  void try_watchers(G4int a, G4int z, G4bool if_nucl);
  void analyse(const G4CollisionOutput& output, const G4InuclParticle& bullet);
  void printResults();
  void printResultsSimple();
  void handleWatcherStatistics();
  void printBulletNtuple(std::ostream& outfile, const G4InuclParticle &particle, G4int Id, G4int nDaughters) const;
  void printParticleNtuple(std::ostream& outfile, const G4InuclElementaryParticle &particle, G4int Id, G4int nDaughters) const;
  void printResultsNtuple();

private: 

  G4int verboseLevel;
  std::ofstream outFile;
  G4bool outputSet;
  G4int pos;

  G4int bulletPDG;
  G4String bulletType;
  G4double bulletKineticEnergy;
  G4double eventNumber;
  G4double Multiplicity;
  G4double ProtonNumber;
  G4double NeutronNumber;
  G4double PionNumber;
  G4double PhotonNumber;
  G4double NucleonKinEnergy;
  G4double ProtonKinEnergy;
  G4double NeutronKinEnergy;
  G4double PionKinEnergy;
  G4double PhotonKinEnergy;
  G4double ExitationEnergy;
  G4double OutgoingNuclei;

  G4double fissy_prob;

  G4double PionPl;
  G4double PionMin;
  G4double Pion0;
  G4double A;
  G4double Z;

  std::vector<G4NuclWatcher> ana_watchers;
  G4double inel_csec;
  G4bool withNuclei;

  std::map<G4int, std::pair<std::string,G4int> > modelCounterMap;
};        

#endif // MYG4ANALYSER2_HH
