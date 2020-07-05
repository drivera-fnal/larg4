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
// 20121017  M. Kelsey -- Add local cache of Bertini pointer
// 20200704  D. Rivera -- Copy based on :
//    geant4.10.03.p03/source/processes/hadronic/stopping/include/G4HadronicAbsorptionBertini.hh
//
//----------------------------------------------------------------------------

#ifndef MyG4HadronicAbsorptionBertini_h
#define MyG4HadronicAbsorptionBertini_h 1

// Class Description:
//
// Intermediate base (or concrete) class for hadronic absorption at rest. 
// Physics lists should reference the concrete subclasses for pi-, K-, Sigma-

#include "Geant4/globals.hh"
#include "Geant4/G4HadronStoppingProcess.hh"
#include <iosfwd>

class G4VParticleChange;
class G4ParticleDefinition;
//class G4CascadeInterface;
class MyG4CascadeInterface;
class G4Track;


class MyG4HadronicAbsorptionBertini : public G4HadronStoppingProcess { 
public:
  // May instantiate this class directly for all three pi-, K-, Sigma-
  MyG4HadronicAbsorptionBertini(G4ParticleDefinition* pdef=0);
  virtual ~MyG4HadronicAbsorptionBertini() {;}
  
  G4bool IsApplicable(const G4ParticleDefinition&);

  void ProcessDescription(std::ostream& outFile) const;

private:
  // hide assignment operator as private 
  MyG4HadronicAbsorptionBertini& operator=(const MyG4HadronicAbsorptionBertini&);
  MyG4HadronicAbsorptionBertini(const MyG4HadronicAbsorptionBertini&);
  
private:
  G4ParticleDefinition* pdefApplicable;
  MyG4CascadeInterface* theCascade;
};

#endif

