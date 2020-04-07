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
// $Id: G4InuclEvaporation.hh 66241 2012-12-13 18:34:42Z gunter $
// Defines an interface to evaporation models of Bertini cascase (BERT)
// based on INUCL code.
//
// 20100405  M. Kelsey -- Pass const-ref std::vector<>
// 20100517  M. Kelsey -- MakeG4EvaporationInuclCollider a data member.
// 20100520  M. Kelsey -- Clean up interface
// 20200329  D. Rivera -- Copy based on:
//    geant4.10.03.p03/source/processes/hadronic/models/cascade/cascade/include/G4InuclEvaporation.hh

#ifndef MYG4INUCLEVAPORATION_h
#define MYG4INUCLEVAPORATION_h 1

#include "Geant4/globals.hh"
#include "Geant4/G4VEvaporation.hh"
#include "Geant4/G4Fragment.hh"

class G4EvaporationInuclCollider;


class MyG4InuclEvaporation : public G4VEvaporation {
public:
  MyG4InuclEvaporation();
  ~MyG4InuclEvaporation();

private:
  MyG4InuclEvaporation(const MyG4InuclEvaporation &right);

  const MyG4InuclEvaporation & operator=(const MyG4InuclEvaporation &right);
  G4bool operator==(const MyG4InuclEvaporation &right) const;
  G4bool operator!=(const MyG4InuclEvaporation &right) const;

public:
  G4FragmentVector * BreakItUp(const G4Fragment &theNucleus);
      
  void setVerboseLevel( const G4int verbose );

private:
  G4int verboseLevel;
  G4EvaporationInuclCollider* evaporator;
};

#endif
