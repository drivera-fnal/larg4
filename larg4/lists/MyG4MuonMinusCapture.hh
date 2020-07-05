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
// $Id: G4MuonMinusCapture.hh 66367 2012-12-18 09:18:08Z gcosmo $
//
//---------------------------------------------------------------------
//
// GEANT4 Class header file
//
// File name:     G4MuonMinusCapture
//
// Author V.Ivanchenko 25 April 2012 on base of G4MuonMinusCaptureAtRest
//
//
// Class Description:
//
// Stopping of mu-
//
// Modifications: 
//   20121003 K. Genser - Changed the constructor argument type
//   20200704 D. Rivera -- Copy based on :
//    geant4.10.03.p03/source/processes/hadronic/stopping/include/G4MuonMinusCapture.hh
//
//------------------------------------------------------------------------

#ifndef MyG4MuonMinusCapture_h
#define MyG4MuonMinusCapture_h 1
 
#include "Geant4/globals.hh"
#include "Geant4/G4HadronStoppingProcess.hh"
#include "Geant4/G4ParticleDefinition.hh"

class G4HadronicInteraction;

class MyG4MuonMinusCapture : public G4HadronStoppingProcess 
{ 
public:
 
  explicit MyG4MuonMinusCapture(G4HadronicInteraction* hiptr=0);

  ~MyG4MuonMinusCapture();

  G4bool IsApplicable(const G4ParticleDefinition&);

  void ProcessDescription(std::ostream& outFile) const;

private:

  // hide assignment operator as private 
  MyG4MuonMinusCapture& operator=(const MyG4MuonMinusCapture &right);
  MyG4MuonMinusCapture(const MyG4MuonMinusCapture& );

};

#endif
 





