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
// $Id: G4BertiniPionBuilder.hh 66892 2013-01-17 10:57:59Z gunter $
//
//---------------------------------------------------------------------------
//
// ClassName:   MyG4BertiniPionBuilder
//
// Author: 2010 G.Folger
//  devired from G4BertiniPiKBuilder
//
// Modified:
// 30.03.2009 V.Ivanchenko create cross section by new
// 09.12.2019 D. Rivera copy based on:
//    geant4.10.03.p03/source/physics_lists/builders/include/G4BertiniPionBuilder.hh
//
//----------------------------------------------------------------------------
//
#ifndef MyG4BertiniPionBuilder_h
#define MyG4BertiniPionBuilder_h 1

#include "Geant4/globals.hh"

#include "Geant4/G4HadronElasticProcess.hh"
#include "Geant4/G4HadronFissionProcess.hh"
#include "Geant4/G4HadronCaptureProcess.hh"
#include "Geant4/G4NeutronInelasticProcess.hh"
#include "Geant4/G4VPionBuilder.hh"

#include "Geant4/G4PiNuclearCrossSection.hh"
//#include "Geant4/G4CascadeInterface.hh"

#include "MyG4CascadeInterface.hh"

class MyG4BertiniPionBuilder : public G4VPionBuilder
{
  public: 
    MyG4BertiniPionBuilder();
    virtual ~MyG4BertiniPionBuilder();

  public: 
    virtual void Build(G4HadronElasticProcess * aP);
    virtual void Build(G4PionPlusInelasticProcess * aP);
    virtual void Build(G4PionMinusInelasticProcess * aP);
    
    void SetMinEnergy(G4double aM) {theMin = aM;}
    void SetMaxEnergy(G4double aM) {theMax = aM;}

  private:
    G4PiNuclearCrossSection* thePiData;
    MyG4CascadeInterface * theModel;    
    G4double theMin;
    G4double theMax;

};
#endif

