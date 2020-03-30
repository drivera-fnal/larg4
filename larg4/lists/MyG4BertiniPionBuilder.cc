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
// $Id: MyG4BertiniPionBuilder.cc 83699 2014-09-10 07:18:25Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   MyG4BertiniPionBuilder
//
// Author: 2010 G.Folger
//  devired from G4BertiniPiKBuilder
//
// Modified:
// 02.04.2009 V.Ivanchenko remove add cross section, string builderis reponsible 
//
//----------------------------------------------------------------------------
//
#include "MyG4BertiniPionBuilder.hh"
#include "Geant4/G4SystemOfUnits.hh"
#include "Geant4/G4ParticleDefinition.hh"
#include "Geant4/G4ParticleTable.hh"
#include "Geant4/G4ProcessManager.hh"
#include "Geant4/G4CrossSectionDataSetRegistry.hh"

MyG4BertiniPionBuilder::
MyG4BertiniPionBuilder() 
 {
   thePiData = (G4PiNuclearCrossSection*)G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4PiNuclearCrossSection::Default_Name());
   theMin = 0*GeV;
   theMax = 9.9*GeV;
   theModel = new MyG4CascadeInterface;
   theModel->SetMinEnergy(theMin);
   theModel->SetMaxEnergy(theMax); 
 }

MyG4BertiniPionBuilder::~MyG4BertiniPionBuilder() 
{
}

void MyG4BertiniPionBuilder::
Build(G4PionPlusInelasticProcess * aP)
 {
   aP->RegisterMe(theModel);
   theModel->SetMinEnergy(theMin);
   theModel->SetMaxEnergy(theMax);
 }

void MyG4BertiniPionBuilder::
Build(G4PionMinusInelasticProcess * aP)
 {
   aP->RegisterMe(theModel);
   theModel->SetMinEnergy(theMin);
   theModel->SetMaxEnergy(theMax);
 }

void MyG4BertiniPionBuilder::
Build(G4HadronElasticProcess * ) {}

         