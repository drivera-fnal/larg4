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
// GEANT4 tag $Name: $
//
//---------------------------------------------------------------------------
//
// ClassName:   MyG4HyperonFTFPBuilder
//
// Author: 2012 G.Folger
//
// Modified:
// 30.03.2020 D. Rivera copy based on:
//    geant4.10.03.p03/source/physics_lists/builders/include/G4HyperonFTFPBuilder.hh
//
//----------------------------------------------------------------------------
//
#ifndef MyG4HyperonFTFPBuilder_h
#define MyG4HyperonFTFPBuilder_h 1

#include "Geant4/globals.hh"

#include "Geant4/G4LambdaInelasticProcess.hh"
#include "Geant4/G4AntiLambdaInelasticProcess.hh"
#include "Geant4/G4SigmaPlusInelasticProcess.hh"
#include "Geant4/G4SigmaMinusInelasticProcess.hh"
#include "Geant4/G4AntiSigmaPlusInelasticProcess.hh"
#include "Geant4/G4AntiSigmaMinusInelasticProcess.hh"
#include "Geant4/G4XiZeroInelasticProcess.hh"
#include "Geant4/G4XiMinusInelasticProcess.hh"
#include "Geant4/G4AntiXiZeroInelasticProcess.hh"
#include "Geant4/G4AntiXiMinusInelasticProcess.hh"
#include "Geant4/G4OmegaMinusInelasticProcess.hh"
#include "Geant4/G4AntiOmegaMinusInelasticProcess.hh"

#include "Geant4/G4TheoFSGenerator.hh"
#include "Geant4/G4GeneratorPrecompoundInterface.hh"
#include "Geant4/G4FTFModel.hh"
#include "Geant4/G4LundStringFragmentation.hh"
#include "Geant4/G4ExcitedStringDecay.hh"
//#include "Geant4/G4CascadeInterface.hh"
#include "MyG4CascadeInterface.hh"

#include "Geant4/G4ChipsHyperonInelasticXS.hh"
class MyG4HyperonFTFPBuilder
{
  public:
    MyG4HyperonFTFPBuilder();
    virtual ~MyG4HyperonFTFPBuilder();

  public:
    void Build();

  private:

    G4TheoFSGenerator * HyperonFTFP;
    G4TheoFSGenerator * AntiHyperonFTFP;
    G4GeneratorPrecompoundInterface * theCascade;
    G4FTFModel * theStringModel;
    G4ExcitedStringDecay * theStringDecay;
    G4LundStringFragmentation * theLund;
    //<--G4CascadeInterface * theBertini;
    MyG4CascadeInterface * theBertini;
    G4LambdaInelasticProcess*  theLambdaInelastic;
    G4AntiLambdaInelasticProcess*  theAntiLambdaInelastic;
    G4SigmaMinusInelasticProcess*  theSigmaMinusInelastic;
    G4AntiSigmaMinusInelasticProcess*  theAntiSigmaMinusInelastic;
    G4SigmaPlusInelasticProcess*  theSigmaPlusInelastic;
    G4AntiSigmaPlusInelasticProcess*  theAntiSigmaPlusInelastic;
    G4XiZeroInelasticProcess*  theXiZeroInelastic;
    G4AntiXiZeroInelasticProcess*  theAntiXiZeroInelastic;
    G4XiMinusInelasticProcess*  theXiMinusInelastic;
    G4AntiXiMinusInelasticProcess*  theAntiXiMinusInelastic;
    G4OmegaMinusInelasticProcess*  theOmegaMinusInelastic;
    G4AntiOmegaMinusInelasticProcess*  theAntiOmegaMinusInelastic;

    //  G4QHadronInelasticDataSet * theCHIPSInelastic;
    G4VCrossSectionDataSet* theCHIPSInelastic;
    G4bool wasActivated;
};
#endif
