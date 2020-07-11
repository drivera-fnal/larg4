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
// 04.04.2020 D. Rivera -- copy based on :
//    geant4.10.03.p03/source/physics_lists/constructors/gamma_lepto_nuclear/include/G4BertiniElectronNuclearBuilder.hh
//
//----------------------------------------------------------------------------
//

#ifndef MyG4BertiniElectroNuclearBuilder_h
#define MyG4BertiniElectroNuclearBuilder_h 1

#include "Geant4/globals.hh"
#include "Geant4/G4ios.hh"

#include "Geant4/G4TheoFSGenerator.hh"
#include "Geant4/G4GeneratorPrecompoundInterface.hh"
#include "Geant4/G4QGSModel.hh"
#include "Geant4/G4GammaParticipants.hh"
#include "Geant4/G4QGSMFragmentation.hh"
#include "Geant4/G4ExcitedStringDecay.hh"

//#include "Geant4/G4CascadeInterface.hh"
#include "MyG4CascadeInterface.hh"
#include "Geant4/G4ElectroVDNuclearModel.hh"
//#include "G4ElectroNuclearReaction.hh"
#include "Geant4/G4PhotoNuclearProcess.hh"
#include "Geant4/G4ElectronNuclearProcess.hh"
#include "Geant4/G4PositronNuclearProcess.hh"

//A. Dotti (June2013): No need to change this class for MT
// Since each thread owns its own instance (created by G4EmExtraPhysics)

class MyG4BertiniElectroNuclearBuilder 
{
  public: 
    MyG4BertiniElectroNuclearBuilder();
    virtual ~MyG4BertiniElectroNuclearBuilder();

  public: 
    virtual void Build();

  protected:
    G4PhotoNuclearProcess * thePhotoNuclearProcess;
    G4ElectronNuclearProcess * theElectronNuclearProcess;
    G4PositronNuclearProcess * thePositronNuclearProcess;
    G4ElectroVDNuclearModel * theElectroReaction;
    //G4CascadeInterface * theGammaReaction;  
    MyG4CascadeInterface * theGammaReaction;  
    
    G4TheoFSGenerator * theModel;
    G4GeneratorPrecompoundInterface * theCascade;
    G4QGSModel< G4GammaParticipants > * theStringModel;
    G4QGSMFragmentation * theFragmentation;
    G4ExcitedStringDecay * theStringDecay;
    G4bool wasActivated;
};

// 2002 by J.P. Wellisch

#endif

