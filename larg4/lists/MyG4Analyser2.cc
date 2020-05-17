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
// $Id: MyG4Analyser.cc 66241 2012-12-13 18:34:42Z gunter $
//
// 20100726  M. Kelsey -- Use references for fetched lists
// 20101010  M. Kelsey -- Migrate to integer A and Z
// 20101019  M. Kelsey -- CoVerity report, unitialized constructor
// 20191209  D. Rivera -- Copy based on:
//    geant4.10.03.p03/source/processes/hadronic/models/cascade/cascade/src/G4Analyser.cc
// 20200129  D. Rivera -- Fixed pi-, pi0, and pi+ counters
// 20200324  D. Rivera -- Added photon counters, and photon averageKE info
// 20200325  D. Rivera -- Added Model counter

#include "MyG4Analyser2.hh"
#include <cmath>
#include <iomanip>

MyG4Analyser2::MyG4Analyser2()
  : verboseLevel(0), outputSet(false), bulletPDG(0), bulletType("unknown"),
    bulletKineticEnergy(0.0), eventNumber(0.0), Multiplicity(0.0),
    ProtonNumber(0.0), NeutronNumber(0.0),
    PionNumber(0.0), PhotonNumber(0.0),
    NucleonKinEnergy(0.0), ProtonKinEnergy(0.0),
    NeutronKinEnergy(0.0), PionKinEnergy(0.0), PhotonKinEnergy(0.0),
    ExitationEnergy(0.0), OutgoingNuclei(0.0), fissy_prob(0.0), PionPl(0.0),
    PionMin(0.0), Pion0(0.0), A(0.0), Z(0.0),
    inel_csec(0.0), withNuclei(true) {
    // -- enable analyzer with Nuclei enabled
  if (verboseLevel > 3) {
    G4cout << " >>> MyG4Analyser2::MyG4Analyser2" << G4endl;
  }

  for (auto const &[i,model] : modelMap){
    modelCounterMap.emplace(i,std::make_pair(model, 0));
    G4cout << "Entry " << i << " : " << model << G4endl;
  }
}

void MyG4Analyser2::setOutputFile(G4String outFileName) {
  outFile.open(outFileName, std::ios::out | std::ios::app );
  outputSet = true;

  G4cerr << "INSIDE setOutputFile! " << G4endl;
  pos = outFile.tellp(); // -- position in the stream
  G4cerr << "File is open: " << outFile.is_open() << " Position : " << pos << G4endl;

  // -- Header
  if (pos==0)
  {
    G4cerr << "INSIDE setOutputFile if loop!" << G4endl;
    outFile << "#       PDG"
            << " Mdl"
            << "   Id"
            << "  MId"
            << "  Gen"
            << " Daug"
            << " Zone"
            << "   In"
            << "             KE"
            << "              E"
            << "             Px"
            << "             Py"
            << "             Pz"
            << "             Vx"
            << "             Vy"
            << "             Vz"
            << std::endl;
    pos = outFile.tellp();
  }
  G4cerr << "File is open: " << outFile.is_open() << " Position : " << pos << G4endl;
}

void MyG4Analyser2::
printLineBreak() {

  std::streambuf *psbuf, *backup;

  backup = std::cout.rdbuf(); // -- back up cout's streambuf

  psbuf = outFile.rdbuf();    // -- get file's streambuf
  std::cout.rdbuf(psbuf);     // -- assign streambuf to cout

  std::cout << std::setfill('#') << std::setw(165) << '#' << std::endl; // -- written to the file

  std::cout.rdbuf(backup);    // -- restore cout's original streambuf

}

void MyG4Analyser2::setInelCsec(G4double csec, G4bool withn) {

  if (verboseLevel > 3) {
    G4cout << " >>> MyG4Analyser2::setInelCsec" << G4endl;
  }

  inel_csec = csec; // mb
  withNuclei = withn;

  if (verboseLevel > 3) {
    G4cout << " total inelastic " << inel_csec << G4endl;
  }
}

void MyG4Analyser2::setWatchers(const std::vector<G4NuclWatcher>& watchers) {

  if (verboseLevel > 3) {
    G4cout << " >>> MyG4Analyser2::setWatchers" << G4endl;
  }

  ana_watchers = watchers;
  if (verboseLevel > 3) {
    G4cout << " watchers set " << watchers.size() << G4endl;
  }
}

void MyG4Analyser2::try_watchers(G4int a, G4int z, G4bool if_nucl) {

  if (verboseLevel > 3) {
    G4cout << " >>> MyG4Analyser2::try_watchers" << G4endl;
  }

  for (G4int iw = 0; iw < G4int(ana_watchers.size()); iw++) { 

    if (if_nucl) {

      if (ana_watchers[iw].look_forNuclei()) ana_watchers[iw].watch(a, z); 

    } else {

      if (!ana_watchers[iw].look_forNuclei()) ana_watchers[iw].watch(a, z); 
    }
  }
}


void MyG4Analyser2::analyse(const G4CollisionOutput& output, const G4InuclParticle& bullet) {

  if (verboseLevel > 3) {
    G4cout << " >>> MyG4Analyser2::analyse" << G4endl;
  }

  bulletKineticEnergy = bullet.getKineticEnergy();
  bulletType = bullet.getDefinition()->GetParticleName();
  bulletPDG  = bullet.getDefinition()->GetPDGEncoding();

  const std::vector<G4InuclElementaryParticle>& particles =
    output.getOutgoingParticles();
  Multiplicity += particles.size();

  // -- Add entry for primary bullet particle  (ostream, partile, mId, nDaughters
  if (outputSet) printBulletNtuple(outFile, bullet, 0, particles.size());


  if (withNuclei) { // -- include nuclei analysis
    const std::vector<G4InuclNuclei>& nuclei = output.getOutgoingNuclei();

    //    if (nuclei.size() >= 0) {
    if (nuclei.size() > 0) { // -- if there are no nuclei.. this gets skipped?

      // -- process particles first
      for (G4int i = 0; i < G4int(particles.size()); i++) {
        G4int ap = 0;
        G4int zp = 0;
        G4InuclParticle::Model model = particles[i].getModel();
        G4int modelId = (G4int)model;
        modelCounterMap[modelId].second +=1;
        if (outputSet) printParticleNtuple(outFile, particles[i], (i+1));

        if (particles[i].nucleon()) { // -- is a nucleon
          NucleonKinEnergy += particles[i].getKineticEnergy();

          if (particles[i].type() == 1) {
            zp = 1;
            ap = 1;
            ProtonNumber += 1.0;
            ProtonKinEnergy += particles[i].getKineticEnergy();

          } else {
            ap = 1;
            zp = 0;
            NeutronNumber += 1.0;
            NeutronKinEnergy += particles[i].getKineticEnergy();
          };  

        } else if (particles[i].pion()) {
          PionKinEnergy += particles[i].getKineticEnergy();
          PionNumber += 1.0;
          ap = 0;

          if (particles[i].type() == 3) {
            zp = 1;
            PionPl += 1.0;

          } else if (particles[i].type() == 5) {  
            zp = -1;
            PionMin += 1.0;

          } else if (particles[i].type() == 7) { 
            zp = 0;
            Pion0 += 1.0;
          };
        } else if (particles[i].isPhoton()) {
          PhotonKinEnergy += particles[i].getKineticEnergy();
          PhotonNumber += 1.0;
          ap = 0;
          zp = 0;
        }; // -- 
        try_watchers(ap, zp, false);
      };

      G4int nbig = 0;
      OutgoingNuclei += nuclei.size();

      for (G4int in = 0; in < G4int(nuclei.size()); in++)
      {
        ExitationEnergy += nuclei[in].getExitationEnergy();

        G4int a = nuclei[in].getA();
        G4int z = nuclei[in].getZ();

        if (in == 0) {  // -- why only the first entry in the nuclei vector?
          A = a;
          Z = z;
        };

        if (a > 10) nbig++;
        try_watchers(a, z, true);

        if (outputSet) printNucleiNtuple(outFile, nuclei[in], (in+1));
      };

      if (nbig > 1) fissy_prob += 1.0;
      eventNumber += 1.0;
      /*
      const std::vector<G4InuclElementaryParticle>& particles =
        output.getOutgoingParticles();
      Multiplicity += particles.size();

      // -- Add entry for primary bullet particle  (ostream, partile, mId, nDaughters
      if (outputSet) printBulletNtuple(outFile, bullet, 0, particles.size());
      */
      // -- loop over all outgoing particles

    } else {
      G4cout << " >> MyG4Analyser2::analyse "
             << "withNuclei = true, but vector of outgoing nuclei is empty..." << G4endl;
    };

  } else { // -- don't do nuclei analysis
    eventNumber += 1.0;
    /*
    const std::vector<G4InuclElementaryParticle>& particles =
      output.getOutgoingParticles();
    Multiplicity += particles.size();
    */
    for (G4int i = 0; i < G4int(particles.size()); i++) {
      G4InuclParticle::Model model = particles[i].getModel();
      G4int modelId = (G4int)model;
      modelCounterMap[modelId].second +=1;

      if (particles[i].nucleon()) {
        NucleonKinEnergy += particles[i].getKineticEnergy();

        if (particles[i].type() == 1) {
          ProtonNumber += 1.0;
          ProtonKinEnergy += particles[i].getKineticEnergy();

        } else {
          NeutronNumber += 1.0;
          NeutronKinEnergy += particles[i].getKineticEnergy();
        }

      } else if (particles[i].pion()) {
        PionNumber += 1.0;
        PionKinEnergy += particles[i].getKineticEnergy();

        if (particles[i].type() == 3) {
          PionPl += 1.0;

        } else if (particles[i].type() == 5) {  
          PionMin += 1.0;

        } else if (particles[i].type() == 7) { 
          Pion0 += 1.0;
        }
      } else if (particles[i].isPhoton()) {
        PhotonNumber +=1.0;
        PhotonKinEnergy += particles[i].getKineticEnergy();
      }
    }
  }
}

void MyG4Analyser2::printResultsSimple() {

  if (verboseLevel > 3) {
    G4cout << " >>> MyG4Analyser2::printResultsSimple" << G4endl;
  }

  G4cout << " Number of events " << eventNumber << G4endl
         << "  multiplicity " << Multiplicity << G4endl
         << "  proton number " << ProtonNumber << G4endl
         << "  neutron number " << NeutronNumber << G4endl
         << "  total nucleon Ekin " << NucleonKinEnergy << G4endl
         << "  total proton Ekin " << ProtonKinEnergy << G4endl
         << "  total neutron Ekin " << NeutronKinEnergy << G4endl
         << "  pion number " << PionNumber << G4endl
         << "  total pion Ekin " << PionKinEnergy << G4endl
         << "  photon number " << PhotonNumber << G4endl
         << "  total photon Ekin " << PhotonKinEnergy << G4endl;

  if (withNuclei) {
    G4cout                 
      << "  total Exitation Energy " << ExitationEnergy << G4endl
      << "  num of fragments " << OutgoingNuclei << G4endl;

    G4cout 
      << " fission prob. " << fissy_prob << " c.sec "
      << inel_csec * fissy_prob << G4endl;
  }
}

void MyG4Analyser2::printResults() {

  if (verboseLevel > 3) {
    G4cout << " >>> MyG4Analyser2::printResults" << G4endl;
  }

  G4cout << " Event Number " << eventNumber << G4endl
         << " incident particle type " << bulletType << " [" << bulletPDG << "]" <<  G4endl
         << " incident particle kinetic energy " << bulletKineticEnergy << G4endl
         << "  multiplicity " << Multiplicity << G4endl
         << "  proton number " << ProtonNumber << G4endl
         << "  neutron number " << NeutronNumber << G4endl
         << "  total nucleon Ekin " << NucleonKinEnergy << G4endl
         << "  total proton Ekin " << ProtonKinEnergy << G4endl
         << "  total neutron Ekin " << NeutronKinEnergy << G4endl
         << "  pion number " << PionNumber << G4endl
         << "  total pion Ekin " << PionKinEnergy << G4endl
         << "  pi+ number " << PionPl << G4endl
         << "  pi- number " << PionMin << G4endl
         << "  pi0 number " << Pion0 << G4endl
         << "  photon number " << PhotonNumber << G4endl
         << "  total photon Ekin " << PhotonKinEnergy << G4endl;
                   
  if (withNuclei) {
    G4cout
      << "  A[0] " << A << G4endl                 
      << "  Z[0] " << Z << G4endl                 
      << "  total excitation Energy " << ExitationEnergy << G4endl
      << "  num of fragments " << OutgoingNuclei << G4endl;

    G4cout
      << " fission prob. " << fissy_prob << " c.sec "
      << inel_csec * fissy_prob << G4endl;

    handleWatcherStatistics();
  }
  G4cout << "//////////-Model Statistics-//////////" << G4endl;
  for(auto const &[m,info] : modelCounterMap){
    G4cout << "Model " << info.first << " : " << info.second << G4endl;
  }
}

void MyG4Analyser2::handleWatcherStatistics() {

  if (verboseLevel > 3) {
    G4cout << " >>> MyG4Analyser2::handleWatcherStatistics" << G4endl;
  }

  // const G4double small = 1.0e-10;

  if (verboseLevel > 3) {
    G4cout << " >>>Isotope analysis:" << G4endl;
  };

  G4double fgr = 0.0;
  G4double averat = 0.0;
  G4double ave_err = 0.0;
  G4double gl_chsq = 0.0;
  G4double tot_exper = 0.0;
  G4double tot_exper_err = 0.0;
  G4double tot_inucl = 0.0;
  G4double tot_inucl_err = 0.0;
  G4double checked = 0.0;

  for (G4int iw = 0; iw < G4int(ana_watchers.size()); iw++) {
    ana_watchers[iw].setInuclCs(inel_csec, G4int(eventNumber));
    ana_watchers[iw].print();

    if (ana_watchers[iw].to_check()) {
      std::pair<G4double, G4double> rat_err = ana_watchers[iw].getAverageRatio();
      averat += rat_err.first;
      ave_err += rat_err.second;
      gl_chsq += ana_watchers[iw].getChsq();   
      std::pair<G4double, G4double> cs_err = ana_watchers[iw].getExpCs();
      tot_exper += cs_err.first;
      tot_exper_err += cs_err.second;
      std::pair<G4double, G4double> inucl_cs_err = ana_watchers[iw].getInuclCs();
      tot_inucl += inucl_cs_err.first;
      tot_inucl_err += inucl_cs_err.second;
      G4double iz_checked = ana_watchers[iw].getNmatched();

      if (iz_checked > 0.0) {
        fgr += ana_watchers[iw].getLhood();
        checked += iz_checked;    
      };
    };
  };

  if (checked > 0.0) {
    gl_chsq = std::sqrt(gl_chsq) / checked;
    averat /= checked;
    ave_err /= checked;
    fgr = std::pow(10.0, std::sqrt(fgr / checked)); 
  };

  if (verboseLevel > 3) {
    G4cout << " total exper c.s. " << tot_exper << " err " << tot_exper_err <<
      " tot inucl c.s. " << tot_inucl << " err " << tot_inucl_err << G4endl;
    G4cout << " checked total " << checked << " lhood " << fgr << G4endl
           << "  ratio " << averat << " err " << ave_err << G4endl
           << " global chsq " << gl_chsq << G4endl;
  }
}

void MyG4Analyser2::printParticleNtuple(std::ostream& outfile, const G4InuclElementaryParticle &particle, G4int Id=-9) const {
  if (verboseLevel > 3) {
    G4cout << " >>> MyG4Analyser2::printParticleNtuple" << G4endl;
  }

  // -- reconstruct the generation based on the model
  int Gen, modl = particle.getModel();
  if (modl < 2) Gen = 0; // -- Default(0), bullet(1), or target(1)
  else if (modl == 3) Gen = 1; // -- EPCollider(3)
  else Gen = 2; // -- INCascader(4) or above, the true generation is lost at this point, but I've
                //    attemped to fix the over-use of EPCollider model (vs. the INCascader) within
                //    MyG4NucleiModel by re-classifying all EPCollider secondaries with generation
                //    > 1 as INCascader particles. This will only affect the Generation value for
                //    G4CascadParticles. Generation of models > 3 default to 2, but are not really
                //    accurate.

  // Create one line of ASCII data.
  outfile <<
    std::setw(11) << particle.getDefinition()->GetPDGEncoding() <<  /*PDG*/
    std::setw(4)  << particle.getModel() <<                         /*Modl*/
    std::setw(5)  << Id <<                                          /*Id*/
    std::setw(5)  << 0 <<                                           /*MId*/
    std::setw(5)  << Gen <<                                         /*Gen*/
    std::setw(5)  << 0 <<                                           /*Daug*/
    std::setw(5)  << 3 <<                                           /*Zone*/
    std::setw(5)  << 0 <<                                           /*In*/
    std::setw(15) << particle.getKineticEnergy() <<                 /*KE*/
    std::setw(15) << (particle.getMomentum())[3] <<                 /*E*/
    std::setw(15) << (particle.getMomentum())[0] <<                 /*Px*/
    std::setw(15) << (particle.getMomentum())[1] <<                 /*Py*/
    std::setw(15) << (particle.getMomentum())[2] <<                 /*Pz*/
    std::setw(15) << -99. <<                                        /*Vx*/
    std::setw(15) << -99. <<                                        /*Vy*/
    std::setw(15) << -99.                                           /*Vz*/
    << std::endl;
}

void MyG4Analyser2::printBulletNtuple(std::ostream& outfile, const G4InuclParticle &particle, G4int Id=-9, G4int nDaughters=-1) const {
  if (verboseLevel > 3) {
    G4cout << " >>> MyG4Analyser2::printBulletNtuple" << G4endl;
  }

  // Create one line of ASCII data.
  // all bullet particles are considered Generation 0
  outfile <<
    std::setw(11) << particle.getDefinition()->GetPDGEncoding() <<  /*PDG*/
    std::setw(4)  << particle.getModel() <<                         /*Modl*/
    std::setw(5)  << Id <<                                          /*Id*/
    std::setw(5)  << -1 <<                                          /*MId*/
    std::setw(5)  << 0 <<                                           /*Gen*/
    std::setw(5)  << nDaughters <<                                  /*Daug*/
    std::setw(5)  << 3 <<                                           /*Zone*/
    std::setw(5)  << 0 <<                                           /*In*/
    std::setw(15) << particle.getKineticEnergy() <<                 /*KE*/
    std::setw(15) << (particle.getMomentum())[3] <<                 /*E*/
    std::setw(15) << (particle.getMomentum())[0] <<                 /*Px*/
    std::setw(15) << (particle.getMomentum())[1] <<                 /*Py*/
    std::setw(15) << (particle.getMomentum())[2] <<                 /*Pz*/
    std::setw(15) << -99. <<                                        /*Vx*/
    std::setw(15) << -99. <<                                        /*Vy*/
    std::setw(15) << -99.                                           /*Vz*/
    << std::endl;
}

void MyG4Analyser2::printNucleiNtuple(std::ostream& outfile, const G4InuclNuclei &particle, G4int Id=-9) const {
  if (verboseLevel > 3) {
    G4cout << " >>> MyG4Analyser2::printNucleiNtuple" << G4endl;
  }

  // Create one line of ASCII data.
  // all bullet particles are considered Generation 0
  outfile <<
    std::setw(11) << particle.getDefinition()->GetPDGEncoding() <<  /*PDG*/
    std::setw(4)  << particle.getModel() <<                         /*Modl*/
    std::setw(5)  << Id <<                                          /*Id*/
    std::setw(5)  << -1 <<                                          /*MId*/
    std::setw(5)  << 0 <<                                           /*Gen*/
    std::setw(5)  << 0 <<                                           /*Daug*/
    std::setw(5)  << 3 <<                                           /*Zone*/
    std::setw(5)  << 0 <<                                           /*In*/
    std::setw(15) << particle.getKineticEnergy() <<                 /*KE*/
    std::setw(15) << (particle.getMomentum())[3] <<                 /*E*/
    std::setw(15) << (particle.getMomentum())[0] <<                 /*Px*/
    std::setw(15) << (particle.getMomentum())[1] <<                 /*Py*/
    std::setw(15) << (particle.getMomentum())[2] <<                 /*Pz*/
    std::setw(15) << -99. <<                                        /*Vx*/
    std::setw(15) << -99. <<                                        /*Vy*/
    std::setw(15) << -99.                                           /*Vz*/
    << std::endl;
}

void MyG4Analyser2::printResultsNtuple() {

  if (verboseLevel > 3) {
    G4cout << " >>> MyG4Analyser2::printResultsNtuple" << G4endl;
  }

  // Create one line of ASCII data. 
  // Several runs should create ntuple for data-analysis 
  G4cout <<
    std::setw(15) << bulletPDG <<
    std::setw(15) << bulletKineticEnergy <<
    std::setw(15) << int(eventNumber+0.01) <<
    std::setw(15) << Multiplicity << 
    std::setw(15) << ProtonNumber <<
    std::setw(15) << NeutronNumber << " " <<
    std::setw(15) << NucleonKinEnergy << " " <<
    std::setw(15) << ProtonKinEnergy << " " <<
    std::setw(15) << NeutronKinEnergy << " " <<
    std::setw(15) << PionNumber << " " <<
    std::setw(15) << PionKinEnergy << " " <<
    std::setw(15) << PhotonNumber << " " <<
    std::setw(15) << PhotonKinEnergy << G4endl;
}
