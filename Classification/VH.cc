// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "HepMC/GenParticle.h"
#include <iostream>

using HepMC::GenParticle;
using HepMC::GenVertex;

namespace Rivet {
  const double lep_cone_size = 0.1;

  const double jet_min_pT = 30*GeV;
  const double jet_max_eta = 4.7;

  bool cuts_0 [1];
  bool cuts_1 [3];
  bool cuts_2 [4];

  double m_sumw;

  /// @brief Add a short analysis description here
  class VH : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(VH);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections

      // the basic final-state projection: 
      // all final-state particles within 
      // the given eta acceptance
      const FinalState fs;
      declare(fs, "fs");

      m_sumw = 0.0;

      // the final-state particles declared above are clustered using FastJet with
      // the anti-kT algorithm and a jet-radius parameter 0.4
      // muons and neutrinos are excluded from the clustering
      //FastJets jets(fs, FastJets::ANTIKT, 0.4, JetAlg::Muons::NONE, JetAlg::Invisibles::NONE);
      //declare(jets, "jets");

      // FinalState of prompt photons and bare muons and electrons in the event
      PromptFinalState photons(Cuts::abspid == PID::PHOTON);
      PromptFinalState bare_leps(Cuts::abspid == PID::MUON || Cuts::abspid == PID::ELECTRON);
      declare(bare_leps, "bare_leps"); 
      // dress the prompt bare leptons with prompt photons within dR < 0.1
      DressedLeptons dressed_leps(photons, bare_leps, lep_cone_size);
      declare(dressed_leps, "leptons");

      //missing momentum
      MissingMomentum Met(fs);
      declare(Met, "Met");

      //book histograms
      book(H_PT, "H_PT", {0.0, 20.0, 45.0, 80.0, 120.0, 200.0});
      book(N_jets, "N_jets", {0, 1, 2, 3, 4});
      book(acc_H_PT, "acc_H_PT", {0.0, 20.0, 45.0, 80.0, 120.0, 200.0});
      book(acc_N_jets, "acc_N_jets", {0, 1, 2, 3, 4});

      book(acceptance, "acceptance", {0,1,2,3,4,5,6,7});
      book(lepton_no, "lepton_no", {0,1,2,3,4,5,6,7,8,9,10});
      book(V_mass, "V_mass", 50, 0, 200);
      book(V_PT, "V_PT", {0, 75, 150, 250, 400});
      book(acc_V_PT, "acc_V_PT", {0, 75, 150, 250, 400});
    }

    bool originateFrom(const Particle& p, const Particles& ptcls ) {
      const GenVertex* prodVtx = p.genParticle()->production_vertex();
      if (prodVtx == nullptr) return false;
      for (auto ancestor : Rivet::HepMCUtils::particles(prodVtx, HepMC::ancestors)){ 
          for ( auto part:ptcls ) 
              if ( ancestor==part.genParticle() ) return true;
          }
          return false; 
      }
     
    bool originateFrom(const Particle& p, const Particle& p2){
        Particles ptcls = {p2}; return originateFrom(p,ptcls);
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double weight = event.weights()[0];
      m_sumw += weight;

      //grab higgs
      Particle higgs;
      for (const GenParticle *ptcl : Rivet::HepMCUtils::particles(event.genEvent())){
          if (!PID::isHiggs(ptcl->pdg_id())) continue;
          higgs = Particle(ptcl);
      }

      //grav vector boson
      Particle V;
      for (const GenParticle *ptcl : Rivet::HepMCUtils::particles(event.genEvent())){
          if (((ptcl->pdg_id()) != 23) && ((ptcl->pdg_id()) != 24) && ((ptcl->pdg_id()) != -24)) continue;
          V = Particle(ptcl);
      }

      // grab jets not including particles from decay  
      Particles fs = apply<FinalState>(event, "fs").particles();
      Particles hadrons;
      for (const Particle &p : fs) {
          if (originateFrom(p,higgs) || originateFrom(p,V)) continue;
          hadrons += p;
      }
      FinalState fps_temp;
      FastJets jets(fps_temp, FastJets::ANTIKT, 0.4);
      jets.calc(hadrons); 
      Jets jets30 = jets.jetsByPt(Cuts::pT>jet_min_pT && Cuts::abseta<jet_max_eta);

      //fill histograms before acceptance cuts
      H_PT->fill(higgs.pT()/GeV);
      V_PT->fill(V.pT()/GeV);
      N_jets->fill(jets30.size());

      //start acceptance
      //Particles leptons = apply<DressedLeptons>(event, "leptons").particlesByPt(Cuts::pT>20*GeV && Cuts::abseta<2.5);
      Particles temp_leptons = apply<DressedLeptons>(event, "leptons").particlesByPt(); 
      lepton_no->fill(temp_leptons.size());

      //apply mLL cut
      FourMomentum LL = temp_leptons[0].momentum() + temp_leptons[1].momentum();
      if ( (LL.mass()<65*GeV) || (LL.mass()>95*GeV) ) vetoEvent;
 
      acc_H_PT->fill(higgs.pT()/GeV);
      acc_V_PT->fill(V.pT()/GeV);                                                                   
      acc_N_jets->fill(jets30.size());
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      double sf = m_sumw>0?1.0/m_sumw:1.0;
      for (auto hist:{H_PT, acc_H_PT, N_jets, acc_N_jets, V_PT, acc_V_PT}) scale(hist, sf);
    }
    //@}

    
    /// @name Histograms
    //@{
    //Histo1DPtr leptons_PT;
    Histo1DPtr H_PT;
    Histo1DPtr N_jets;
    Histo1DPtr acc_H_PT;
    Histo1DPtr acc_N_jets;
    Histo1DPtr acceptance;
    Histo1DPtr lepton_no;

    Histo1DPtr V_mass;
    Histo1DPtr V_PT;
    Histo1DPtr acc_V_PT;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(VH);


}
