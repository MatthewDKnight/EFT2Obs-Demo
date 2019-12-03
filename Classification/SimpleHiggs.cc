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

  const double leading_lep_min_pT = 25*GeV;
  const double trailing_lep_min_pT = 13*GeV;
  const double lep_max_eta = 2.5*GeV;
  const double dilep_min_mass = 12*GeV;
  const double dilep_min_pT = 30*GeV;
  const double trailing_lep_min_mT = 30*GeV;
  const double higgs_min_mT = 60*GeV;

  const double jet_min_pT = 30*GeV;
  const double jet_max_eta = 4.7;

  const double lep_cone_size = 0.1;

  const bool acceptance_on = 1;

  bool cuts [8];

  double m_sumw;

  /// @brief Add a short analysis description here
  class SimpleHiggs : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(SimpleHiggs);


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

      book(leading_lep_PT, "leading_lep_PT", 25,0,150);
      book(trailing_lep_PT, "trailing_lep_PT", 25,0,75);
      book(leading_lep_eta, "leading_lep_eta", 25,0,5);
      book(trailing_lep_eta, "trailing_lep_eta", 25,0,5);
      book(dilep_mass, "dilep_mass", 25,0,100);
      book(dilep_PT, "dilep_PT", 25,0,200);
      book(trailing_lep_mT_hist, "trailing_lep_mT", 25,0,120);
      book(higgs_mT_hist, "H_mT", 25,0,140);
      book(dphi_hist, "delta_phi", 25,0,3.14);
      book(cos_dphi, "cos_dphi", 25,-1,1);  

      book(acceptance, "acceptance", {0,1,2,3,4,5,6,7});
      book(lepton_no, "lepton_no", {0,1,2,3,4,5,6,7,8,9,10});
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

    void setCuts(const Event& event, Particles leptons, FourMomentum LL, double trailing_lep_mT, double higgs_mT) {
	//find extra variables for cutting                                      
	/*             
        Particles leptons = apply<DressedLeptons>(event, "leptons").particlesByPt(); 
        FourMomentum LL = (leptons[0].momentum() + leptons[1].momentum());
                                                                                     
        FourMomentum p_miss = apply<MissingMomentum>(event, "Met").missingMomentum();
                                                                                     
        double dphi = deltaPhi(leptons[1], p_miss);
        double trailing_lep_mT = sqrt(2*leptons[1].pT()*p_miss.pT() * (1-cos(dphi)));
                                                                                     
        dphi = deltaPhi(LL, p_miss);
        double higgs_mT = sqrt(2*LL.pT()*p_miss.pT() * (1-cos(dphi)));               
	*/

	//set cuts	
	cuts[0] = (leptons[0].charge() == leptons[1].charge()) || (leptons[0].abspid() == leptons[1].abspid());
	cuts[1] = leptons[0].pT() < leading_lep_min_pT;
	cuts[2] = leptons[1].pT() < trailing_lep_min_pT;
	cuts[3] = (leptons[0].abseta() > lep_max_eta || leptons[1].abseta() > lep_max_eta);
	cuts[4] = LL.mass() < dilep_min_mass;
	cuts[5] = LL.pT() < dilep_min_pT;
	cuts[6] = trailing_lep_mT < trailing_lep_min_mT;
	cuts[7] = higgs_mT < higgs_min_mT;
    }

    //goes through cuts, if event failed any cuts, vetoEvent
    void applyVeto() {
	for (bool cut : cuts) {
		if (cut) vetoEvent;
	}
    }

    bool checkCuts() {
	for (bool cut : cuts) {
		if (cut) return true;
	}
	return false;
    }

    //checks all cuts with one exception
    bool checkCutsException(int exception_index) {
    	for (int i = 0; i < 8; i++) {
		if (i != exception_index) {
			if (cuts[i] == true) return true;
		}
	}
	return false;	
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

      // grab jets not including particles from decay  
      Particles fs = apply<FinalState>(event, "fs").particles();
      Particles hadrons;
      for (const Particle &p : fs) {
          if (originateFrom(p,higgs)) continue;
          hadrons += p;
      }
      FinalState fps_temp;
      FastJets jets(fps_temp, FastJets::ANTIKT, 0.4);
      jets.calc(hadrons); 
      Jets jets30 = jets.jetsByPt(Cuts::pT>jet_min_pT && Cuts::abseta<jet_max_eta);
      
      //fill histograms before acceptance cuts
      H_PT->fill(higgs.pT()/GeV);
      N_jets->fill(jets30.size());

      //find extra variables for cutting
      Particles leptons = apply<DressedLeptons>(event, "leptons").particlesByPt();      
      FourMomentum LL = (leptons[0].momentum() + leptons[1].momentum());

      FourMomentum p_miss = apply<MissingMomentum>(event, "Met").missingMomentum();     
                                                                                   
      double dphi = deltaPhi(leptons[1], p_miss);
      double trailing_lep_mT = sqrt(2*leptons[1].pT()*p_miss.pT() * (1-cos(dphi)));
                                                                                   
      dphi = deltaPhi(LL, p_miss);
      double higgs_mT = sqrt(2*LL.pT()*p_miss.pT() * (1-cos(dphi)));               



      //angular separation of two leptons
      dphi = deltaPhi(leptons[0], leptons[1]);

      //find which cuts event passes/fails
      setCuts(event, leptons, LL, trailing_lep_mT, higgs_mT);

      //fill each histogram with events that pass all but the one cut associated with variable
      if (!checkCutsException(1)) leading_lep_PT->fill(leptons[0].pT()/GeV);
      if (!checkCutsException(2)) trailing_lep_PT->fill(leptons[1].pT()/GeV);
      if (!checkCutsException(3)) leading_lep_eta->fill(leptons[0].abseta());
      if (!checkCutsException(3)) trailing_lep_eta->fill(leptons[1].abseta());
      if (!checkCutsException(4)) dilep_mass->fill(LL.mass()/GeV);
      if (!checkCutsException(5)) dilep_PT->fill(LL.pT()/GeV);
      if (!checkCutsException(6)) trailing_lep_mT_hist->fill(trailing_lep_mT/GeV);
      if (!checkCutsException(7)) higgs_mT_hist->fill(higgs_mT/GeV);
      dphi_hist->fill(dphi);
      cos_dphi->fill(cos(dphi));

      //vetoEvent if event fails any of the cuts
      if (checkCuts()) vetoEvent;    
      
      //fill histograms (if event not been vetoed)       
      acc_H_PT->fill(higgs.pT()/GeV);                                                                   
      acc_N_jets->fill(jets30.size());                                                
    }

    /// Normalise histograms etc., after the run
    void finalize() {
      double sf = m_sumw>0?1.0/m_sumw:1.0;
      for (auto hist:{H_PT, acc_H_PT, N_jets, acc_N_jets}) scale(hist, sf);

      /* 
      normalize(H_PT); // normalize to unity
      scale(H_PT, crossSection()/picobarn/sumOfWeights()); // norm to cross section
      normalize(N_jets);
      scale(N_jets, crossSection()/picobarn/sumOfWeights());

      normalize(acc_H_PT); // normalize to unity
      scale(acc_H_PT, crossSection()/picobarn/sumOfWeights()); // norm to cross section
      normalize(acc_N_jets);
      scale(acc_N_jets, crossSection()/picobarn/sumOfWeights());
      
     
      higgs_PT->normalize(crossSection()/picobarn, true);
      N_jets->normalize(crossSection()/picobarn, true);

      higgs_PT_acceptance->normalize(crossSection()/picobarn, true);
      N_jets_acceptance->normalize(crossSection()/picobarn, true);
      
      H_PT->normalize(82.52, true);
      N_jets->normalize(82.52, true);

      acc_H_PT->normalize(82.52, true);
      acc_N_jets->normalize(82.52, true);
      */
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

    Histo1DPtr leading_lep_PT;
    Histo1DPtr trailing_lep_PT;
    Histo1DPtr leading_lep_eta;
    Histo1DPtr trailing_lep_eta;
    Histo1DPtr dilep_mass;
    Histo1DPtr dilep_PT;
    Histo1DPtr trailing_lep_mT_hist;
    Histo1DPtr higgs_mT_hist;
    Histo1DPtr dphi_hist;
    Histo1DPtr cos_dphi;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(SimpleHiggs);


}
