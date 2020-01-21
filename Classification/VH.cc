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

  const bool acceptance_on = 0;

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

    /*
    void setCuts(const Event& event, Particles leptons, FourMomentum LL, double trailing_lep_mT, double higgs_mT) {
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
    */   

    void LepCuts_0(const Event& event, FourMomentum p_miss) {
	//find extra variables for cutting                                      

	//set cuts
	cuts_0[0] = p_miss.pT() < 170*GeV;	
    }


    void LepCuts_1(const Event& event, FourMomentum p_miss, Particles temp_leptons, double dphi) {
	//find extra variables for cutting                                      

	//set cuts
	FourMomentum V = p_miss + temp_leptons[0];
	cuts_1[0] = V.pT() < 100*GeV;

	for (const Particle p : temp_leptons) {
          if (((p.pid()) == 11) || ((p.pid()) == -11)) {
		cuts_1[1] = temp_leptons[0].pT() < 30*GeV;
          } else if (((p.pid()) == 13) || ((p.pid()) == -13)) {
		cuts_1[1] = temp_leptons[0].pT() < 25*GeV;
	  } else {
		cuts_1[1] = 1;
	  }
	cuts_1[2] = dphi > 2; 
        }
     }

    void LepCuts_2(const Event& event, Particles temp_leptons) {
	//find extra variables for cutting                                      

	//set cuts
	bool pair_found = 0;
      	Particle leading_lep;
      	Particle trailing_lep;
      	if (temp_leptons.size()>1) {
        	leading_lep = temp_leptons[0];
        	for (const Particle p: temp_leptons) {
          		if (!((p.pid())==(-leading_lep.pid()))) continue;
          		trailing_lep = p;
          		pair_found = 1;
          		break;
        	}	
      	}

	if (pair_found) {
	  cuts_2[0] = 0;
          FourMomentum LL = (leading_lep.momentum() + trailing_lep.momentum());
          V_mass->fill(LL.mass());
        
	  for (const Particle p : temp_leptons){
	    if (((p.pid()) == 13) || ((p.pid()) == -13)) {
	      cuts_2[1] = (LL.pT() < 50*GeV || LL.pT() > 150*GeV);
	    } else if (((p.pid()) == 11) || ((p.pid()) == -11)) {
	      cuts_2[1] = LL.pT() < 150*GeV;
	    } else {
	      cuts_2[1]= 1;
	    }
	  cuts_2[2]= (leading_lep.pT() < 20* GeV) || (trailing_lep.pT() < 20*GeV);
	  cuts_2[3] = (LL.mass() < 75*GeV) || (LL.mass() > 105*GeV);
	  }
	} else {
	  cuts_2[0] = 1;
	}
	 
    }

    /*
    //goes through cuts, if event failed any cuts, vetoEvent
    void applyVeto() {
        for (bool cut : cuts) {
		if (cut) vetoEvent;
	}
    }
    */

    bool checkCuts_0() {
	for (bool cut : cuts_0) {
		if (cut) return true;
	}
	return false;
    }

    bool checkCuts_1() {
        for (bool cut : cuts_1) {
                if (cut) return true;
        }
        return false;
    }

    bool checkCuts_2() {
        for (bool cut : cuts_2) {
                if (cut) return true;
        }
        return false;
    }

    /*
    //checks all cuts with one exception
    bool checkCutsException(int exception_index) {
    	for (int i = 0; i < 8; i++) {
		if (i != exception_index) {
			if (cuts[i] == true) return true;
		}
	}
	return false;	
    }
    */

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
      Particles temp_leptons = apply<PromptFinalState>(event, "bare_leps").particlesByPt(Cuts::pT>20*GeV && Cuts::abseta<2.5);  

      lepton_no->fill(temp_leptons.size());
      								  
      //fill histograms (if event not been vetoed)       
      if (temp_leptons.size() == 0) {
	FourMomentum p_miss = apply<MissingMomentum>(event, "Met").missingMomentum();
        LepCuts_0(event, p_miss);
	if (checkCuts_0()) vetoEvent;
        acceptance->fill(0);
      } else if (temp_leptons.size() == 1) {
	FourMomentum p_miss = apply<MissingMomentum>(event, "Met").missingMomentum();
	double dphi = deltaPhi(temp_leptons[0], p_miss);
	LepCuts_1(event, p_miss, temp_leptons, dphi);
	if (checkCuts_1()) vetoEvent;
        acceptance->fill(1);
      } else {
	LepCuts_2(event, temp_leptons);
	if (checkCuts_2()) vetoEvent;
        acceptance->fill(2);
      } 


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
    Histo1DPtr V_mass;
    Histo1DPtr V_PT;
    Histo1DPtr acc_V_PT;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(VH);


}
