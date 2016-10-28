// -*- C++ -*-
//
// Package:    Acceptance/TagAndProbeAnalyzer
// Class:      TagAndProbeAnalyzer
// 
/**\class TagAndProbeAnalyzer TagAndProbeAnalyzer.cc Acceptance/TagAndProbeAnalyzer/plugins/TagAndProbeAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Tyler Henry Ruggles
//         Created:  Mon, 29 Feb 2016 10:32:40 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

// Added
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TTree.h"
#include "TH1D.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Math/interface/deltaR.h"

// Trigger stuff...
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/Tau.h"

#include <algorithm>
//#include <map>

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class TagAndProbeAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit TagAndProbeAnalyzer(const edm::ParameterSet&);
      ~TagAndProbeAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<std::vector<reco::GenJet>> genHadronicTausToken_;
      edm::EDGetTokenT<std::vector<reco::GenJet>> genElectronicTausToken_;
      edm::EDGetTokenT<std::vector<reco::GenJet>> genMuonicTausToken_;
      edm::EDGetTokenT<std::vector<PileupSummaryInfo>> puToken_;
      edm::EDGetTokenT<std::vector<pat::Tau>> tauToken_;
      edm::EDGetTokenT<std::vector<pat::Muon>> muonToken_;
      edm::EDGetTokenT<std::vector<pat::Electron>> electronToken_;
      edm::EDGetTokenT<std::vector<pat::Jet>> jetToken_;
      edm::EDGetTokenT<std::vector<pat::MET>> metToken_;
      edm::EDGetTokenT<std::vector<reco::Vertex>> vertexToken_;
      edm::EDGetTokenT<edm::TriggerResults> triggerToken1_;
      edm::EDGetTokenT<edm::TriggerResults> triggerToken2_;
      edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjectsToken_;
      edm::EDGetTokenT<BXVector<l1t::Tau>> stage2TauToken1_;
      edm::EDGetTokenT<BXVector<l1t::Tau>> stage2TauToken2_;
      edm::EDGetTokenT<std::vector<reco::GenParticle>> genToken_;
      // l1 extras
      //edm::EDGetTokenT<edm::TriggerResults> triggerToken_;

      TTree *tree;
      TH1D *nEvents;
      TH1D *cutFlow;
      double eventD;
      float run, lumi, nTruePU, nvtx, nvtxCleaned, IsoMu20, IsoMu22, IsoMu24,
        IsoMu27, IsoMu21MediumIsoTau32, IsoMu19MediumIsoTau32,IsoMu19MediumCombinedIsoTau32,IsoMu21MediumCombinedIsoTau32, IsoMu19LooseIsoTau20_SingleL1, IsoMu19LooseIsoTau20, IsoMu21LooseIsoTau20_SingleL1, IsoMu17LooseIsoTau20_SingleL1, IsoMu17LooseIsoTau20, TrigPass, mPt, mEta, mPhi,
        tPt, tEta, tPhi, tMVAIsoVLoose, tMVAIsoLoose, tMVAIsoMedium, 
        tMVAIsoTight, tMVAIsoVTight, m_vis, transMass, SS,
        leptonDR, mTrigMatch, tTrigMatch, mL1Match, tL1Match,
        t_gen_match;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
TagAndProbeAnalyzer::TagAndProbeAnalyzer(const edm::ParameterSet& iConfig) :
    genHadronicTausToken_(consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("hadronSrc"))),
    genElectronicTausToken_(consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("tauElectronSrc"))),
    genMuonicTausToken_(consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("tauMuonSrc"))),
    puToken_(consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("puSrc"))),

    tauToken_(consumes<std::vector<pat::Tau>>(iConfig.getParameter<edm::InputTag>("tauSrc"))),
    muonToken_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muonSrc"))),
    electronToken_(consumes<std::vector<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electronSrc"))),
    jetToken_(consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jetSrc"))),
    metToken_(consumes<std::vector<pat::MET>>(iConfig.getParameter<edm::InputTag>("metSrc"))),
    vertexToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("pvSrc"))),
    triggerToken1_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerSrc1"))),
    triggerToken2_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerSrc2"))),
    triggerObjectsToken_(consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerObjectsSrc"))),
    stage2TauToken1_(consumes<BXVector<l1t::Tau>>(iConfig.getParameter<edm::InputTag>("stage2TauSrc1"))),
    stage2TauToken2_(consumes<BXVector<l1t::Tau>>(iConfig.getParameter<edm::InputTag>("stage2TauSrc2"))),
    genToken_(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genSrc")))
{
   //now do what ever initialization is needed
   //usesResource("TFileService");
   edm::Service<TFileService> fs;
   TFileDirectory subDir = fs->mkdir( "tagAndProbe" );
   nEvents = subDir.make<TH1D>("nEvents","nEvents",1,-0.5,0.5);
   cutFlow = subDir.make<TH1D>("cutFlow","cutFlow",10,-0.5,9.5);
   tree = subDir.make<TTree>("Ntuple","My T-A-P Ntuple");
   tree->Branch("run",&run,"run/F");
   tree->Branch("lumi",&lumi,"lumi/F");
   tree->Branch("eventD",&eventD,"eventD/D");
   tree->Branch("nTruePU",&nTruePU,"nTruePU/F");
   tree->Branch("nvtx",&nvtx,"nvtx/F");
   tree->Branch("nvtxCleaned",&nvtxCleaned,"nvtxCleaned/F");
   tree->Branch("IsoMu20",&IsoMu20,"IsoMu20/F");
   tree->Branch("IsoMu22",&IsoMu22,"IsoMu22/F");
   tree->Branch("IsoMu24",&IsoMu24,"IsoMu24/F");
   tree->Branch("IsoMu27",&IsoMu27,"IsoMu27/F");
   tree->Branch("IsoMu21MediumIsoTau32",&IsoMu21MediumIsoTau32,"IsoMu21MediumIsoTau32/F");
   tree->Branch("IsoMu19MediumIsoTau32",&IsoMu19MediumIsoTau32,"IsoMu19MediumIsoTau32/F");
   tree->Branch("IsoMu19MediumCombinedIsoTau32",&IsoMu19MediumCombinedIsoTau32,"IsoMu19MediumCombinedIsoTau32/F");
   tree->Branch("IsoMu21MediumCombinedIsoTau32",&IsoMu21MediumCombinedIsoTau32,"IsoMu21MediumCombinedIsoTau32/F");
   tree->Branch("IsoMu19LooseIsoTau20_SingleL1", &IsoMu19LooseIsoTau20_SingleL1, "IsoMu19LooseIsoTau20_SingleL1/F");
   tree->Branch("IsoMu19LooseIsoTau20", &IsoMu19LooseIsoTau20, "IsoMu19LooseIsoTau20/F");
   tree->Branch("IsoMu21LooseIsoTau20_SingleL1", &IsoMu21LooseIsoTau20_SingleL1, "IsoMu21LooseIsoTau20_SingleL1/F");
   tree->Branch("IsoMu17LooseIsoTau20_SingleL1", &IsoMu17LooseIsoTau20_SingleL1,"IsoMu17LooseIsoTau20_SingleL1/F");
   tree->Branch("IsoMu17LooseIsoTau20", &IsoMu17LooseIsoTau20, "IsoMu17LooseIsoTau20/F");
   tree->Branch("TrigPass",&TrigPass,"TrigPass/F");
   tree->Branch("mPt",&mPt,"mPt/F");
   tree->Branch("mEta",&mEta,"mEta/F");
   tree->Branch("mPhi",&mPhi,"mPhi/F");
   tree->Branch("mTrigMatch",&mTrigMatch,"mTrigMatch/F");
   tree->Branch("mL1Match",&mL1Match,"mL1Match/F");
   tree->Branch("tPt",&tPt,"tPt/F");
   tree->Branch("tEta",&tEta,"tEta/F");
   tree->Branch("tPhi",&tPhi,"tPhi/F");
   tree->Branch("t_gen_match",&t_gen_match,"t_gen_match/F");
   tree->Branch("tMVAIsoVLoose",&tMVAIsoVLoose,"tMVAIsoVLoose/F");
   tree->Branch("tMVAIsoLoose",&tMVAIsoLoose,"tMVAIsoLoose/F");
   tree->Branch("tMVAIsoMedium",&tMVAIsoMedium,"tMVAIsoMedium/F");
   tree->Branch("tMVAIsoTight",&tMVAIsoTight,"tMVAIsoTight/F");
   tree->Branch("tMVAIsoVTight",&tMVAIsoVTight,"tMVAIsoVTight/F");
   tree->Branch("tTrigMatch",&tTrigMatch,"tTrigMatch/F");
   tree->Branch("tL1Match",&tL1Match,"tL1Match/F");
   tree->Branch("leptonDR",&leptonDR,"leptonDR/F");
   tree->Branch("m_vis",&m_vis,"m_vis/F");
   tree->Branch("transMass",&transMass,"transMass/F");
   tree->Branch("SS",&SS,"SS/F");

}


TagAndProbeAnalyzer::~TagAndProbeAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
TagAndProbeAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;

    // First thing, fill the nEvents
    nEvents->Fill(0.);
    cutFlow->Fill(0., 1.);

    edm::Handle<std::vector<reco::Vertex>> vertices;   
    iEvent.getByToken(vertexToken_, vertices);
    if (vertices->empty()) return; // skip the event if no PV found
    const reco::Vertex &PV = vertices->front();
    if (PV.ndof() < 4) return; // bad vertex
    nvtx = vertices.product()->size();
    for (const reco::Vertex &vertex : *vertices)
        if (!vertex.isFake()) ++nvtxCleaned;

    // Get the number of true events
    // This is used later for pile up reweighting
    edm::Handle<std::vector<PileupSummaryInfo>> puInfo;   
    iEvent.getByToken(puToken_, puInfo);
    if (puInfo.isValid()) {
        if (puInfo->size() > 0) {
            nTruePU = puInfo->at(1).getTrueNumInteractions();
        }
    }

    cutFlow->Fill(1., 1.); // Good vertex



    edm::Handle<std::vector<pat::Muon>> muons; 
    iEvent.getByToken(muonToken_, muons);
    // Storage for the "best" muon
    pat::Muon bestMuon;
    int passingMuons = 0;

    for (const pat::Muon &mu : *muons) {
        if (mu.pt() < 25 || fabs(mu.eta()) > 2.1 || !mu.isMediumMuon()) continue;
        float mIso = (mu.pfIsolationR04().sumChargedHadronPt
            + TMath::Max(0., mu.pfIsolationR04().sumNeutralHadronEt
            + mu.pfIsolationR04().sumPhotonEt
            - 0.5*mu.pfIsolationR04().sumPUPt))
            /mu.pt();
        if (mIso > 0.1) continue;
        passingMuons++;
        bestMuon = mu;
    }
    // Require strictly 1 muon
    if (passingMuons == 0) return;
    cutFlow->Fill(2., 1.);
    // Extra lepton veto (muons)
    if (passingMuons > 1) return;
    cutFlow->Fill(3., 1.);


    edm::Handle<std::vector<pat::Electron>> electrons;   
    iEvent.getByToken(electronToken_, electrons);
    int passingElectrons = 0;
    for (const pat::Electron &el : *electrons) {
        if (el.pt() < 20 || fabs(el.eta()) > 2.1 || 
            el.electronID("mvaEleID-Spring15-25ns-nonTrig-V1-wp90") < 0.5) continue;
        float eIso = (el.pfIsolationVariables().sumChargedHadronPt + TMath::Max(
            el.pfIsolationVariables().sumNeutralHadronEt +
            el.pfIsolationVariables().sumPhotonEt -
            0.5 * el.pfIsolationVariables().sumPUPt, 0.0)) / el.pt();
        if (eIso > 0.1) continue;
        passingElectrons++;
    }
    // Extra lepton veto (electrons)
    if (passingElectrons > 0) return;
    cutFlow->Fill(4., 1.);


    edm::Handle<std::vector<pat::Tau>> taus;   
    iEvent.getByToken(tauToken_, taus);
    // Storage for the "best" muon
    pat::Tau bestTau;
    int passingTaus = 0;
    for (const pat::Tau &tau : *taus) {
        if (tau.pt() < 20 || fabs(tau.eta()) > 2.1 || tau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT") < 0.5
            || tau.tauID("decayModeFinding") < 0.5 || tau.tauID("againstElectronTightMVA6") < 0.5
            || tau.tauID("againstMuonTight3") < 0.5 || fabs(tau.charge())!=1) continue;
        // No extra lepton vetoes rely on tau number
        // so, if multiple taus, choose highest pt one.
        // Do trigger matching later
        if (passingTaus == 0) bestTau = tau;
        else if (tau.pt() > bestTau.pt())
            bestTau = tau;
        ++passingTaus;
    }
    // Tau study so...
    if (passingTaus == 0) return;
    cutFlow->Fill(5., 1.);


    // Check for non-overlapping bjets
    // using Medium CISV value of 0.8
    edm::Handle<std::vector<pat::Jet>> jets;   
    iEvent.getByToken(jetToken_, jets);
    bool btagged = false;
    for (const pat::Jet &j : *jets) {
        if (j.pt() < 20 || fabs(j.eta()) > 2.4 || j.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") < 0.8) continue;
        if (deltaR(j, bestMuon) > 0.5) continue;
        if (deltaR(j, bestTau) > 0.5) continue;
        btagged = true;
    }
    if (btagged) return;
    cutFlow->Fill(6., 1.);


    // Save our best tau and muon variables
    mPt = bestMuon.pt();
    mEta = bestMuon.eta();
    mPhi = bestMuon.phi();
    tPt = bestTau.pt();
    tEta = bestTau.eta();
    tPhi = bestTau.phi();
    tMVAIsoVLoose = bestTau.tauID("byVLooseIsolationMVArun2v1DBoldDMwLT");
    tMVAIsoLoose = bestTau.tauID("byLooseIsolationMVArun2v1DBoldDMwLT");
    tMVAIsoMedium = bestTau.tauID("byMediumIsolationMVArun2v1DBoldDMwLT");
    tMVAIsoTight = bestTau.tauID("byTightIsolationMVArun2v1DBoldDMwLT");
    tMVAIsoVTight = bestTau.tauID("byVTightIsolationMVArun2v1DBoldDMwLT");
    leptonDR = deltaR( bestMuon, bestTau );


    // Get MET for transverse mass calculation 
    edm::Handle<std::vector<pat::MET>> mets;   
    iEvent.getByToken(metToken_, mets);
    const pat::MET &met = mets->front();
    transMass = TMath::Sqrt( 2. * bestMuon.pt() * met.pt() * (1. - TMath::Cos( bestMuon.phi() - met.phi())));


    // Get Visible Mass
    TLorentzVector l1 = TLorentzVector( 0., 0., 0., 0. );
    l1.SetPtEtaPhiM( bestMuon.pt(), bestMuon.eta(),
        bestMuon.phi(), bestMuon.mass() );
    TLorentzVector l2 = TLorentzVector( 0., 0., 0., 0. );
    l2.SetPtEtaPhiM( bestTau.pt(), bestTau.eta(),
        bestTau.phi(), bestTau.mass() );
    m_vis = (l1 + l2).M();


    // Same sign comparison
    if (bestMuon.charge() + bestTau.charge() == 0) SS = 0;
    else SS = 1;


    eventD = iEvent.eventAuxiliary().event();
    lumi = iEvent.eventAuxiliary().luminosityBlock();
    run = iEvent.eventAuxiliary().run();

    //std::cout << "Run: " <<run<< "  Evt: " <<eventD<< "  Lumi: " <<lumi<<std::endl;
    printf("Run: %.0f    Evt: %.0f   Lumi: %.0f\n", run, eventD, lumi);
    


    edm::Handle<edm::TriggerResults> trigger1;   
    iEvent.getByToken(triggerToken1_, trigger1);
    edm::Handle<edm::TriggerResults> trigger2;   
    iEvent.getByToken(triggerToken2_, trigger2);
    edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
    iEvent.getByToken(triggerObjectsToken_, triggerObjects);

    IsoMu20 = 0.;
    IsoMu22 = 0.;
    IsoMu24 = 0.;
    IsoMu27 = 0.;
    IsoMu21MediumIsoTau32 = 0.;
    TrigPass = 0.;

    edm::Handle<edm::TriggerResults> triggerResults;
    if (trigger1.isValid()) {
        //std::cout << "Trigger Objects   HLT is valid" << std::endl;
        triggerResults = trigger1;
    }
    else if (trigger2.isValid()) {
        //std::cout << "Trigger Objects   HLT2 is valid" << std::endl;
        triggerResults = trigger2;
    }

    std::vector<std::string> usedPaths;
    const edm::TriggerNames &names = iEvent.triggerNames(*triggerResults);
    // See https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD2014#Trigger
    for (unsigned int i = 0, n = triggerResults->size(); i < n; ++i) {
        //std::cout << names.triggerName(i) << std::endl;
        if (names.triggerName(i).find("HLT_IsoMu20_v") != std::string::npos) {
            if (triggerResults->accept(i)) IsoMu20 = 1; TrigPass = 1; usedPaths.push_back( names.triggerName(i) );
        }
        if (names.triggerName(i).find("HLT_IsoMu22_v") != std::string::npos) {
            if (triggerResults->accept(i)) IsoMu22 = 1; TrigPass = 1; usedPaths.push_back( names.triggerName(i) );
        }
        if (names.triggerName(i).find("HLT_IsoMu24_v") != std::string::npos) {
            if (triggerResults->accept(i)) IsoMu24 = 1; TrigPass = 1; usedPaths.push_back( names.triggerName(i) );
        }
        if (names.triggerName(i).find("HLT_IsoMu27_v") != std::string::npos) {
            if (triggerResults->accept(i)) IsoMu27 = 1; TrigPass = 1; usedPaths.push_back( names.triggerName(i) );
        }
        if (names.triggerName(i).find("HLT_IsoMu21_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg_v") != std::string::npos) {
            if (triggerResults->accept(i)) IsoMu21MediumIsoTau32 = 1; TrigPass = 1; usedPaths.push_back( names.triggerName(i) );
        }
        //adding more filters for tauPOG
	if (names.triggerName(i).find("HLT_IsoMu19_eta2p1_MediumIsoPFTau32_Trk1_eta2p1_Reg_v") != std::string::npos) {
            if (triggerResults->accept(i)) IsoMu19MediumIsoTau32 = 1; TrigPass = 1; usedPaths.push_back( names.triggerName(i) );
        }
        if (names.triggerName(i).find("HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1_v") != std::string::npos) {
            if (triggerResults->accept(i)) IsoMu19LooseIsoTau20_SingleL1 = 1; TrigPass = 1; usedPaths.push_back( names.triggerName(i) );
        }
        if (names.triggerName(i).find("HLT_IsoMu19_eta2p1_LooseIsoPFTau20_v") != std::string::npos) {
            if (triggerResults->accept(i)) IsoMu19LooseIsoTau20 = 1; TrigPass = 1; usedPaths.push_back( names.triggerName(i) );
        }
        if (names.triggerName(i).find("HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1_v") != std::string::npos) {
            if (triggerResults->accept(i))  IsoMu21LooseIsoTau20_SingleL1 = 1; TrigPass = 1; usedPaths.push_back( names.triggerName(i) );
        }
        if (names.triggerName(i).find("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_SingleL1_v") != std::string::npos) {
            if (triggerResults->accept(i))  IsoMu17LooseIsoTau20_SingleL1 = 1; TrigPass = 1; usedPaths.push_back( names.triggerName(i) );
        }
        if (names.triggerName(i).find("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v") != std::string::npos) {
            if (triggerResults->accept(i))  IsoMu17LooseIsoTau20 = 1; TrigPass = 1; usedPaths.push_back( names.triggerName(i) );
        }
        if (names.triggerName(i).find("HLT_IsoMu19_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v") != std::string::npos) { 
            if (triggerResults->accept(i))  IsoMu19MediumCombinedIsoTau32 = 1; TrigPass = 1; usedPaths.push_back( names.triggerName(i) );
        }
        if (names.triggerName(i).find("HLT_IsoMu21_eta2p1_MediumCombinedIsoPFTau32_Trk1_eta2p1_Reg_v") != std::string::npos) {
            if (triggerResults->accept(i))  IsoMu21MediumCombinedIsoTau32 = 1; TrigPass = 1; usedPaths.push_back( names.triggerName(i) );
        }
    }



    // Do trigger object matching
    // for the moment, just record the number
    // of times our 'best' objects match
    // this can be expanded later to indivual trigs if necessary
    mTrigMatch = 0;
    tTrigMatch = 0;
    for (pat::TriggerObjectStandAlone obj : *triggerObjects) { // note: not "const &" since we want to call unpackPathNames
        obj.unpackPathNames(names);
        std::vector<std::string> pathNamesLast = obj.pathNames(true);
        // pathNamesLast = vector of flags, if this object was used in the final 
        // filter of a succeeding HLT path resp. in a succeeding 
        // condition of a succeeding L1 algorithm
        for (unsigned h = 0, n = pathNamesLast.size(); h < n; ++h) {
            if (std::find( usedPaths.begin(), usedPaths.end(), pathNamesLast[h]) != usedPaths.end()) {
                //std::cout << " ---  " << pathNamesLast[h] << std::endl;
                //std::cout << "\tTrigger object:  pt " << obj.pt() << ", eta " << obj.eta() << ", phi " << obj.phi() << std::endl;
                float drMu = deltaR( bestMuon, obj );
                float drTau = deltaR( bestTau, obj );
                //std::cout << "\tbestMuon dR: " << drMu << std::endl;
                //std::cout << "\tbestTau dR: " << drTau << std::endl;
                if (drMu < 0.5) ++mTrigMatch;
                if (drTau < 0.5) ++tTrigMatch;
            }
        }
    }


    // Do l1extra object matching
    // Make sure to only consider l1 objects from the
    // intended BX.  That is the size(0) and
    // at(0,i) notation
    edm::Handle<BXVector<l1t::Tau>> l1Taus1; 
    iEvent.getByToken(stage2TauToken1_, l1Taus1);
    edm::Handle<BXVector<l1t::Tau>> l1Taus2; 
    iEvent.getByToken(stage2TauToken2_, l1Taus2);

    edm::Handle<BXVector<l1t::Tau>> l1Taus; 
    if (l1Taus1.isValid()) {
        //std::cout << "Trigger Objects   HLT is valid" << std::endl;
        l1Taus = l1Taus1;
    }
    else if (l1Taus2.isValid()) {
        //std::cout << "Trigger Objects   HLT2 is valid" << std::endl;
        l1Taus = l1Taus2;
    }
    
    tL1Match = 0;
    if (l1Taus.isValid()) {
        std::cout << "L1 Extras is valid" << std::endl;
        for (size_t i = 0; i < l1Taus->size(0); ++i) {
            const l1t::Tau &l1Tau = l1Taus->at(0,i);
            // skip l1Tau if it's low pt b/c the trigger we want
            // actual results for is seeded by
            // L1_DoubleIsoTau28
            if (l1Tau.hwIso()<1 || l1Tau.pt()<27.5) continue; // hardware Iso bit
            float drTau = deltaR( bestTau, l1Tau );
            //std::cout << " - " << i << " L1Tau pt: " << l1Tau.pt() 
            //<< " Iso: " << l1Tau.hwIso() << " dr: " << drTau << std::endl;
            if (drTau < 0.5) ++tL1Match;
        }    
    } // end l1Taus



    // Check for overlapping Gen Taus
    // with reconstructed gen taus of 3 types
    // and normal gen particles
    edm::Handle<std::vector<reco::GenJet>> genHTaus;   
    iEvent.getByToken(genHadronicTausToken_, genHTaus);
    edm::Handle<std::vector<reco::GenJet>> genETaus;   
    iEvent.getByToken(genElectronicTausToken_, genETaus);
    edm::Handle<std::vector<reco::GenJet>> genMTaus;   
    iEvent.getByToken(genMuonicTausToken_, genMTaus);
    edm::Handle<std::vector<reco::GenParticle>> genParticles; 
    iEvent.getByToken(genToken_, genParticles);
    t_gen_match = -1;
    if (genParticles.isValid()) {
        // Find the closest gen particle to our candidate
        if ( genParticles->size() > 0 ) {
            reco::GenParticle closest = genParticles->at(0);
            float closestDR = 999;
            // The first two codes are based off of matching to true electrons/muons
            // Find the closest gen particle...
            for(size_t m = 0; m != genParticles->size(); ++m) {
                reco::GenParticle genp = genParticles->at(m);
                float tmpDR = deltaR( bestTau.p4(), genp.p4() );
                if ( tmpDR < closestDR ) { closest = genp; closestDR = tmpDR; }
            }
            float genID = abs(closest.pdgId());

            // Loop over all versions of gen taus and find closest one
            float closestDR_HTau = 999;
            float closestDR_ETau = 999;
            float closestDR_MTau = 999;
            if ( genHTaus->size() > 0 ) {
                for (size_t j = 0; j != genHTaus->size(); ++j) {
                    float tmpDR = deltaR( bestTau.p4(), genHTaus->at(j).p4() );
                    if (tmpDR < closestDR_HTau) closestDR_HTau = tmpDR;
                }
            }
            if ( genETaus->size() > 0 ) {
                for (size_t j = 0; j != genETaus->size(); ++j) {
                    float tmpDR = deltaR( bestTau.p4(), genETaus->at(j).p4() );
                    if (tmpDR < closestDR_ETau) closestDR_ETau = tmpDR;
                }
            }
            if ( genMTaus->size() > 0 ) {
                for (size_t j = 0; j != genMTaus->size(); ++j) {
                    float tmpDR = deltaR( bestTau.p4(), genMTaus->at(j).p4() );
                    if (tmpDR < closestDR_MTau) closestDR_MTau = tmpDR;
                }
            }

            // Now return the value based on which object is closer, the closest
            // single gen particle, or the rebuild gen taus
            // The first two codes are based off of matching to true electrons/muons
            float closestGetTau = TMath::Min(closestDR_ETau, closestDR_MTau);
            if (closestDR_HTau < closestGetTau) closestGetTau = closestDR_HTau;

            // Make sure we don't overwrite a proper value
            if (closestDR < closestGetTau && genID == 11 && closest.pt() > 8
                    && closest.statusFlags().isPrompt() && closestDR < 0.2 )
                        t_gen_match = 1.0;
            else if (closestDR < closestGetTau && genID == 13 && closest.pt() > 8
                    && closest.statusFlags().isPrompt() && closestDR < 0.2 )
                        t_gen_match = 2.0;
            // Other codes based off of not matching previous 2 options
            // as closest gen particle, retruns based on closest rebuilt gen tau
            else if (closestDR_ETau < 0.2 && closestDR_ETau < TMath::Min(closestDR_MTau, 
                    closestDR_HTau)) t_gen_match = 3.0;
            else if (closestDR_MTau < 0.2 && closestDR_MTau < TMath::Min(closestDR_ETau, 
                    closestDR_HTau)) t_gen_match = 4.0;
            else if (closestDR_HTau < 0.2 && closestDR_HTau < TMath::Min(closestDR_ETau, 
                    closestDR_MTau)) t_gen_match = 5.0;
            else t_gen_match = 6.0; // No match, return 6 for "fake tau"
        }
    }



    // Denominator section first, just check if there's the correct #
    // of the associated leptons
//    if (hTaus->size() >= 2) TauTauD = 1;
//    if (hTaus->size() >= 1) {
//        if(eTaus->size() >= 1) ETauD = 1;
//        if(mTaus->size() >= 1) MuTauD = 1;
//    }
//    if(eTaus->size() >= 1 && mTaus->size() >= 1) EMuD = 1;
//    if(mTaus->size() >= 2) MuMuD = 1;
//    
//    float nLeptons = 0;
//    nLeptons += hTaus->size() + eTaus->size() + mTaus->size();
//    if (nLeptons >= 3) threeLeptons = 1;
//    nLooseTaus = hTaus->size();
//    nLooseElec = eTaus->size();
//    nLooseMu = mTaus->size();
//
//
//    size_t ETau_e = 0;
//    size_t ETau_t = 0;
//    size_t MuTau_m = 0;
//    size_t MuTau_t = 0;
//    size_t EMu_e = 0;
//    size_t EMu_m = 0;
//    size_t TauTau_t1 = 0;
//    size_t TauTau_t2 = 0;
//    size_t MuMu_m1 = 0;
//    size_t MuMu_m2 = 0;
//
//    //uint32_t nHTaus = hTaus->size();
//    //if (hTaus->size() > 0)// std::cout << " --- N Hadronic Taus: "<<nHTaus<<std::endl;
//    //{
//    //const std::vector<reco::GenJet> nTaus_ = hTaus.product();
//    std::vector< float > pts;
//    for (const reco::GenJet &tau : *hTaus) {
//        pts.push_back( tau.pt() );
//        if ( TMath::Abs(tau.eta()) < 2.3 && tau.pt() > 20 ) {
//            MuTau_t += 1;
//            ETau_t += 1;}
//        if ( TMath::Abs(tau.eta()) < 2.1 && tau.pt() > 40 ) TauTau_t1 += 1;
//        if ( TMath::Abs(tau.eta()) < 2.1 && tau.pt() > 30 ) TauTau_t2 += 1;
//    }
//
//    if ( pts.size() > 0 ) tauPt1 = pts.at(0);
//    if ( pts.size() > 1 ) tauPt2 = pts.at(1);
//    if ( pts.size() > 2 ) tauPt3 = pts.at(2);
//
//    //}
//    //uint32_t nETaus = eTaus->size();
//    //if (nETaus > 0)// std::cout << " ### N Electronic Taus: "<<nETaus<<std::endl;
//    //{
//    for (const reco::GenJet &ele : *eTaus) {
//        if ( TMath::Abs(ele.eta()) < 2.5 && ele.pt() > 13 ) EMu_e += 1;
//        if ( TMath::Abs(ele.eta()) < 2.1 && ele.pt() > 24 ) ETau_e += 1;
//    //}
//    }
//    //uint32_t nMTaus = mTaus->size();
//    //if (nMTaus > 0)// std::cout << " *** N Muonic Taus: "<<nMTaus<<std::endl;
//    //{
//    for (const reco::GenJet &mu : *mTaus) {
//        if ( TMath::Abs(mu.eta()) < 2.4 && mu.pt() > 10 ) EMu_m += 1;
//        if ( TMath::Abs(mu.eta()) < 2.1 && mu.pt() > 19 ) MuTau_m += 1;
//        if ( TMath::Abs(mu.eta()) < 2.4 && mu.pt() > 20 ) MuMu_m1 += 1; // # of mu passing "leading" cut
//        if ( TMath::Abs(mu.eta()) < 2.4 && mu.pt() < 20 && mu.pt() > 10 ) MuMu_m2 += 1; // # of mu passing only "trailing" cut
//    }
//    //}
//    
//    // Check if we have matches
//    if (ETau_e > 0 && ETau_t > 0) ETauPass = 1;
//    if (MuTau_m > 0 && MuTau_t > 0) MuTauPass = 1;
//    if (EMu_e > 0 && EMu_m > 0) EMuPass = 1;
//    if (TauTau_t1 > 1) TauTauPass = 1;
//    if (TauTau_t1 == 1 && TauTau_t2 > 0) TauTau4030Pass = 1;
//    if (MuMu_m1 > 1) MuMuPass = 1;
//    if (MuMu_m1 > 0 && MuMu_m2 > 0) MuMuPass = 1;
//    
//
//    if (invmass.size() == 2) {
//        //std::cout << "Len InvMass: "<<invmass.size()<<std::endl;
//        TLorentzVector diLep = invmass[0];
//        diLep += invmass[1];
//        //std::cout << "m(ll) " << diLep.M() << std::endl;
//        genMass = diLep.M();
//
//    }
    
 
    
    

    //LogInfo("Demo") << "number of gen taus "<<nGenTaus;
    //std::cout << genTaus << std::endl;

    tree->Fill();


#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
TagAndProbeAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TagAndProbeAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TagAndProbeAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TagAndProbeAnalyzer);
