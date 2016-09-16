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
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include <map>

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
      // l1 extras
      //edm::EDGetTokenT<edm::TriggerResults> triggerToken_;

      TTree *tree;
      TH1D *nEvents;
      //float run, lumi;
      double eventD;
      float nTruePU;
      float run, lumi, nvtx, nvtxCleaned, IsoMu20, IsoMu22, IsoMu24,
        IsoMu27, IsoMu21MediumIsoTau32, TrigPass, mPt, mEta, mPhi,
        tPt, tEta, tPhi, tMVAIsoLoose, tMVAIsoMedium, tMVAIsoTight,
        tMVAIsoVTight, m_vis, mt, SS;
      //std::map< std::string, float > varMap;
      //std::string run_s = "run";
      //std::string lumi_s = "lumi";
      //varMap[run_s] = float(run);
      //varMap[lumi_s] = float(lumi);
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
    triggerToken2_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerSrc2")))
{
   //now do what ever initialization is needed
   //usesResource("TFileService");
   edm::Service<TFileService> fs;
   TFileDirectory subDir = fs->mkdir( "tagAndProbe" );
   nEvents = subDir.make<TH1D>("nEvents","nEvents",1,-0.5,0.5);
   tree = subDir.make<TTree>("Ntuple","My T-A-P Ntuple");
   tree->Branch("run",&run,"run/F");
   tree->Branch("lumi",&lumi,"lumi/F");
   tree->Branch("eventD",&eventD,"eventD/D");
   //for (auto& pair : varMap ) {
   //   tree->Branch( pair.first,&pair.second,pair.first+"/F");
   //}

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
    edm::Handle<std::vector<reco::GenJet>> hTaus;   
    iEvent.getByToken(genHadronicTausToken_, hTaus);
    edm::Handle<std::vector<reco::GenJet>> eTaus;   
    iEvent.getByToken(genElectronicTausToken_, eTaus);
    edm::Handle<std::vector<reco::GenJet>> mTaus;   
    iEvent.getByToken(genMuonicTausToken_, mTaus);
    edm::Handle<std::vector<PileupSummaryInfo>> puInfo;   
    iEvent.getByToken(puToken_, puInfo);
    edm::Handle<edm::TriggerResults> trigger1;   
    iEvent.getByToken(triggerToken1_, trigger1);
    edm::Handle<edm::TriggerResults> trigger2;   
    iEvent.getByToken(triggerToken2_, trigger2);

    edm::Handle<std::vector<reco::Vertex>> vertices;   
    iEvent.getByToken(vertexToken_, vertices);
    if (vertices->empty()) return; // skip the event if no PV found
    const reco::Vertex &PV = vertices->front();


    edm::Handle<std::vector<pat::Muon>> muons;   
    iEvent.getByToken(muonToken_, muons);
    for (const pat::Muon &mu : *muons) {
        if (mu.pt() < 5 || !mu.isLooseMuon()) continue;
        printf("muon with pt %4.1f, dz(PV) %+5.3f, POG loose id %d, tight id %d\n",
                mu.pt(), mu.muonBestTrack()->dz(PV.position()), mu.isLooseMuon(), mu.isTightMuon(PV));
    }


    edm::Handle<std::vector<pat::Electron>> electrons;   
    iEvent.getByToken(electronToken_, electrons);
    for (const pat::Electron &el : *electrons) {
        if (el.pt() < 5) continue;
        printf("elec with pt %4.1f, supercluster eta %+5.3f, sigmaIetaIeta %.3f (%.3f with full5x5 shower shapes), pass conv veto %d\n",
                    el.pt(), el.superCluster()->eta(), el.sigmaIetaIeta(), el.full5x5_sigmaIetaIeta(), el.passConversionVeto());
    }


    edm::Handle<std::vector<pat::Tau>> taus;   
    iEvent.getByToken(tauToken_, taus);
    for (const pat::Tau &tau : *taus) {
        if (tau.pt() < 20) continue;
        printf("tau  with pt %4.1f, dxy signif %.1f, ID(byMediumCombinedIsolationDeltaBetaCorr3Hits) %.1f, lead candidate pt %.1f, pdgId %d \n",
                    tau.pt(), tau.dxy_Sig(), tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"), tau.leadCand()->pt(), tau.leadCand()->pdgId());
    }


    edm::Handle<std::vector<pat::Jet>> jets;   
    iEvent.getByToken(jetToken_, jets);
    int ijet = 0;
    for (const pat::Jet &j : *jets) {
        if (j.pt() < 20) continue;
        printf("jet  with pt %5.1f (raw pt %5.1f), eta %+4.2f, btag CSV %.3f, CISV %.3f, pileup mva disc %+.2f\n",
            j.pt(), j.pt()*j.jecFactor("Uncorrected"), j.eta(), std::max(0.f,j.bDiscriminator("combinedSecondaryVertexBJetTags")), std::max(0.f,j.bDiscriminator("combinedInclusiveSecondaryVertexBJetTags")), j.userFloat("pileupJetId:fullDiscriminant"));
        if ((++ijet) == 1) { // for the first jet, let's print the leading constituents
            std::vector<reco::CandidatePtr> daus(j.daughterPtrVector());
            std::sort(daus.begin(), daus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) { return p1->pt() > p2->pt(); }); // the joys of C++11
            for (unsigned int i2 = 0, n = daus.size(); i2 < n && i2 <= 3; ++i2) {
                const pat::PackedCandidate &cand = dynamic_cast<const pat::PackedCandidate &>(*daus[i2]);
                printf("         constituent %3d: pt %6.2f, dz(pv) %+.3f, pdgId %+3d\n", i2,cand.pt(),cand.dz(PV.position()),cand.pdgId());
            }
        }
    }
 
    edm::Handle<std::vector<pat::MET>> mets;   
    iEvent.getByToken(metToken_, mets);
    const pat::MET &met = mets->front();
    printf("MET: pt %5.1f, phi %+4.2f, sumEt (%.1f). genMET %.1f. MET with JES up/down: %.1f/%.1f\n",
        met.pt(), met.phi(), met.sumEt(),
        met.genMET()->pt(),
        met.shiftedPt(pat::MET::JetEnUp), met.shiftedPt(pat::MET::JetEnDown));

    printf("\n");

    run = -1.0;
    lumi = -1.0;
    eventD = -1.0;
   //for (auto& pair : varMap ) {
   //   pair.second = -1.0;
   //}

    // First thing, fill the nEvents
    nEvents->Fill(0.);

    //std::cout << iEvent.eventAuxiliary().event() << std::endl;
    run = iEvent.eventAuxiliary().run();
    lumi = iEvent.eventAuxiliary().luminosityBlock();
    //varMap["run"] = iEvent.eventAuxiliary().run();
    //varMap["lumi"] = iEvent.eventAuxiliary().luminosityBlock();
    eventD = iEvent.eventAuxiliary().event();
    //std::cout << "Run: " <<varMap["run"]<< "Evt: " <<eventD<< "Lumi: " <<varMap["lumi"]<<std::endl;
    std::cout << "Run: " <<run<< "Evt: " <<eventD<< "Lumi: " <<lumi<<std::endl;

    // Get the number of true events
    // This is used later for pile up reweighting
    if (puInfo->size() > 0) {
        //std::cout<<"pu size = "<<puInfo->size()<<std::endl;
        //std::cout<<puInfo->at(1).getTrueNumInteractions()<<std::endl;
        nTruePU = puInfo->at(1).getTrueNumInteractions();
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
