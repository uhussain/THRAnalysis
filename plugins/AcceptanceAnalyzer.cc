// -*- C++ -*-
//
// Package:    Acceptance/AcceptanceAnalyzer
// Class:      AcceptanceAnalyzer
// 
/**\class AcceptanceAnalyzer AcceptanceAnalyzer.cc Acceptance/AcceptanceAnalyzer/plugins/AcceptanceAnalyzer.cc

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
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"



//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class AcceptanceAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit AcceptanceAnalyzer(const edm::ParameterSet&);
      ~AcceptanceAnalyzer();

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
      edm::EDGetTokenT<LHEEventProduct> lheToken_;
      TTree *tree;
      float genMass, ETauPass, MuTauPass, EMuPass, TauTauPass;
      float ETauD, MuTauD, EMuD, TauTauD;
      float threeLeptons, nLooseTaus, nLooseElec, nLooseMu;
      float nTruePU;
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
AcceptanceAnalyzer::AcceptanceAnalyzer(const edm::ParameterSet& iConfig) :
    genHadronicTausToken_(consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("hadronSrc"))),
    genElectronicTausToken_(consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("electronSrc"))),
    genMuonicTausToken_(consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("muonSrc"))),
    puToken_(consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("puSrc"))),
    lheToken_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheSrc")))
{
   //now do what ever initialization is needed
   //usesResource("TFileService");
   edm::Service<TFileService> fs;
   TFileDirectory subDir = fs->mkdir( "events" );
   tree = subDir.make<TTree>("Ntuple","My Analyzer Ntuple");
   tree->Branch("genMass",&genMass,"genMass/F");
   tree->Branch("ETauPass",&ETauPass,"ETauPass/F");
   tree->Branch("ETauD",&ETauD,"ETauD/F");
   tree->Branch("MuTauPass",&MuTauPass,"MuTauPass/F");
   tree->Branch("MuTauD",&MuTauD,"MuTauD/F");
   tree->Branch("EMuPass",&EMuPass,"EMuPass/F");
   tree->Branch("EMuD",&EMuD,"EMuD/F");
   tree->Branch("TauTauPass",&TauTauPass,"TauTauPass/F");
   tree->Branch("TauTauD",&TauTauD,"TauTauD/F");
   tree->Branch("threeLeptons",&threeLeptons,"threeLeptons/F");
   tree->Branch("nLooseTaus",&nLooseTaus,"nLooseTaus/F");
   tree->Branch("nLooseElec",&nLooseElec,"nLooseElec/F");
   tree->Branch("nLooseMu",&nLooseMu,"nLooseMu/F");
   tree->Branch("nTruePU",&nTruePU,"nTruePU/F");

}


AcceptanceAnalyzer::~AcceptanceAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
AcceptanceAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
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
    edm::Handle<LHEEventProduct> lheProd;   
    iEvent.getByToken(lheToken_, lheProd);


    genMass = -1.0;
    ETauPass = 0;
    MuTauPass = 0;
    EMuPass = 0;
    TauTauPass = 0;
    ETauD = 0;
    MuTauD = 0;
    EMuD = 0;
    TauTauD = 0;
    threeLeptons = 0;
    nLooseTaus = 0;
    nLooseElec = 0;
    nLooseMu = 0;
    nTruePU = -1.0;

    // Get the number of true events
    // This is used later for pile up reweighting
    if (puInfo->size() > 0) {
        //std::cout<<"pu size = "<<puInfo->size()<<std::endl;
        //std::cout<<puInfo->at(1).getTrueNumInteractions()<<std::endl;
        nTruePU = puInfo->at(1).getTrueNumInteractions();
    }

    // Denominator section first, just check if there's the correct #
    // of the associated leptons
    if (hTaus->size() >= 2) TauTauD = 1;
    if (hTaus->size() >= 1) {
        if(eTaus->size() >= 1) ETauD = 1;
        if(mTaus->size() >= 1) MuTauD = 1;
    }
    if(eTaus->size() >= 1 && mTaus->size() >= 1) EMuD = 1;
    
    float nLeptons = 0;
    nLeptons += hTaus->size() + eTaus->size() + mTaus->size();
    if (nLeptons >= 3) threeLeptons = 1;
    nLooseTaus = hTaus->size();
    nLooseElec = eTaus->size();
    nLooseMu = mTaus->size();


    size_t ETau_e = 0;
    size_t ETau_t = 0;
    size_t MuTau_m = 0;
    size_t MuTau_t = 0;
    size_t EMu_e = 0;
    size_t EMu_m = 0;
    size_t TauTau_t = 0;

    //const lhef::HEPEUP lhe = lheProd.product()->hepeup();
    //uint32_t nHTaus = hTaus->size();
    //if (hTaus->size() > 0)// std::cout << " --- N Hadronic Taus: "<<nHTaus<<std::endl;
    //{
    //const std::vector<reco::GenJet> nTaus_ = hTaus.product();
    for (const reco::GenJet &tau : *hTaus) {
        if ( TMath::Abs(tau.eta()) < 2.3 && tau.pt() > 20 ) {
            MuTau_t += 1;
            ETau_t += 1;}
        if ( TMath::Abs(tau.eta()) < 2.1 && tau.pt() > 40 ) TauTau_t += 1;
    }
    //}
    //uint32_t nETaus = eTaus->size();
    //if (nETaus > 0)// std::cout << " ### N Electronic Taus: "<<nETaus<<std::endl;
    //{
    for (const reco::GenJet &ele : *eTaus) {
        if ( TMath::Abs(ele.eta()) < 2.5 && ele.pt() > 13 ) EMu_e += 1;
        if ( TMath::Abs(ele.eta()) < 2.1 && ele.pt() > 24 ) ETau_e += 1;
    //}
    }
    //uint32_t nMTaus = mTaus->size();
    //if (nMTaus > 0)// std::cout << " *** N Muonic Taus: "<<nMTaus<<std::endl;
    //{
    for (const reco::GenJet &mu : *mTaus) {
        if ( TMath::Abs(mu.eta()) < 2.4 && mu.pt() > 10 ) EMu_m += 1;
        if ( TMath::Abs(mu.eta()) < 2.1 && mu.pt() > 19 ) MuTau_m += 1;
    }
    //}
    
    // Check if we have matches
    if (ETau_e > 0 && ETau_t > 0) ETauPass = 1;
    if (MuTau_m > 0 && MuTau_t > 0) MuTauPass = 1;
    if (EMu_e > 0 && EMu_m > 0) EMuPass = 1;
    if (TauTau_t > 1) TauTauPass = 1;
    

    //std::cout << lheProd << std::endl;
    //lhe
    //std::cout << lheProd.isValid() << std::endl;
    //lhef::HEPEUP lhe;
    //if (lheProd.isValid()) {
    //  lhe = lheProd->hepeup();
    //}
    const lhef::HEPEUP lhe = lheProd.product()->hepeup();
    //std::cout << lhe.ISTUP[0] << std::endl;
    std::vector<int> outgoing;
    std::vector<TLorentzVector> invmass;
    for (uint32_t i = 0; i < lhe.ISTUP.size(); ++i) {
        if (lhe.ISTUP[i]) {
            int Id = TMath::Abs( lhe.IDUP[i] );
            if (Id==21||Id==1||Id==2||Id==3||Id==4||Id==5)
                outgoing.push_back( TMath::Abs( Id ));
            int Id2 = TMath::Abs( lhe.IDUP[i] );
            if (Id2==11||Id2==13||Id2==15) {
                TLorentzVector l = TLorentzVector( lhe.PUP[i][0],
                                                   lhe.PUP[i][1],
                                                   lhe.PUP[i][2],
                                                   lhe.PUP[i][3]);
                invmass.push_back( l );
            }
        }
    }
    if (invmass.size() == 2) {
        //std::cout << "Len InvMass: "<<invmass.size()<<std::endl;
        TLorentzVector diLep = invmass[0];
        diLep += invmass[1];
        //std::cout << "m(ll) " << diLep.M() << std::endl;
        genMass = diLep.M();

    }
    
 
    
    

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
AcceptanceAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
AcceptanceAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
AcceptanceAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(AcceptanceAnalyzer);
