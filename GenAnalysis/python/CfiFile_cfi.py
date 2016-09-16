import FWCore.ParameterSet.Config as cms

demo = cms.EDAnalyzer('AcceptanceAnalyzer',
    hadronSrc = cms.InputTag('tauGenJetsSelectorAllHadrons'),
    electronSrc = cms.InputTag('tauGenJetsSelectorElectrons'),
    muonSrc = cms.InputTag('tauGenJetsSelectorMuons'),
    puSrc = cms.InputTag('slimmedAddPileupInfo'),
    lheSrc = cms.InputTag('externalLHEProducer')
)
