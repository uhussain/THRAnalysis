import FWCore.ParameterSet.Config as cms

tagAndProbe = cms.EDAnalyzer('TagAndProbeAnalyzer',
    hadronSrc = cms.InputTag('tauGenJetsSelectorAllHadrons'),
    tauElectronSrc = cms.InputTag('tauGenJetsSelectorElectrons'),
    tauMuonSrc = cms.InputTag('tauGenJetsSelectorMuons'),
    puSrc = cms.InputTag('slimmedAddPileupInfo'),
    tauSrc = cms.InputTag('slimmedTaus'),
    muonSrc = cms.InputTag('slimmedMuons'),
    electronSrc = cms.InputTag('slimmedElectrons'),
    jetSrc = cms.InputTag('slimmedJets'),
    metSrc = cms.InputTag('slimmedMETs'),
    pvSrc = cms.InputTag('offlineSlimmedPrimaryVertices'),
    triggerSrc1 = cms.InputTag("TriggerResults","","HLT"),
    triggerSrc2 = cms.InputTag("TriggerResults","","HLT2"),
)
