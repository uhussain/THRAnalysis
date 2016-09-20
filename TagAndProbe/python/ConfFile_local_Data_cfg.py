import FWCore.ParameterSet.Config as cms

process = cms.Process("TagAndProbe")

process.load("FWCore.MessageService.MessageLogger_cfi")

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(1000)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'root://eoscms//eos/cms/store/data/Run2016C/SingleMuon/MINIAOD/PromptReco-v2/000/275/657/00000/6A8D3886-763B-E611-8C1F-02163E01261C.root',
    )
)



import HLTrigger.HLTfilters.triggerResultsFilter_cfi as hlt
process.hltFilter = hlt.triggerResultsFilter.clone( 
        hltResults = cms.InputTag( "TriggerResults","","HLT"), 
        triggerConditions = ( 'HLT_IsoMu20_v*', 'HLT_IsoMu22_v*',\
            'HLT_IsoMu24_v*', 'HLT_IsoMu27_v*'),
        l1tResults = '', 
        throw = False 
        )



process.load("THRAnalysis.TagAndProbe.CfiFile_cfi")




process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string('ttree.root')
                                   )



process.p = cms.Path(
            process.hltFilter*
            process.tagAndProbe)

#process.e = cms.EndPath(process.out)

#print process.dumpPython()



