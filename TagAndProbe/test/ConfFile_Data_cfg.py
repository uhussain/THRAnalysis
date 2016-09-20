import FWCore.ParameterSet.Config as cms


process = cms.Process("TagAndProbe")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(5000)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring($inputFileNames)
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
                                       fileName = cms.string("$outputFileName")
                                   )



process.p = cms.Path(
            process.hltFilter*
            process.tagAndProbe)

#print process.dumpPython()



