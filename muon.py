import FWCore.ParameterSet.Config as cms

process = cms.Process("muonexttocvs")



process.load("FWCore.MessageService.MessageLogger_cfi")



process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )



process.source = cms.Source("PoolSource",

    fileNames = cms.untracked.vstring(
'root://eospublic.cern.ch//eos/opendata/cms/Run2011A/DoubleMu/AOD/12Oct2013-v1/10000/000D143E-9535-E311-B88B-002618943934.root',

        'root://eospublic.cern.ch//eos/opendata/cms/Run2011A/ElectronHad/AOD/12Oct2013-v1/20001/001F9231-F141-E311-8F76-003048F00942.root'

    )

)


#needed to get the actual prescale values used from the global tag

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.GlobalTag.connect = cms.string('sqlite_file:/cvmfs/cms-opendata-conddb.cern.ch/FT_53_LV5_AN1_RUNA.db')

process.GlobalTag.globaltag = 'FT_53_LV5_AN1::All'

process.muonextractorToCsv = cms.EDAnalyzer('MuonObjectInfoExtractorToCsv',
    processName = cms.string("HLT"),
                              triggerName = cms.string("@"),         
                              datasetName = cms.string("Mu"),           
                              triggerResults = cms.InputTag("TriggerResults","","HLT"),
                              triggerEvent   = cms.InputTag("hltTriggerSummaryAOD","","HLT"), 
InputCollection = cms.InputTag("muons"),

maxNumberMuons = cms.untracked.int32(3)#default is 5

)





process.p = cms.Path(process.muonextractorToCsv)
