import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('python')

options.register('outputFilename', 'BtagEffiecncy.root',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Output file name"
                 )

options.register('process', 'QCD',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "MC-simulated event type"
)

options.register('baseJetCollection','slimmedJets',
                 VarParsing.multiplicity.singleton,
                 VarParsing.varType.string,
                 "Base jet collection"
                 )


process = cms.Process("Demo")


process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.MagneticField_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')




process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = ''
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
fileNames = cms.untracked.vstring('/store/mc/RunIISummer16MiniAODv2/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/50000/36CDAE89-B3BE-E611-B022-0025905B8604.root'),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
    )


if options.process == "QCD":
    process.source.fileNames = [
        # /QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_magnetOn_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM
        'root://cmseos.fnal.gov//store/user/cmsdas/2017/short_exercises/BTagging/PUFlat0to70_80X_mcRun2_asymptotic_2016_TrancheIV_v4-v1/50000/00BC8956-278B-E611-99AD-0CC47A4D763C.root'
    ]


if options.process == "DATA":
    process.source.fileNames = [
        # /QCD_Pt-15to7000_TuneCUETP8M1_Flat_13TeV_pythia8/RunIIFall15MiniAODv2-PU25nsData2015v1_magnetOn_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM
        'root://cms-xrd-global.cern.ch//store/data/Run2016B/SingleElectron/MINIAOD/18Apr2017_ver2-v1/00000/00ED7FAC-3540-E711-954F-D067E5F9217B.root'
    ]

## Output file
process.TFileService = cms.Service("TFileService",
   fileName = cms.string(options.outputFilename.replace('.root','_' + options.process + '.root'))
)


pileupProductName = "slimmedAddPileupInfo"


#tfile service
#process.TFileService = cms.Service("TFileService",
#                                   fileName = cms.string(options.outFilename)
#                                   )
#################################################
## Update PAT jets
#################################################

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

## b-tag discriminators
bTagDiscriminators = [
    'pfTrackCountingHighEffBJetTags',
    'pfTrackCountingHighPurBJetTags',
    'pfJetProbabilityBJetTags',
    'pfJetBProbabilityBJetTags',
    'pfSimpleSecondaryVertexHighEffBJetTags',
    'pfSimpleSecondaryVertexHighPurBJetTags',
    'pfCombinedSecondaryVertexV2BJetTags',
    'pfCombinedInclusiveSecondaryVertexV2BJetTags',
    'pfCombinedMVAV2BJetTags'
]

from PhysicsTools.PatAlgos.tools.jetTools import *
## Update the slimmedJets in miniAOD: corrections from the chosen Global Tag are applied and the b-tag discriminators are re-evaluated
updateJetCollection(
    process,
    jetSource = cms.InputTag('slimmedJets'),
    jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    btagDiscriminators = bTagDiscriminators
)

#################################################


process.demo = cms.EDAnalyzer('BtagEfficiency',
				#jets                = cms.InputTag("slimmedJets"), 
				jets                = cms.InputTag("slimmedJetsPuppi"), 
				#jets                = cms.InputTag("slimmedJetsAK8PFCHSSoftDropPacked"), 
				#jets                = cms.InputTag("slimmedJetsAK8"), 
				vertices            = cms.InputTag("offlineSlimmedPrimaryVertices"),
				secVertices         = cms.InputTag("slimmedSecondaryVertices"),
				genEventInfoProduct = cms.InputTag('generator'),
				pileup              = cms.InputTag( pileupProductName ),
                                rho                 = cms.InputTag("fixedGridRhoFastjetAll"),
                                beamSpot            = cms.InputTag('offlineBeamSpot'),
				electrons           = cms.InputTag("slimmedElectrons"),
				generator           = cms.InputTag("slimmedGenJets"),

				trackLabel          = cms.InputTag('packedPFCandidates'), 
				tracks              = cms.InputTag("packedPFCandidates"),
		bDiscriminators = cms.vstring(      # list of b-tag discriminators to access
				         'pfTrackCountingHighEffBJetTags',
				         'pfTtrackCountingHighPurBJetTags',
				         'pfJetProbabilityBJetTags',
					 'pfJetBProbabilityBJetTags',
				         'pfSimpleSecondaryVertexHighEffBJetTags',
				         'pfSimpleSecondaryVertexHighPurBJetTags',
				         'pfCombinedSecondaryVertexV2BJetTags',
				         'pfCombinedInclusiveSecondaryVertexV2BJetTags',
					 'pfCombinedMVAV2BJetTags'
    						)
)


#process.p = cms.Path(process.demo)
process.p = cms.Path(process.TrackRefitter*process.demo)
