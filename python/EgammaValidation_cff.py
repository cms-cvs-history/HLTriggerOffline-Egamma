import FWCore.ParameterSet.Config as cms

pathsElectron = ['veryHighEtDQM',
                 'highEtDQM',
                 'singleElectronRelaxedDQM',
                 'singleElectronDQM',
                 'singleElectronRelaxedLargeWindowDQM', 
                 'singleElectronLargeWindowDQM',
                 'doubleElectronRelaxedDQM',
                 'doubleElectronDQM']

pathsPhoton = ['veryHighEtDQM',
               'highEtDQM',
               'singlePhotonRelaxedDQM',
               'singlePhotonDQM',
               'doublePhotonRelaxedDQM',
               'doublePhotonDQM']

#define common modules
leptons = cms.EDFilter("PdgIdAndStatusCandViewSelector",
    status = cms.vint32(1),
    src = cms.InputTag("genParticles"),
    pdgId = cms.vint32(11)
)
cut = cms.EDFilter("EtaPtMinCandViewSelector",
    src = cms.InputTag("leptons"),
    etaMin = cms.double(-2.5),
    etaMax = cms.double(2.5),
    ptMin = cms.double(2.0)
)

#define sequences/noncommon modules
selZ = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("cut"),
    minNumber = cms.uint32(2)
)
Zseq='*('

selW = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("cut"),
    minNumber = cms.uint32(1)
)
Wseq='*('

selPJ = cms.EDFilter("CandViewCountFilter",
    src = cms.InputTag("cutPhoton"),
    minNumber = cms.uint32(1)
)
PJseq='*('

###########################################################
#  Electron DQM

first= True
#load modules
for trig in pathsElectron:
    if not first:
        Zseq=Zseq+'+'
        Wseq=Wseq+'+'
    first= False    

    imp = 'from HLTriggerOffline.Egamma.' + trig + '_cfi import *'
    exec imp

    #clone for Z
    clon = trig + '_Z = ' + trig + '.clone()'
    exec clon
    #adjust MC match pid
    mcmatch = trig + '_Z.pdgGen=11'
    exec mcmatch
    Zseq=Zseq + trig + '_Z' 

    #clone for W
    clon = trig + '_W = ' + trig + '.clone()'
    exec clon
    #adjust MC match pid
    mcmatch = trig + '_W.pdgGen=11'
    exec mcmatch
    Wseq=Wseq + trig + '_W'

Zseq=Zseq + ')'
Wseq=Wseq + ')'
###########################################################


###########################################################
#  Photon DQM

first= True
#load modules
for trig in pathsPhoton:
    if not first:
       PJseq=PJseq+'+'
    first= False
    
    imp = 'from HLTriggerOffline.Egamma.' + trig + '_cfi import *'
    exec imp

    #clone for Photon+Jet
    clon = trig + '_PJ = ' + trig + '.clone()'
    exec clon
    #adjust MC match pid
    mcmatch = trig + '_PJ.pdgGen=22'
    exec mcmatch
    PJseq=PJseq + trig + '_PJ'

PJseq=PJseq + ')'
###########################################################


###########################################################
# Electron DQM
#scom = 'egammavalZee = cms.Sequence(leptons*cut*selZ' + Zseq +')'
scom = 'egammavalZee = cms.Sequence(leptons*cut' + Zseq +')'
exec scom
#scom = 'egammavalWenu = cms.Sequence(leptons*cut*selW' + Wseq +')'
scom = 'egammavalWenu = cms.Sequence(leptons*cut' + Wseq +')'
exec scom

# Photon DQM
leptons.pdgId=cms.vint32(22)
#scom = 'egammavalPhotonJet = cms.Sequence(leptons*cut*selPJ' + PJseq +')'
scom = 'egammavalPhotonJet = cms.Sequence(leptons*cut' + PJseq +')'
exec scom
