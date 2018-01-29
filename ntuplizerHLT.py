#!/usr/bin/python
import ROOT
import itertools
import resource
from array import array
from math import sqrt, pi, log10, log, exp
# load FWlite python libraries
from DataFormats.FWLite import Handle, Events
from utils import deltaR,SetVariable,DummyClass,productWithCheck,checkTriggerIndex

#ROOT.gROOT.LoadMacro("/scratch/sdonato/NtupleForPaolo/CMSSW_8_0_3_patch1/src/DataFormats/L1Trigger/interface/EtSumHelper.h")

Handle.productWithCheck = productWithCheck

maxJets         = 50
bunchCrossing   = 0
pt_min          = 20

def FillVector(source,variables,minPt=pt_min):
    variables.num[0] = 0
    for obj in source.productWithCheck():
        if obj.pt()<minPt: continue
        if variables.num[0]<len(variables.pt):
            for (name,var) in variables.__dict__.items():
                if name == "pt" :           var[variables.num[0]] = obj.pt()
                elif name == "eta" :        var[variables.num[0]] = obj.eta()
                elif name == "phi" :        var[variables.num[0]] = obj.phi()
                elif name == "mass" :       var[variables.num[0]] = obj.mass()
                elif name == "softdropmass":var[variables.num[0]] = obj.userFloat('ak8PFJetsCHSValueMap:ak8PFJetsCHSSoftDropMass')
                elif name == "neHadEF" :    var[variables.num[0]] = obj.neutralHadronEnergyFraction()
                elif name == "neEmEF" :     var[variables.num[0]] = obj.neutralEmEnergyFraction()
                elif name == "chHadEF" :    var[variables.num[0]] = obj.chargedHadronEnergyFraction()
                elif name == "chEmEF" :     var[variables.num[0]] = obj.chargedEmEnergyFraction()
                elif name == "muEF" :       var[variables.num[0]] = obj.muonEnergyFraction()
                elif name == "mult" :       var[variables.num[0]] = obj.chargedMultiplicity()+obj.neutralMultiplicity();
                elif name == "neMult" :     var[variables.num[0]] = obj.neutralMultiplicity()
                elif name == "chMult" :     var[variables.num[0]] = obj.chargedMultiplicity()
                elif name == "doubleBTag" and getattr(obj, 'bDiscriminator', 0):
                                            var[variables.num[0]] = obj.bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags")
            variables.num[0] += 1


def FillBtag(btags_source, jets, jet_btags):
    for i in range(jets.num[0]):
        jet_btags[i] = -20.
        dRmax = 0.3
        btags = btags_source.productWithCheck()
        for j in range(0,btags.size()):
            jetB = btags.key(j).get()
            dR = deltaR(jetB.eta(),jetB.phi(),jets.eta[i],jets.phi[i])
            if dR<dRmax:
                jet_btags[i] = max(-1.,btags.value(j))
                dRmax = dR

def Matching(phi, eta, jets):
    index = -1
    for i in range(jets.num[0]):
        dRmax = 0.3
        dR = deltaR(eta,phi,jets.eta[i],jets.phi[i])
        if dR<dRmax:
            index = i
            dRmax = dR
    return index


def getVertex(vertex_source):
    vertices = vertex_source.productWithCheck()
    if vertices.size()>0:
        return vertices.at(0).z()
    else:
        return -1000

def WithFallback(product,method="pt"):
    if product.size()>0:
        return getattr(product[0],method)()
    else:
        return -10

def BookVector(tree,name="vector",listMembers=[]):
    obj = DummyClass()
    obj.num   = SetVariable(tree,name+'_num' ,'I')
    for member in listMembers:
        if "match" in name:
            setattr(obj,member,SetVariable(tree,name+'_'+member  ,'I',name+'_num',maxJets))
        else:
            setattr(obj,member,SetVariable(tree,name+'_'+member  ,'F',name+'_num',maxJets))
    return obj

    ##########################################################################

def launchNtupleFromHLT(fileOutput,filesInput, secondaryFiles, maxEvents,preProcessing=True, firstEvent=0):
    bunchCrossing   = 12
    print "filesInput: ",filesInput
    print "fileOutput: ",fileOutput
    print "secondaryFiles: ",secondaryFiles
    print "maxEvents: ",maxEvents
    print "preProcessing: ",preProcessing
    print "firstEvent: ",firstEvent

    isMC = True  # 'SIM' in filesInput[0]

    ## Pre-processing
    if preProcessing:
        from PhysicsTools.Heppy.utils.cmsswPreprocessor import CmsswPreprocessor
        from PhysicsTools.HeppyCore.framework.config import MCComponent
        cmsRun_config = "HLT_tunedcontent_dump.py"
        preprocessor = CmsswPreprocessor(cmsRun_config)
        cfg = MCComponent("OutputHLT",filesInput, secondaryfiles=secondaryFiles)
        print "Run cmsswPreProcessing using:"
        print cfg.name
        print cfg.files
        print cfg.secondaryfiles
        print
        try:
            preprocessor.run(cfg,".",firstEvent,maxEvents)
        except:
            print "cmsswPreProcessing failed!"
            print "cat cmsRun_config.py"
            config = file(cmsRun_config)
            print config.read()
            print "cat cmsRun.log"
            log = file("cmsRun.log")
            print log.read()
            preprocessor.run(cfg,".",firstEvent,maxEvents)
            raise Exception("CMSSW preprocessor failed!")

    f = ROOT.TFile(fileOutput,"recreate")
    tree = ROOT.TTree("tree","tree")

    if preProcessing:
        fwLiteInputs = ["cmsswPreProcessing.root"]
    else:
        fwLiteInputs = filesInput
    if not filesInput:
        exit(-1)
    import os.path
    if not os.path.isfile(fwLiteInputs[0]):
        raise Exception( fwLiteInputs[0] + " does not exist.")
    events = Events (fwLiteInputs)

    ### list of input variables ###

    # triggerBits, triggerBitLabel = Handle("edm::TriggerResults"), ("TriggerResults::MYHLT")
    triggerBits, triggerBitLabel                    = Handle("edm::TriggerResults"), ("TriggerResults")  # ("TriggerResults::TEST20180119174441")
    pileUp_source, pileUp_label                     = Handle("vector<PileupSummaryInfo>"), ("addPileupInfo")
    generator_source, generator_label               = Handle("GenEventInfoProduct"), ("generator")
    genParticles_source, genParticles_label         = Handle("vector<reco::GenParticle>"), ("genParticles")

    # NOPE, not there yet  /  ak8-ify!!!!    caloJets_source, caloJets_label    = Handle("vector<reco::CaloJet>"), ("hltAK4CaloJetsCorrectedIDPassed")
    # NOPE, not there yet  /  ak8-ify!!!!    calobtag_source, calobtag_label    = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>"), ("hltCombinedSecondaryVertexBJetTagsCalo")

    PixelVertices_source, PixelVertices_label       = Handle("vector<reco::Vertex>"), ("hltPixelVertices")
    VerticesPF_source, VerticesPF_label             = Handle("vector<reco::Vertex>"), ("hltVerticesPF")
    VerticesL3_source, VerticesL3_label             = Handle("vector<reco::Vertex>"), ("hltVerticesL3")

    hltPFJetsAK8_source, hltPFJetsAK8_label         = Handle("vector<reco::PFJet>"),  ('hltPFJetForBtagAK8')  #  ('hltPFJetForBtagAK8','','TEST20180119174441')
    hltbtags_source, hltbtag_label                  = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>"), ("hltBoostedDBJetTagsPFAK8")

    offlinePfJets_source, offlinePfJets_label       = Handle("vector<pat::Jet>"), ("slimmedJetsAK8")
    #offlinePfJets_source, offlinePfJets_label       = Handle("vector<reco::Jet>"), ("ak8PFJetsCHS")
    #hltbtags_source, hltbtag_label                  = Handle("edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>"), ("boostedDoubleSecondaryVertexBJetTagsPFAK8???")

    genJetsNoNu_source, genJetsNoNu_label           = Handle("vector<reco::GenJet>"), ("ak8GenJetsNoNu")
    genJets_source, genJets_label                   = Handle("vector<reco::GenJet>"), ("ak8GenJets")


    ### create output variables ###

    evt         = SetVariable(tree,'evt')
    lumi        = SetVariable(tree,'lumi')
    run         = SetVariable(tree,'run')

    hltpfJetsAK8    = BookVector(tree,"hltpfJetsAK8",['pt','eta','phi','mass','matchOff','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','doubleBTag'])
    pfJetsAK8       = BookVector(tree,"pfJetsAK8",['pt','eta','phi','mass','softdropmass','matchGen','neHadEF','neEmEF','chHadEF','chEmEF','muEF','mult','neMult','chMult','doubleBTag'])

    if isMC:
        pu                  = SetVariable(tree,'pu')
        ptHat               = SetVariable(tree,'ptHat')
        maxPUptHat          = SetVariable(tree,'maxPUptHat')
        genJets             = BookVector(tree,"genJets",['pt','eta','phi','mass','mcFlavour','mcPt'])
        genJetsNoNu         = BookVector(tree,"genJets",['pt','eta','phi','mass','mcFlavour','mcPt'])
        PixelVertices       = SetVariable(tree,'PixelVertices')
        VerticesPF          = SetVariable(tree,'VerticesPF')
        VerticesL3          = SetVariable(tree,'VerticesL3')
        trueVertex          = SetVariable(tree,'trueVertex')

    f.cd()

    ##get trigger names
    events.to(0)
    try:
        event = next(iter(events))
    except StopIteration:
        f.Write()
        f.Close()
        return

    event.getByLabel(triggerBitLabel, triggerBits)
    names = event.object().triggerNames(triggerBits.product())
    triggerNames = names.triggerNames()
    for name in triggerNames:
        name = name.split("_v")[0]
    nTriggers = len(triggerNames)
    triggerVars = {}
    for trigger in triggerNames:
        triggerVars[trigger]=array( 'i', [ 0 ] )
        tree.Branch( trigger, triggerVars[trigger], trigger+'/O' )

    ##event loop
    events.to(0)
    for iev,event in enumerate(events):
        if iev>maxEvents and maxEvents>=0:
            break

        run[0]  = event.eventAuxiliary().run()
        lumi[0] = event.eventAuxiliary().luminosityBlock()
        evt[0]  = event.eventAuxiliary().event()

        event.getByLabel(triggerBitLabel, triggerBits)
        event.getByLabel(offlinePfJets_label, offlinePfJets_source)
        event.getByLabel(hltPFJetsAK8_label, hltPFJetsAK8_source)
        event.getByLabel(hltbtag_label, hltbtags_source)
        if isMC:
            event.getByLabel(generator_label, generator_source)
            event.getByLabel(pileUp_label, pileUp_source)
            event.getByLabel(genJetsNoNu_label, genJetsNoNu_source)
            event.getByLabel(genJets_label, genJets_source)
            event.getByLabel(PixelVertices_label, PixelVertices_source)
            event.getByLabel(VerticesPF_label, VerticesPF_source)
            event.getByLabel(VerticesL3_label, VerticesL3_source)
            event.getByLabel(genParticles_label, genParticles_source)

        ####################################################

        FillVector(hltPFJetsAK8_source, hltpfJetsAK8)
        FillBtag(hltbtags_source, hltpfJetsAK8, hltpfJetsAK8.doubleBTag)
        FillVector(offlinePfJets_source, pfJetsAK8, 15)
        FillVector(genJets_source,genJets,15)
        FillVector(genJetsNoNu_source,genJetsNoNu,15)

        if isMC:
            PixelVertices[0] = getVertex(PixelVertices_source)
            VerticesPF[0] = getVertex(VerticesPF_source)
            VerticesL3[0] = getVertex(VerticesL3_source)
            trueVertex[0] = genParticles_source.productWithCheck().at(2).vertex().z()


            for i in range(hltpfJetsAK8.num[0]):
                hltpfJetsAK8.matchOff[i] = Matching(hltpfJetsAK8.phi[i],hltpfJetsAK8.eta[i],pfJetsAK8)
                hltpfJetsAK8.matchGen[i] = Matching(hltpfJetsAK8.phi[i],hltpfJetsAK8.eta[i],genJetsNoNu)

            for i in range(pfJetsAK8.num[0]):
                pfJetsAK8.matchGen[i] = Matching(pfJetsAK8.phi[i],pfJetsAK8.eta[i],genJetsNoNu)

            for i in range(genJetsNoNu.num[0]):
                genJetsNoNu.mcFlavour[i] = -100
                genJetsNoNu.mcPt[i] = -100

            # for genParticle in genParticles_source.productWithCheck():
            #     if genParticle.pt()<5:
            #         continue
            #     if not (abs(genParticle.pdgId()) in [21,1,2,3,4,5,11,13,15]):
            #         continue
            #     if genParticle.numberOfMothers() and genParticle.mother().pt()>5 and (abs(genParticle.mother().pdgId()) in [21,1,2,3,4,5,11,13]):
            #         continue
            #     if evt[0]==7826939:
            #         print "genParticle:"
            #         print genParticle.pt(),genParticle.eta(),genParticle.phi(),genParticle.pdgId()
            #         print "genJetsNoNu:"
            #     for i in range(genJetsNoNu.num[0]):
            #         if genParticle.pt()<0.2*genJetsNoNu.pt[i]:
            #             continue
            #         if deltaR(genParticle.eta(),genParticle.phi(),genJetsNoNu.eta[i],genJetsNoNu.phi[i])<0.4:
            #             if evt[0]==7826939:
            #                 print genJetsNoNu.pt[i],genJetsNoNu.eta[i],genJetsNoNu.eta[i],genJetsNoNu.mcFlavour[i]
            #                 print "not (int(abs(genJetsNoNu.mcFlavour[i])) in [5,4,3,2,1]):",not (int(abs(genJetsNoNu.mcFlavour[i])) in [5,4,3,2,1])
            #             if abs(genParticle.pdgId())==5:
            #                 if genJetsNoNu.mcFlavour[i]!=5 or genParticle.pt()>genJetsNoNu.mcPt[i]:
            #                     genJetsNoNu.mcFlavour[i] = genParticle.pdgId()
            #                     genJetsNoNu.mcPt[i]      = genParticle.pt()
            #             elif abs(genParticle.pdgId())==4 and not int(abs(genJetsNoNu.mcFlavour[i])) in [5]:
            #                 if genJetsNoNu.mcFlavour[i]!=4 or genParticle.pt()>genJetsNoNu.mcPt[i]:
            #                     genJetsNoNu.mcFlavour[i] = genParticle.pdgId()
            #                     genJetsNoNu.mcPt[i]      = genParticle.pt()
            #             elif abs(genParticle.pdgId())==3 and not int(abs(genJetsNoNu.mcFlavour[i])) in [5,4]:
            #                 if genJetsNoNu.mcFlavour[i]!=3 or genParticle.pt()>genJetsNoNu.mcPt[i]:
            #                     genJetsNoNu.mcFlavour[i] = genParticle.pdgId()
            #                     genJetsNoNu.mcPt[i]      = genParticle.pt()
            #             elif abs(genParticle.pdgId())==2 and not int(abs(genJetsNoNu.mcFlavour[i])) in [5,4,3]:
            #                 if genJetsNoNu.mcFlavour[i]!=2 or genParticle.pt()>genJetsNoNu.mcPt[i]:
            #                     genJetsNoNu.mcFlavour[i] = genParticle.pdgId()
            #                     genJetsNoNu.mcPt[i]      = genParticle.pt()
            #             elif abs(genParticle.pdgId())==1 and not int(abs(genJetsNoNu.mcFlavour[i])) in [5,4,3,2]:
            #                 if genJetsNoNu.mcFlavour[i]!=1 or genParticle.pt()>genJetsNoNu.mcPt[i]:
            #                     genJetsNoNu.mcFlavour[i] = genParticle.pdgId()
            #                     genJetsNoNu.mcPt[i]      = genParticle.pt()
            #             elif abs(genParticle.pdgId())==21 and not (int(abs(genJetsNoNu.mcFlavour[i])) in [5,4,3,2,1]):
            #                 if genJetsNoNu.mcFlavour[i]!=21 or genParticle.pt()>genJetsNoNu.mcPt[i]:
            #                     genJetsNoNu.mcFlavour[i] = genParticle.pdgId()
            #                     genJetsNoNu.mcPt[i]      = genParticle.pt()
            #             elif abs(genParticle.pdgId()) in [11,13] and not int(abs(genJetsNoNu.mcFlavour[i])) in [5,4,3,2,1,21]:
            #                 if not (genJetsNoNu.mcFlavour[i] in [11,13]) or genParticle.pt()>genJetsNoNu.mcPt[i]:
            #                     genJetsNoNu.mcFlavour[i] = genParticle.pdgId()
            #                     genJetsNoNu.mcPt[i]      = genParticle.pt()
            #             elif abs(genParticle.pdgId()) in [15] and not int(abs(genJetsNoNu.mcFlavour[i])) in [5,4,3,2,1,21,11,13]:
            #                 if not (genJetsNoNu.mcFlavour[i] in [15]) or genParticle.pt()>genJetsNoNu.mcPt[i]:
            #                     genJetsNoNu.mcFlavour[i] = genParticle.pdgId()
            #                     genJetsNoNu.mcPt[i]      = genParticle.pt()
            #             elif abs(genParticle.pdgId()) in [22] and not int(abs(genJetsNoNu.mcFlavour[i])) in [5,4,3,2,1,21,11,13,15]:
            #                 if not (genJetsNoNu.mcFlavour[i] in [22]) or genParticle.pt()>genJetsNoNu.mcPt[i]:
            #                     genJetsNoNu.mcFlavour[i] = genParticle.pdgId()
            #                     genJetsNoNu.mcPt[i]      = genParticle.pt()
            #             if evt[0]==7826939:
            #                 print "newFlav:",genJetsNoNu.mcFlavour[i]

            if bunchCrossing>=pileUp_source.productWithCheck().size() or pileUp_source.productWithCheck().at(bunchCrossing).getBunchCrossing()!=0:
                found=False
                for bunchCrossing in range(pileUp_source.productWithCheck().size()):
                    if pileUp_source.productWithCheck().at(bunchCrossing).getBunchCrossing() == 0 :
                        found=True
                        break
                if not found:
                    Exception("Check pileupSummaryInfos!")
                print "I'm using bunchCrossing=",bunchCrossing
            pu[0] = pileUp_source.productWithCheck().at(bunchCrossing).getTrueNumInteractions()
            ptHat[0]    = generator_source.product().qScale()

            maxPUptHat[0] = -1
            for ptHat_ in pileUp_source.productWithCheck().at(bunchCrossing).getPU_pT_hats():
                maxPUptHat[0] = max(maxPUptHat[0],ptHat_)

        # end if isMC:

        names = event.object().triggerNames(triggerBits.product())
        for i,triggerName in enumerate(triggerNames):
            index = names.triggerIndex(triggerName)
#            print "index=",index
            if checkTriggerIndex(triggerName,index,names.triggerNames()):
                triggerVars[triggerName][0] = triggerBits.product().accept(index)
#                print "acc:",triggerBits.product().accept(index)
            else:
                triggerVars[triggerName][0] = 0

        if iev%10==1:
            print "Event: ",iev," done."
        tree.Fill()

    f.Write()
    f.Close()

if __name__ == "__main__":
    #filesInput = ["root://eoscms.cern.ch//eos/cms/store/mc/PhaseIFall16DR/GluGluToRSGravitonToHHTo4B_M-450_narrow_13TeV-madgraph/GEN-SIM-RAW/FlatPU28to62HcalNZSRAW_90X_upgrade2017_realistic_v6_C1-v1/80000/F8161EEB-9810-E711-A85C-FA163E0B564E.root"]
    filesInput = ["root://eoscms.cern.ch//eos/cms/store/mc/PhaseIFall16DR/GluGluToRSGravitonToHHTo4B_M-450_narrow_13TeV-madgraph/GEN-SIM-RAW/FlatPU28to62HcalNZSRAW_90X_upgrade2017_realistic_v6_C1-v1/80000/EAACC9D6-0F11-E711-A9B1-FA163EDAFEAB.root"]
    filesInput = ["root://eoscms.cern.ch//eos/cms/store/data/Run2017C/HLTPhysics8/RAW/v1/000/301/959/00000/70AA2AD0-C08B-E711-9723-02163E011C82.root"]
#    filesInput = ["file:/mnt/t3nfs01/data01/shome/sdonato/QCD470_GEN-SIM-RAW_PhaseI_83X_FlatPU28to62.root"]
    secondaryFiles = []
    fileOutput = "tree.root"
    maxEvents = 10000000

    filesInput = ["/pnfs/desy.de/cms/tier2/store/user/htholen/HLT_EDMTuple_DoubleBTag_Hbb_Signal_v0p5/GluGluHToBB_M125_13TeV_powheg_pythia8/HLT_EDMTuple_DoubleBTag_Hbb_Signal_v0p5_GluGluHToBB_M125_13TeV_powheg_pythia8/180126_120400/0000/HLT_tunedcontent_v2-1_79.root"]
    filesInput = ["edm_tuple.root"]

    launchNtupleFromHLT(fileOutput,filesInput,secondaryFiles,maxEvents, preProcessing=False)
