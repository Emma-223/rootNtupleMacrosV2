#!/usr/bin/env python

from plot_class import *
from ROOT import *

File_preselection = GetFile(os.environ['LQDATA'] + "//eejj_analysis/eejj//scaled_output_cutTable_lq_eejj/analysisClass_lq_eejj_plots.root"     )
File_TTBar        = GetFile(os.environ['LQDATA'] + "//eejj_analysis/eejj-ttbar//output_cutTable_lq_eejj/analysisClass_lq_eejj_TTBar_plots.root")
File_QCD          = GetFile(os.environ['LQDATA'] + "//eejj_analysis/eejj_qcd//output_cutTable_lq_eejj/analysisClass_lq_eejj_QCD_plots.root"    )


#### Common values for plots:
zUncBand="no"
makeRatio=1
makeNSigma=1

pt_rebin = 2

histoBaseName = "histo1D__SAMPLE__cutHisto_allOtherCuts___________VARIABLE"
histoBaseName_userDef = "histo1D__SAMPLE__VARIABLE"

# samplesForStackHistos_other = [ "PhotonJets", "WJet_Madgraph", "SingleTop" ]
# samplesForStackHistosQCD     = ["DATA"]
# samplesForStackHistos_ZJets  = [ "TTbar_Madgraph", "ZJet_Madgraph" ]
# samplesForStackHistos = samplesForStackHistos_other + samplesForStackHistos_ZJets

samplesForStackHistos_other = [ "PhotonJets", "WJet_Madgraph", "SingleTop" ]
samplesForStackHistosQCD     = ["DATA"]
samplesForStackHistosTTBar   = ["DATA"]
samplesForStackHistos_ZJets  = [ "ZJet_Madgraph" ]
samplesForStackHistos = samplesForStackHistos_other + samplesForStackHistos_ZJets

keysStack             = [ "#gamma + jets","W/W* + jets" ,"single top"   ,"QCD multijets", "t#bar{t}"           ,  "Z/Z* + jets"  ]
stackColorIndexes     = [3               ,2             , 4             ,7              , 92                   ,  6              ]
stackFillStyleIds     = [3354            ,3395          , 3345          ,3345           , 3354                 ,  3345           ]

stackColorIndexes.reverse()
stackFillStyleIds.reverse()

samplesForHistos = ["LQ_M450"      , "LQ_M550"      , "LQ_M650"      ]
keys             = ["LQ, M=450 GeV", "LQ, M=550 GeV", "LQ, M=650 GeV"]

sampleForDataHisto = "DATA"

QCDScale   = 1.0
TTBarScale = 0.49

#--- nEle_PtCut_IDISO_noOvrlp ---

def makeDefaultPlot ( variableName, histoBaseName, 
                      samplesForHistos, keys,
                      samplesForStackHistos, keysStack,
                      sampleForDataHisto,
                      zUncBand, makeRatio ) :
    plot                   = Plot() 
    plot.histosStack       = generateHistoList( histoBaseName, samplesForStackHistos_other, variableName, File_preselection) + generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDScale) + generateHistoList( histoBaseName, samplesForStackHistosTTBar, variableName, File_TTBar, TTBarScale) + generateHistoList( histoBaseName, samplesForStackHistos_ZJets, variableName, File_preselection) 
    # plot.histosStack       = generateHistoList( histoBaseName, samplesForStackHistos_other, variableName, File_preselection) + generateHistoList( histoBaseName, samplesForStackHistosQCD, variableName, File_QCD, QCDScale) +  generateHistoList( histoBaseName, samplesForStackHistos_ZJets, variableName, File_preselection) 
    plot.keysStack         = keysStack
    plot.histos            = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
    plot.keys              = keys
    plot.addZUncBand       = zUncBand
    plot.makeRatio         = makeRatio
    plot.makeNSigma        = makeNSigma
    plot.histodata         = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)
    plot.ytit              = "Events"
    plot.ylog              = "no"
    plot.name              = variableName
    plot.stackColorIndexes = stackColorIndexes
    plot.stackFillStyleIds = stackFillStyleIds 
    plot.gif_folder        = "gif_eejj_scaled/"
    plot.eps_folder        = "eps_eejj_scaled/"
    plot.suffix            = "eejj"
    plot.lumi_pb           = "4623"
    
    return plot

plots = []

# plots.append ( makeDefaultPlot ( "nEle"      , histoBaseName, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
# plots[-1].ymax  = 10000000
# plots[-1].ymin  = 1e-1
# plots[-1].xmin  = -0.5
# plots[-1].xmax  = 6.5
# plots[-1].ylog  = "yes"
# plots[-1].xtit  = "Number of electrons (cut)"


plots.append ( makeDefaultPlot ( "GeneratorWeight"           ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio ))
plots[-1].xtit = "Generator level weight"
plots[-1].ymin  = 1e-1
plots[-1].makeNSigma = 0
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "PileupWeight"           ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio ))
plots[-1].xtit = "Pileup weight"
plots[-1].ymin  = 1e-1
plots[-1].makeNSigma = 0
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Ele1_Pt"  , histoBaseName, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = pt_rebin
plots[-1].ymax  = 20000000
plots[-1].ymin  = 1e-1
plots[-1].xmin  = 0
plots[-1].xmax  = 500
plots[-1].ylog  = "yes"
plots[-1].xtit  = "1st Electron p_{T} [GeV] (cut)"

plots.append ( makeDefaultPlot ( "Ele1_Eta"	 ,  histoBaseName, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "1st Electron #eta (cut)"   
plots[-1].ymax = 20000000
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Ele2_Pt"  , histoBaseName, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = pt_rebin
plots[-1].ymax  = 2000000
plots[-1].ymin  = 1e-1
plots[-1].xmin  = 0
plots[-1].xmax  = 500
plots[-1].ylog  = "yes"
plots[-1].xtit  = "2nd Electron p_{T} [GeV] (cut)"

plots.append ( makeDefaultPlot ( "Ele2_Eta"	 ,  histoBaseName, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "2nd Electron #eta (cut)"   
plots[-1].ymax = 2000000
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"


plots.append ( makeDefaultPlot ( "nJet"         , histoBaseName, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit  = "Number of jets"
plots[-1].ymax  = 200000000
plots[-1].ymin  = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xmin  = -0.5
plots[-1].xmax  = 10.5

plots.append ( makeDefaultPlot ( "Jet1_Pt"     , histoBaseName, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 5
plots[-1].ymax  = 2000000
plots[-1].ymin  = 1e-1
plots[-1].xmin  = 0
plots[-1].xmax  = 1000
plots[-1].ylog  = "yes"
plots[-1].xtit  = "1st Jet p_{T} [GeV] (cut) "

plots.append ( makeDefaultPlot ( "Jet2_Pt"     , histoBaseName, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = pt_rebin
plots[-1].ymax  = 100000
plots[-1].ymin  = 1e-1
plots[-1].xmin  = 0
plots[-1].xmax  = 1000
plots[-1].ylog  = "yes"
plots[-1].xtit  = "2nd Jet p_{T} [GeV] (cut)"

plots.append ( makeDefaultPlot ( "Jet1_Eta"         ,  histoBaseName, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].rebin = 5 
plots[-1].ylog  = "yes"
plots[-1].xtit = "1st Jet #eta (cut)"

plots.append ( makeDefaultPlot ( "Jet2_Eta"         ,  histoBaseName, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].rebin = 5 
plots[-1].ylog  = "yes"
plots[-1].xtit = "2nd Jet #eta (cut)"


# plots.append ( makeDefaultPlot ( "nMuon"      , histoBaseName, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
# plots[-1].ymax  = 10000000
# plots[-1].ymin  = 1e-1
# plots[-1].xmin  = -0.5
# plots[-1].xmax  = 6.5
# plots[-1].ylog  = "yes"
# plots[-1].xtit  = "Number of muons (cut)"

plots.append ( makeDefaultPlot ( "nElectron_PAS"      , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].ymax  = 10000000
plots[-1].ymin  = 1e-1
plots[-1].xmin  = -0.5
plots[-1].xmax  = 6.5
plots[-1].ylog  = "yes"
plots[-1].xtit  = "Number of electrons (preselection)"

plots.append ( makeDefaultPlot ( "nMuon_PAS"      , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].ymax  = 10000000
plots[-1].ymin  = 1e-1
plots[-1].xmin  = -0.5
plots[-1].xmax  = 6.5
plots[-1].ylog  = "yes"
plots[-1].xtit  = "Number of muons (preselection)"

plots.append ( makeDefaultPlot ( "nJet_PAS"         , histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit  = "Number of jets"
plots[-1].ymax  = 200000000
plots[-1].ymin  = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xmin  = -0.5
plots[-1].xmax  = 10.5

plots.append ( makeDefaultPlot ( "Pt1stEle_PAS"  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 2
plots[-1].ymax  = 2000000
plots[-1].ymin  = 1e-1
plots[-1].xmin  = 0
plots[-1].xmax  = 500
plots[-1].ylog  = "yes"
plots[-1].xtit  = "1st Electron p_{T} [GeV] (preselection)"

plots.append ( makeDefaultPlot ( "Eta1stEle_PAS"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "1st Electron #eta (preselection)"   
plots[-1].ymax = 2000000
plots[-1].rebin = 5
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Phi1stEle_PAS"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )   
plots[-1].xtit = "1st Electron #phi (preselection)"
plots[-1].rebin = 4
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Pt2ndEle_PAS"  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = pt_rebin
plots[-1].ymax  = 2000000
plots[-1].ymin  = 1e-1
plots[-1].xmin  = 0
plots[-1].xmax  = 500
plots[-1].ylog  = "yes"
plots[-1].xtit  = "2nd Electron p_{T} [GeV] (preselection)"

plots.append ( makeDefaultPlot ( "Eta2ndEle_PAS"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "2nd Electron #eta (preselection)"   
plots[-1].ymax = 2000000
plots[-1].rebin = 5
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Phi2ndEle_PAS"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )   
plots[-1].xtit = "2nd Electron #phi (preselection)"
plots[-1].rebin = 4
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Charge1stEle_PAS"      ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit = "1st Electron Charge (preselection)"
plots[-1].ymin = 0.0
plots[-1].ymax = 20000.

plots.append ( makeDefaultPlot ( "Charge2ndEle_PAS"      ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) ) 
plots[-1].xtit = "2nd Electron Charge (preselection)"
plots[-1].ymin = 0.0
plots[-1].ymax = 20000.



plots.append ( makeDefaultPlot ( "MTenu_PAS"               ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "M_{T}(e_{1}m PFMET (preselection)) [GeV]"
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmax = 500
plots[-1].xmin = 0
plots[-1].rebin = 2
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "MET_PAS"               ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "PFMET (preselection) [GeV]"
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmax = 500
plots[-1].xmin = 0
plots[-1].ylog  = "yes"


plots.append ( makeDefaultPlot ( "METSig_PAS"               ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "PFMET Significance (preselection) [GeV]"
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmax = 500
plots[-1].xmin = 0
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "METPhi_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "PFMET #phi (preselection)"
plots[-1].rebin = 1
plots[-1].rebin = 4
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"


plots.append ( makeDefaultPlot ( "METCharged_PAS"               ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "Charged PFMET (preselection) [GeV]"
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmax = 500
plots[-1].xmin = 0
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "METChargedPhi_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "Charged PFMET #phi (preselection)"
plots[-1].rebin = 4
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"


plots.append ( makeDefaultPlot ( "METType1_PAS"               ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "Type1 PFMET (preselection) [GeV]"
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmax = 500
plots[-1].xmin = 0
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "METType1Phi_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "Type1 PFMET #phi (preselection)"
plots[-1].rebin = 4
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Pt1stJet_PAS"          ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].xtit = "1st Jet p_{T} (preselection) [GeV]"
plots[-1].rebin = pt_rebin
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmax = 600
plots[-1].xmin = 0
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Pt2ndJet_PAS"          ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = pt_rebin
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmax = 600
plots[-1].xmin = 0
plots[-1].ylog  = "yes"
plots[-1].xtit = "2nd Jet p_{T} (preselection) [GeV]"

plots.append ( makeDefaultPlot ( "Eta1stJet_PAS"         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "1st Jet #eta (preselection)"

plots.append ( makeDefaultPlot ( "Eta2ndJet_PAS"         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "2nd Jet #eta (preselection)"

plots.append ( makeDefaultPlot ( "Phi1stJet_PAS"         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "1st Jet #phi (preselection)"

plots.append ( makeDefaultPlot ( "Phi2ndJet_PAS"	 ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 1
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "2nd Jet #phi (preselection)"

plots.append ( makeDefaultPlot ( "sTlep_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "S_{T} (1st Electron, 2nd Electron) (preselection) [GeV]"

plots.append ( makeDefaultPlot ( "sTjet_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "S_{T} (1st Jet, 2nd Jet) (preselection) [GeV]"

plots.append ( makeDefaultPlot ( "sT_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmax = 1000.
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "S_{T} (1st Electron, 2nd Electron, 1st Jet, 2nd Jet) (preselection) [GeV]"

plots.append ( makeDefaultPlot ( "Mjj_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].ylog  = "yes"
plots[-1].xtit = "Dijet Mass (preselection) [GeV]"

plots.append ( makeDefaultPlot ( "Meejj_PAS"	         ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )    
plots[-1].rebin = 5
plots[-1].ymax = 200000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "Mass_{eejj} (preselection) [GeV]"

plots.append ( makeDefaultPlot ( "Mee_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 20000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(ee) (preselection) [GeV]"

plots.append ( makeDefaultPlot ( "Mee_EBEB_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 3000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].xtit = "M(ee) (preselection, EB-EB) [GeV]"

plots.append ( makeDefaultPlot ( "Mee_EBEE_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 600
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].xtit = "M(ee) (preselection, EB-EE) [GeV]"


plots.append ( makeDefaultPlot ( "Mee_EEEE_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 300
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].xtit = "M(ee) (preselection, EE-EE) [GeV]"

plots.append ( makeDefaultPlot ( "Mee_EB_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 3000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].xtit = "M(ee) (preselection, EB-EE and EB-EB) [GeV]"


plots.append ( makeDefaultPlot ( "Mee_EBEB_80_100_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 3000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].xtit = "M(ee) [80, 100] (preselection, EB-EB) [GeV]"

plots.append ( makeDefaultPlot ( "Mee_EBEE_80_100_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 600
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].xtit = "M(ee) [80, 100] (preselection, EB-EE) [GeV]"


plots.append ( makeDefaultPlot ( "Mee_EEEE_80_100_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 300
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].xtit = "M(ee) [80, 100] (preselection, EE-EE) [GeV]"

plots.append ( makeDefaultPlot ( "Mee_EB_80_100_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 3000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].xtit = "M(ee) [80, 100] (preselection, EB-EE and EB-EB) [GeV]"

plots.append ( makeDefaultPlot ( "Mee_80_100_Preselection" ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].xmin = 70.
plots[-1].xmax = 110.
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(ee) [80, 100] (preselection) [GeV]"


plots.append ( makeDefaultPlot ( "Mee_EBEB_70_110_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 3000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].xtit = "M(ee) [70, 110] (preselection, EB-EB) [GeV]"

plots.append ( makeDefaultPlot ( "Mee_EBEE_70_110_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 600
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].xtit = "M(ee) [70, 110] (preselection, EB-EE) [GeV]"


plots.append ( makeDefaultPlot ( "Mee_EEEE_70_110_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 300
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].xtit = "M(ee) [70, 110] (preselection, EE-EE) [GeV]"

plots.append ( makeDefaultPlot ( "Mee_EB_70_110_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 3000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0.
plots[-1].xmax = 1000.
plots[-1].xtit = "M(ee) [70, 110] (preselection, EB-EE and EB-EB) [GeV]"

plots.append ( makeDefaultPlot ( "Mee_70_110_Preselection" ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].rebin = 4
plots[-1].xmin = 70.
plots[-1].xmax = 110.
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(ee) [70, 110] (preselection) [GeV]"

plots.append ( makeDefaultPlot ( "Ptee_PAS"  , histoBaseName_userDef , samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 2
plots[-1].ymax  = 2000000
plots[-1].ymin  = 1e-1
plots[-1].xmin  = 0
plots[-1].xmax  = 500
plots[-1].ylog  = "yes"
plots[-1].xtit  = "P_{T}(ee) (preselection) [GeV]}"

plots.append ( makeDefaultPlot ( "Me1j1_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin  = 0
plots[-1].xmax  = 1000
plots[-1].rebin = 2
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(e_{1}j_{1}) (preselection) [GeV]"

plots.append ( makeDefaultPlot ( "Me1j2_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin  = 0
plots[-1].xmax  = 1000
plots[-1].rebin = 2

plots[-1].ylog  = "yes"
plots[-1].xtit = "M(e_{1}j_{2}) (preselection) [GeV]"

plots.append ( makeDefaultPlot ( "Me2j1_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin  = 0
plots[-1].xmax  = 1000
plots[-1].rebin = 2

plots[-1].ylog  = "yes"
plots[-1].xtit = "M(e_{2}j_{1}) (preselection) [GeV]"

plots.append ( makeDefaultPlot ( "Me2j2_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "M(e_{2}j_{2}) (preselection) [GeV]"
plots[-1].xmin  = 0
plots[-1].xmax  = 1000
plots[-1].rebin = 2

plots.append ( makeDefaultPlot ( "Mej_selected_avg_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xtit = "M_{ej} (Preselection) [GeV]"
plots[-1].xmin  = 0
plots[-1].xmax  = 1000
plots[-1].rebin = 2

plots.append ( makeDefaultPlot ( "nVertex_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = 4000  
plots[-1].xmin = -0.5
plots[-1].xmax = 30.5
plots[-1].xtit = "n(vertexes) (preselection)"
# plots[-1].ylog  = "yes"
# plots[-1].ymax = 2000000
plots[-1].ymax = 3000

plots.append ( makeDefaultPlot ( "nVertex_good_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = 4000
plots[-1].xmin = -0.5
plots[-1].xmax = 30.5
plots[-1].xtit = "n(good vertexes) (preselection)"
# plots[-1].ylog  = "yes"
plots[-1].ymax = 3000


plots.append ( makeDefaultPlot ( "DR_Ele1Jet1_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].xtit = "#DeltaR(e_{1},j_{1})"
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0
plots[-1].xmax = 6
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "DR_Ele1Jet2_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].xtit = "#DeltaR(e_{1},j_{2})"
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0
plots[-1].xmax = 6
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "DR_Ele2Jet1_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].xtit = "#DeltaR(e_{2},j_{1})"
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0
plots[-1].xmax = 6
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "DR_Ele2Jet2_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].xtit = "#DeltaR(e_{2},j_{2})"
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0
plots[-1].xmax = 6
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "DR_Ele1Ele2_PAS"                ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].xtit = "#DeltaR(e_{1},e_{2})"
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].xmin = 0
plots[-1].xmax = 6
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "minDR_EleJet_PAS"           ,  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].xtit = "Minimum #DeltaR(e_{1},(j_{1}, j_{2}, j_{3})) (cut)"
plots[-1].ymax = 2000000
plots[-1].ymin = 1e-1
plots[-1].ylog  = "yes"
plots[-1].xmin = 0
plots[-1].xmax = 6

Mej_selected_avg_max = 20000
Mej_selected_min_max = 20000
sT_eejj_max          = 20000
Mee_max              = 20000

plots.append ( makeDefaultPlot ( "Mej_selected_avg_LQ250", histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = Mej_selected_avg_max
plots[-1].xtit = "M_{avg}(ej) [GeV], LQ M=250 optimization"
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Mej_selected_avg_LQ350",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = Mej_selected_avg_max
plots[-1].xtit = "M_{avg}(ej) [GeV], LQ M=350 optimization"
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Mej_selected_avg_LQ400",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = Mej_selected_avg_max
plots[-1].xtit = "M_{avg}(ej) [GeV], LQ M=400 optimization"
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Mej_selected_avg_LQ450",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = Mej_selected_avg_max
plots[-1].xtit = "M_{avg}(ej) [GeV], LQ M=450 optimization"
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Mej_selected_avg_LQ500",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = Mej_selected_avg_max
plots[-1].xtit = "M_{avg}(ej) [GeV], LQ M=500 optimization"
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Mej_selected_avg_LQ550",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = Mej_selected_avg_max
plots[-1].xtit = "M_{avg}(ej) [GeV], LQ M=550 optimization"
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Mej_selected_avg_LQ600",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = Mej_selected_avg_max
plots[-1].xtit = "M_{avg}(ej) [GeV], LQ M=600 optimization"
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Mej_selected_avg_LQ650",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = Mej_selected_avg_max
plots[-1].xtit = "M_{avg}(ej) [GeV], LQ M=650 optimization"
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Mej_selected_avg_LQ750",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = Mej_selected_avg_max
plots[-1].xtit = "M_{avg}(ej) [GeV], LQ M=750 optimization"
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Mej_selected_avg_LQ850",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = Mej_selected_avg_max
plots[-1].xtit = "M_{avg}(ej) [GeV], LQ M=850 optimization"
plots[-1].ylog  = "yes"


plots.append ( makeDefaultPlot ( "Mej_selected_min_LQ250", histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = Mej_selected_min_max
plots[-1].xtit = "M_{min}(ej) [GeV], LQ M=250 optimization"
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Mej_selected_min_LQ350",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = Mej_selected_min_max
plots[-1].xtit = "M_{min}(ej) [GeV], LQ M=350 optimization"
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Mej_selected_min_LQ400",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = Mej_selected_min_max
plots[-1].xtit = "M_{min}(ej) [GeV], LQ M=400 optimization"
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Mej_selected_min_LQ450",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = Mej_selected_min_max
plots[-1].xtit = "M_{min}(ej) [GeV], LQ M=450 optimization"
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Mej_selected_min_LQ500",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = Mej_selected_min_max
plots[-1].xtit = "M_{min}(ej) [GeV], LQ M=500 optimization"
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Mej_selected_min_LQ550",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = Mej_selected_min_max
plots[-1].xtit = "M_{min}(ej) [GeV], LQ M=550 optimization"
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Mej_selected_min_LQ600",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = Mej_selected_min_max
plots[-1].xtit = "M_{min}(ej) [GeV], LQ M=600 optimization"
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Mej_selected_min_LQ650",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = Mej_selected_min_max
plots[-1].xtit = "M_{min}(ej) [GeV], LQ M=650 optimization"
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Mej_selected_min_LQ750",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = Mej_selected_min_max
plots[-1].xtit = "M_{min}(ej) [GeV], LQ M=750 optimization"
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Mej_selected_min_LQ850",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = Mej_selected_min_max
plots[-1].xtit = "M_{min}(ej) [GeV], LQ M=850 optimization"
plots[-1].ylog  = "yes"



plots.append ( makeDefaultPlot ( "sT_eejj_LQ250", histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = sT_eejj_max
plots[-1].xtit = "s_{T} [GeV], LQ M=250 optimization"
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "sT_eejj_LQ350",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = sT_eejj_max
plots[-1].xtit = "s_{T} [GeV], LQ M=350 optimization"
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "sT_eejj_LQ400",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = sT_eejj_max
plots[-1].xtit = "s_{T} [GeV], LQ M=400 optimization"
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "sT_eejj_LQ450",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = sT_eejj_max
plots[-1].xtit = "s_{T} [GeV], LQ M=450 optimization"
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "sT_eejj_LQ500",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = sT_eejj_max
plots[-1].xtit = "s_{T} [GeV], LQ M=500 optimization"
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "sT_eejj_LQ550",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = sT_eejj_max
plots[-1].xtit = "s_{T} [GeV], LQ M=550 optimization"
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "sT_eejj_LQ600",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = sT_eejj_max
plots[-1].xtit = "s_{T} [GeV], LQ M=600 optimization"
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "sT_eejj_LQ650",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = sT_eejj_max
plots[-1].xtit = "s_{T} [GeV], LQ M=650 optimization"
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "sT_eejj_LQ750",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = sT_eejj_max
plots[-1].xtit = "s_{T} [GeV], LQ M=750 optimization"
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "sT_eejj_LQ850",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = sT_eejj_max
plots[-1].xtit = "s_{T} [GeV], LQ M=850 optimization"
plots[-1].ylog  = "yes"



plots.append ( makeDefaultPlot ( "Mee_LQ250", histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = Mee_max
plots[-1].xtit = "M(ee) [GeV], LQ M=250 optimization"
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Mee_LQ350",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = Mee_max
plots[-1].xtit = "M(ee) [GeV], LQ M=350 optimization"
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Mee_LQ400",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = Mee_max
plots[-1].xtit = "M(ee) [GeV], LQ M=400 optimization"
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Mee_LQ450",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = Mee_max
plots[-1].xtit = "M(ee) [GeV], LQ M=450 optimization"
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Mee_LQ500",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = Mee_max
plots[-1].xtit = "M(ee) [GeV], LQ M=500 optimization"
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Mee_LQ550",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = Mee_max
plots[-1].xtit = "M(ee) [GeV], LQ M=550 optimization"
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Mee_LQ600",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = Mee_max
plots[-1].xtit = "M(ee) [GeV], LQ M=600 optimization"
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Mee_LQ650",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = Mee_max
plots[-1].xtit = "M(ee) [GeV], LQ M=650 optimization"
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Mee_LQ750",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = Mee_max
plots[-1].xtit = "M(ee) [GeV], LQ M=750 optimization"
plots[-1].ylog  = "yes"

plots.append ( makeDefaultPlot ( "Mee_LQ850",  histoBaseName_userDef, samplesForHistos, keys, samplesForStackHistos, keysStack, sampleForDataHisto, zUncBand, makeRatio) )
plots[-1].rebin = 1
plots[-1].ymin = 1e-1
plots[-1].ymax = Mee_max
plots[-1].xtit = "M(ee) [GeV], LQ M=850 optimization"
plots[-1].ylog  = "yes"







#--- Generate and print the plots from the list 'plots' define above
c = TCanvas()

fileps = "allPlots_eejj_scaled_analysis.ps"

c.Print(fileps + "[")
for plot in plots:
    plot.Draw(fileps)
c.Print(fileps+"]")
os.system('ps2pdf '+fileps)
os.system('rm '+fileps)