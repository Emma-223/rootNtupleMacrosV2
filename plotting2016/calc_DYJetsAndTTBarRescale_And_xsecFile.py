#!/usr/bin/env python3

##############################################################################
# # USER CODE IS TOWARD THE END OF THE FILE
##############################################################################

##############################################################################
# ############ DON'T NEED TO MODIFY ANYTHING HERE - BEGIN #####################

import sys
import string
import os.path
from ROOT import kTRUE, gROOT, gStyle, TFile, TCanvas, TRandom3, kWhite, TGraphErrors, kBlue, TLegend
import ROOT
import re
import copy
import math
import numpy as np
from tabulate import tabulate

# --- ROOT general options
gROOT.SetBatch(kTRUE)
gStyle.SetOptStat(0)
gStyle.SetPalette(1)
gStyle.SetCanvasBorderMode(0)
gStyle.SetFrameBorderMode(0)
gStyle.SetCanvasColor(kWhite)
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
# --- TODO: WHY IT DOES NOT LOAD THE DEFAULT ROOTLOGON.C ? ---#


def GetFile(filename):
    filename = filename.replace("/eos/cms", "root://eoscms/").replace("/eos/user", "root://eosuser/")
    tfile = TFile.Open(filename)
    if not tfile or tfile.IsZombie():
        print("ERROR: file " + filename + " not found")
        print("exiting...")
        sys.exit(-1)
    return tfile


def GetHisto(histoName, tfile, scale=1):
    histo = tfile.Get(histoName)
    if not histo:
        # in this case, try for one with systematics
        histoWithSystsName = histoName.replace("histo1D", "histo2D")+"WithSystematics"
        histo2D = tfile.Get(histoWithSystsName)
        if not histo2D:
            raise RuntimeError("ERROR: neither histo {} nor histo {} were found in file {}.".format(histoName, histoWithSystsName, tfile.GetName()))
        print("INFO: Using histo named {} and projecting to 1-D".format(histoWithSystsName))
        histo = histo2D.ProjectionX(histoName, 1, 1)  # convert to 1-D nominal hist
    new = copy.deepcopy(histo)
    if scale != 1:
        new.Scale(scale)
    return new


def Create1DSlice(hist, xmin, xmax):
    histName = hist.GetName() + "_" + str(xmin) + "to" + str(xmax)
    minXBin = hist.GetXaxis().FindBin(xmin)
    if xmax == "Inf":
        maxXBin = hist.GetNbinsX()+1  # include overflow
    else:
        maxXBin = hist.GetXaxis().FindBin(xmax-1e-10)  # in case xmax is exactly a bin boundary, go a tiny bit inside the bin to get the right bin number
    # print("INFO: Create1DSlice() - for hist named {}, got xmin={} to bin={} [center={}] and xmax={} to bin={} [center={}]".format(hist.GetName(), xmin, minXBin, hist.GetXaxis().GetBinCenter(minXBin), xmax, maxXBin, hist.GetXaxis().GetBinCenter(maxXBin)), flush=True)
    return hist.ProjectionY(histName, minXBin, maxXBin, "e")


def FillPlotObjects(plotBaseName, suffix, meeMin, meeMax, xmax, datasetName, h_DATA, h_ALLBKG, h_TTbar, h_ZJets, h_SingleTop, h_Diboson, h_QCD, h_WJets):
    plot = Plot()
    plot.histoDATA = copy.deepcopy(h_DATA)
    plot.histoMCall = copy.deepcopy(h_ALLBKG)
    plot.histoTTbar = copy.deepcopy(h_TTbar)
    if doQCD:
        plot.histoQCD = copy.deepcopy(h_QCD)
    plot.histoZJet = copy.deepcopy(h_ZJets)
    # plot.histoWJet = h_WJets # FIXME
    plot.histoSingleTop = copy.deepcopy(h_SingleTop)
    plot.histoDiboson = copy.deepcopy(h_Diboson)
    plot.xmin = meeMin
    plot.xmax = meeMax
    plot.name = plotBaseName + suffix
    plot.fileXsectionNoRescale = xsectionFile
    plot.xminplot = 0
    plot.xmaxplot = xmax
    plot.yminplot = 0
    plot.ymaxplot = 2000
    plot.datasetName = datasetName
    return plot


def SetupPlots(thisHistName, histBaseName, bins=[]):
    ttbarPlotBaseName = histBaseName.replace("BJETBIN1", "gteOneBtaggedJet")
    ttbarPlotBaseName = ttbarPlotBaseName.replace("BJETBIN2", "gteTwoBtaggedJets")
    plotsTTbar = []
    # plotBaseName = histBaseName.replace("BJETBIN1", "noBtaggedJets")
    # plotBaseName = plotBaseName.replace("BJETBIN2", "noBtaggedJets")
    dyjPlotBaseName = histBaseName.replace("_BJETBIN1", "")
    dyjPlotBaseName = dyjPlotBaseName.replace("_BJETBIN2", "")
    plotsDYJets = []
    backgroundNames = ["DATA", allBkg, ttbar, zjet, singletop, diboson, qcd, wjet]
    for idx, plotBaseName in enumerate([ttbarPlotBaseName, dyjPlotBaseName]):
        isDYJ = True
        sampleName = "DYJets"
        if idx == 0:
            sampleName = "TTBar"
            isDYJ = False
        print("INFO: For {}, using plotBaseName: {}".format(sampleName, plotBaseName))
        histDict = {}
        histDict[ttbar] = GetHisto(
            thisHistName.replace("SAMPLE", ttbar) + plotBaseName, File_preselection
        )
        histDict[zjet] = GetHisto(
            thisHistName.replace("SAMPLE", zjet) + plotBaseName,
            File_preselection,
        )
        histDict[singletop] = GetHisto(
            thisHistName.replace("SAMPLE", singletop) + plotBaseName, File_preselection
        )
        histDict[diboson] = GetHisto(
            thisHistName.replace("SAMPLE", diboson) + plotBaseName,
            File_preselection,
        )
        histDict[allBkg] = copy.deepcopy(histDict[diboson].Clone(allBkg))
        histDict[allBkg].Add(histDict[ttbar])
        histDict[allBkg].Add(histDict[zjet])
        histDict[allBkg].Add(histDict[singletop])
        histDict["DATA"] = GetHisto(
            thisHistName.replace("SAMPLE", "DATA") + plotBaseName, File_preselection
            )
        histDict[qcd] = None
        histDict[wjet] = None
        if doQCD:
            histDict[qcd] = GetHisto(
                thisHistName.replace("SAMPLE", qcd)
                #+ plotBaseName.replace("_btagSFDownShift", "").replace("_btagSFUpShift", ""),
                + plotBaseName,
                File_QCD_preselection,
            )
            # h_QCD = GetHisto("histo1D__QCD_EMEnriched__"+plotBaseName,File_QCD_preselection)
        elif useSingleFakeMC:
            histDict[wjet] = GetHisto(
                thisHistName.replace("SAMPLE", wjet) + plotBaseName,
                File_preselection,
            )
        if len(bins):
            for idx, binCoord in enumerate(bins):
                if idx == len(bins)-1:
                    break
                binHigh = bins[idx+1]
                # print("INFO: make plot from {} from {} to {}; histDict[zjet] has {} entries".format(plotBaseName, binCoord, binHigh, histDict[zjet].GetEntries()), flush=True)
                thisPlotBaseName = plotBaseName + "_" + str(binCoord) + "to" + str(binHigh)
                histsForBinRange = []
                for bkg in backgroundNames:
                    if histDict[bkg] is None:
                        histsForBinRange.append(histDict[bkg])
                        continue
                    sliceHist = Create1DSlice(histDict[bkg], binCoord, binHigh)
                    xmax = histDict[bkg].GetXaxis().GetXmax()
                    histsForBinRange.append(sliceHist)
                if isDYJ:
                    plotsDYJets.append(FillPlotObjects(thisPlotBaseName, "_DYJets", meeMinDYJets, meeMaxDYJets, xmax, zjetDatasetName, *histsForBinRange))
                else:
                    plotsTTbar.append(FillPlotObjects(thisPlotBaseName, "_TTbar", meeMinTTBar, meeMaxTTBar, xmax, ttbarDatasetName, *histsForBinRange))
        else:
            if isDYJ:
                plotsDYJets.append(FillPlotObjects(plotBaseName, "_DYJets", meeMinDYJets, meeMaxDYJets, meePlotMaxDYJets, zjetDatasetName, *[histDict[bkg] for bkg in backgroundNames]))
            else:
                plotsTTbar.append(FillPlotObjects(plotBaseName, "_TTbar", meeMinTTBar, meeMaxTTBar, meePlotMaxTTBar, ttbarDatasetName, *[histDict[bkg] for bkg in backgroundNames]))
    return plotsDYJets, plotsTTbar


def GetIntegralTH1(histo, xmin, xmax, verbose=False):
    # get integral
    axis = histo.GetXaxis()
    bmin = axis.FindBin(xmin)
    bmax = axis.FindBin(xmax)
    bminResidual = (
        histo.GetBinContent(bmin)
        * (xmin - axis.GetBinLowEdge(bmin))
        / axis.GetBinWidth(bmin)
    )
    bmaxResidual = (
        histo.GetBinContent(bmax)
        * (axis.GetBinUpEdge(bmax) - xmax)
        / axis.GetBinWidth(bmax)
    )
    integral = histo.Integral(bmin, bmax) - bminResidual - bmaxResidual
    if verbose:
        print("GetIntegralTH1(" + histo.GetName(), xmin, xmax, ")=", integral)
    return integral


def GetErrorIntegralTH1(histo, xmin, xmax, verbose=False):
    if verbose:
        print("## calculating error for integral of histo " + str(histo))
        print("## in the x range [" + str(xmin) + "," + str(xmax) + "]")
    # get integral
    axis = histo.GetXaxis()
    bmin = axis.FindBin(xmin)
    bmax = axis.FindBin(xmax)
    bminResidual = (
        histo.GetBinContent(bmin)
        * (xmin - axis.GetBinLowEdge(bmin))
        / axis.GetBinWidth(bmin)
    )
    bmaxResidual = (
        histo.GetBinContent(bmax)
        * (axis.GetBinUpEdge(bmax) - xmax)
        / axis.GetBinWidth(bmax)
    )
    integral = histo.Integral(bmin, bmax) - bminResidual - bmaxResidual
    error = 0
    for bin in range(bmin, bmax + 1):
        # print "bin: " +str(bin)
        if bin == bmax and bmaxResidual == histo.GetBinContent(
            bmax
        ):  # skip last bin if out of range
            if verbose:
                print("     --> skip bin: " + str(bin))
        else:
            error = error + histo.GetBinError(bin) ** 2
            # print "error**2 : " + str(error)

    error = math.sqrt(error)
    if verbose:
        print(" ")
    return error


# The Plot class: add members if needed
class Plot:
    def __init__(self):
        histoDATA = ""  # DATA
        histoTTbar = ""  # MCTTbar
        histoMCall = ""  # MCall
        histoQCD = ""  # QCD
        histoZJet = ""
        histoWJet = ""
        histoSingleTop = ""
        histoDiboson = ""
        xtit = ""  # xtitle
        ytit = ""  # ytitle
        xmin = ""  # set xmin to calculate rescaling factor (-- please take into account the bin size of histograms --)
        xmax = ""  # # set xmax to calculate rescaling factor (-- please take into account the bin size of histograms --)
        xminplot = ""  # min x axis range (need to set both min and max. Leave it as is for full range)
        xmaxplot = ""  # max x axis range (need to set both min and max. Leave it as is for full range)
        yminplot = ""  # min y axis range (need to set both min and max. Leave it as is for full range)
        ymaxplot = ""  # max y axis range (need to set both min and max. Leave it as is for full range)
        lpos = (
            ""  # legend position (default = top-right, option="bottom-center", "top-left")
        )
        #    xlog         = "" # log scale of X axis (default = no, option="yes") ### IT SEEMS IT DOES NOT WORK
        ylog = ""  # log scale of Y axis (default = no, option="yes")
        # rebin       = "" # rebin x axis (default = 1, option = set it to whatever you want )
        name = ""  # name of the final plots
        lint = "2.6 fb^{-1}"  # integrated luminosity of the sample ( example "10 pb^{-1}" )
        fileXsectionNoRescale = ""  # cross section file (with no rescale
        writeXSecFile = False
        datasetName = ""  # string for pattern recognition of dataset name (rescaling will be done only on matched datasets)

    def CheckMCDataConsistency(self):
        # checks
        if self.histoMCall.GetNbinsX() != self.histoDATA.GetNbinsX():
            print("WARNING! number of bins is different between DATA and MC")
            print("exiting...")
            sys.exit()
        if self.histoMCall.GetBinWidth(1) != self.histoDATA.GetBinWidth(1):
            print("WARNING! bin width is different between DATA and MC")
            print("exiting...")
            sys.exit()


def GetRTTBarDYJets(plotObjTTBar, plotObjDYJets, randomize=False, verbose=False):

    # integrals: ttbar
    integralDATA_ttbar = GetIntegralTH1(
        plotObjTTBar.histoDATA, plotObjTTBar.xmin, plotObjTTBar.xmax
    )
    ERRintegralDATA_ttbar = GetErrorIntegralTH1(
        plotObjTTBar.histoDATA, plotObjTTBar.xmin, plotObjTTBar.xmax
    )
    integralMCall_ttbar = GetIntegralTH1(
        plotObjTTBar.histoMCall, plotObjTTBar.xmin, plotObjTTBar.xmax
    )
    ERRintegralMCall_ttbar = GetErrorIntegralTH1(
        plotObjTTBar.histoMCall, plotObjTTBar.xmin, plotObjTTBar.xmax
    )
    integralTTbar_ttbar = GetIntegralTH1(
        plotObjTTBar.histoTTbar, plotObjTTBar.xmin, plotObjTTBar.xmax
    )
    ERRintegralTTbar_ttbar = GetErrorIntegralTH1(
        plotObjTTBar.histoTTbar, plotObjTTBar.xmin, plotObjTTBar.xmax
    )
    integralDYJets_ttbar = GetIntegralTH1(
        plotObjTTBar.histoZJet, plotObjTTBar.xmin, plotObjTTBar.xmax
    )
    ERRintegralDYJets_ttbar = GetErrorIntegralTH1(
        plotObjTTBar.histoZJet, plotObjTTBar.xmin, plotObjTTBar.xmax
    )
    if doQCD:
        integralQCD_ttbar = GetIntegralTH1(
            plotObjTTBar.histoQCD, plotObjTTBar.xmin, plotObjTTBar.xmax
        )
        ERRintegralQCD_ttbar = GetErrorIntegralTH1(
            plotObjTTBar.histoQCD, plotObjTTBar.xmin, plotObjTTBar.xmax
        )
    else:
        integralQCD_ttbar = 0
        ERRintegralQCD_ttbar = 0
    # contamination from other backgrounds (except TTbar and DYJets) in the integral range [QCD is not in MCall]
    integralMCothers_ttbar = (
        integralMCall_ttbar - integralTTbar_ttbar - integralDYJets_ttbar
    )
    ERRintegralMCothers_ttbar = math.sqrt(
        ERRintegralMCall_ttbar ** 2 + ERRintegralTTbar_ttbar ** 2
    )

    # integrals: wjets
    integralDATA_dyjets = GetIntegralTH1(
        plotObjDYJets.histoDATA, plotObjDYJets.xmin, plotObjDYJets.xmax
    )
    ERRintegralDATA_dyjets = GetErrorIntegralTH1(
        plotObjDYJets.histoDATA, plotObjDYJets.xmin, plotObjDYJets.xmax
    )
    integralMCall_dyjets = GetIntegralTH1(
        plotObjDYJets.histoMCall, plotObjDYJets.xmin, plotObjDYJets.xmax
    )
    ERRintegralMCall_dyjets = GetErrorIntegralTH1(
        plotObjDYJets.histoMCall, plotObjDYJets.xmin, plotObjDYJets.xmax
    )
    integralTTbar_dyjets = GetIntegralTH1(
        plotObjDYJets.histoTTbar, plotObjDYJets.xmin, plotObjDYJets.xmax
    )
    ERRintegralTTbar_dyjets = GetErrorIntegralTH1(
        plotObjDYJets.histoTTbar, plotObjDYJets.xmin, plotObjDYJets.xmax
    )
    integralDYJets_dyjets = GetIntegralTH1(
        plotObjDYJets.histoZJet, plotObjDYJets.xmin, plotObjDYJets.xmax
    )
    ERRintegralDYJets_dyjets = GetErrorIntegralTH1(
        plotObjDYJets.histoZJet, plotObjDYJets.xmin, plotObjDYJets.xmax
    )
    if doQCD:
        integralQCD_dyjets = GetIntegralTH1(
            plotObjDYJets.histoQCD, plotObjDYJets.xmin, plotObjDYJets.xmax
        )
        ERRintegralQCD_dyjets = GetErrorIntegralTH1(
            plotObjDYJets.histoQCD, plotObjDYJets.xmin, plotObjDYJets.xmax
        )
    else:
        integralQCD_dyjets = 0
        ERRintegralQCD_dyjets = 0
    # contamination from other backgrounds (except dyjets and TTBar) in the integral range [QCD is not in MCall]
    integralMCothers_dyjets = (
        integralMCall_dyjets - integralDYJets_dyjets - integralTTbar_dyjets
    )
    ERRintegralMCothers_dyjets = math.sqrt(
        ERRintegralMCall_dyjets ** 2 + ERRintegralDYJets_dyjets ** 2
    )
    contamination_dyjets = (
        integralMCothers_dyjets + integralTTbar_dyjets + integralQCD_dyjets
    ) / (integralMCall_dyjets + integralQCD_dyjets)

    # randomize
    trand = TRandom3(0)
    integralDATA_ttbar = trand.Gaus(integralDATA_ttbar, ERRintegralDATA_ttbar)
    integralMCall_ttbar = trand.Gaus(integralMCall_ttbar, ERRintegralMCall_ttbar)
    integralTTbar_ttbar = trand.Gaus(integralTTbar_ttbar, ERRintegralTTbar_ttbar)
    integralDYJets_ttbar = trand.Gaus(integralDYJets_ttbar, ERRintegralDYJets_ttbar)
    integralQCD_ttbar = trand.Gaus(integralQCD_ttbar, ERRintegralQCD_ttbar)
    integralMCothers_ttbar = (
        integralMCall_ttbar - integralTTbar_ttbar - integralDYJets_ttbar
    )
    #
    integralDATA_dyjets = trand.Gaus(integralDATA_dyjets, ERRintegralDATA_dyjets)
    integralMCall_dyjets = trand.Gaus(integralMCall_dyjets, ERRintegralMCall_dyjets)
    integralTTbar_dyjets = trand.Gaus(integralTTbar_dyjets, ERRintegralTTbar_dyjets)
    integraldyjets_dyjets = trand.Gaus(integralDYJets_dyjets, ERRintegralDYJets_dyjets)
    integralQCD_dyjets = trand.Gaus(integralQCD_dyjets, ERRintegralQCD_dyjets)

    # solve the system of equations
    # (1) --> dyjets
    # (2) --> ttbar
    rTTBar = (
        integralDATA_dyjets * integralDYJets_ttbar
        - (integralMCothers_dyjets + integralQCD_dyjets) * integralDYJets_ttbar
    )
    rTTBar += (
        integralMCothers_ttbar + integralQCD_ttbar - integralDATA_ttbar
    ) * integraldyjets_dyjets
    if verbose:
        print(integralTTbar_dyjets, integralDYJets_ttbar, integralTTbar_ttbar, integralDYJets_dyjets)
    try:
        rTTBar /= (
            integralTTbar_dyjets * integralDYJets_ttbar
            - integralTTbar_ttbar * integralDYJets_dyjets
        )
    except ZeroDivisionError:
        print("ERROR: ZeroDivisionError: one or more of the integrals above are zero")
        return -1, -1
    rDYJets = (
        integralDATA_dyjets * integralTTbar_ttbar
        - (integralMCothers_dyjets + integralQCD_dyjets) * integralTTbar_ttbar
    )
    rDYJets += (
        integralMCothers_ttbar + integralQCD_ttbar - integralDATA_ttbar
    ) * integralTTbar_dyjets
    rDYJets /= (
        integralTTbar_ttbar * integralDYJets_dyjets
        - integralTTbar_dyjets * integralDYJets_ttbar
    )

    return rTTBar, rDYJets


def CalculateRescaleFactor(plotObjTTBar, plotObjDYJets, fileps):
    # calculate rescaling factor for Z/gamma+jet background and create new cross section file
    canvas = TCanvas()

    plotObjTTBar.CheckMCDataConsistency()
    plotObjDYJets.CheckMCDataConsistency()

    # do the following N times, so we can calculate the stat. uncertainty
    N = 100
    # NB: the commented code below makes a nice progress bar but causes the dict to be undefined...
    # steps = N
    print("Randomizing histos and calculating scale factors:", end=' ')
    # progressString = "0% [" + " " * steps + "] 100%"
    # print progressString,
    # print "\b" * (len(progressString) - 3),
    sys.stdout.flush()
    rTTBarList = []
    rDYJetsList = []
    for i in range(0, N):
        # print ""
        rTTBar, rDYJets = GetRTTBarDYJets(plotObjTTBar, plotObjDYJets, True)
        rTTBarList.append(rTTBar)
        rDYJetsList.append(rDYJets)
        # print "\b.",
        # sys.stdout.flush()
    # print "\b] 100%"
    print("Done.")

    ttMean = 0
    for rtt in rTTBarList:
        ttMean += rtt
    ttMean /= N
    #
    wMean = 0
    for rw in rDYJetsList:
        wMean += rw
    wMean /= N

    rTTBarSigma = 0
    for rtt in rTTBarList:
        rTTBarSigma += pow(rtt - ttMean, 2)
    rTTBarSigma /= N
    rTTBarSigma = math.sqrt(rTTBarSigma)
    #
    rDYJetsSigma = 0
    for rdy in rDYJetsList:
        rDYJetsSigma += pow(rdy - wMean, 2)
    rDYJetsSigma /= N
    rDYJetsSigma = math.sqrt(rDYJetsSigma)

    # integrals: ttbar
    integralDATA_ttbar = GetIntegralTH1(
        plotObjTTBar.histoDATA, plotObjTTBar.xmin, plotObjTTBar.xmax
    )
    ERRintegralDATA_ttbar = GetErrorIntegralTH1(
        plotObjTTBar.histoDATA, plotObjTTBar.xmin, plotObjTTBar.xmax
    )
    integralMCall_ttbar = GetIntegralTH1(
        plotObjTTBar.histoMCall, plotObjTTBar.xmin, plotObjTTBar.xmax
    )
    ERRintegralMCall_ttbar = GetErrorIntegralTH1(
        plotObjTTBar.histoMCall, plotObjTTBar.xmin, plotObjTTBar.xmax
    )
    integralTTbar_ttbar = GetIntegralTH1(
        plotObjTTBar.histoTTbar, plotObjTTBar.xmin, plotObjTTBar.xmax
    )
    # print "plotObjTTBar.histoTTbar=", plotObjTTBar.histoTTbar
    ERRintegralTTbar_ttbar = GetErrorIntegralTH1(
        plotObjTTBar.histoTTbar, plotObjTTBar.xmin, plotObjTTBar.xmax
    )
    integralDYJets_ttbar = GetIntegralTH1(
        plotObjTTBar.histoZJet, plotObjTTBar.xmin, plotObjTTBar.xmax
    )
    ERRintegralDYJets_ttbar = GetErrorIntegralTH1(
        plotObjTTBar.histoZJet, plotObjTTBar.xmin, plotObjTTBar.xmax
    )
    if doQCD:
        integralQCD_ttbar = GetIntegralTH1(
            plotObjTTBar.histoQCD, plotObjTTBar.xmin, plotObjTTBar.xmax
        )
        ERRintegralQCD_ttbar = GetErrorIntegralTH1(
            plotObjTTBar.histoQCD, plotObjTTBar.xmin, plotObjTTBar.xmax
        )
    else:
        integralQCD_ttbar = 0
        ERRintegralQCD_ttbar = 0
    # contamination from other backgrounds (except TTbar and DYJets) in the integral range [QCD is not in MCall]
    integralMCothers_ttbar = (
        integralMCall_ttbar - integralTTbar_ttbar - integralDYJets_ttbar
    )
    ERRintegralMCothers_ttbar = math.sqrt(
        ERRintegralMCall_ttbar ** 2 + ERRintegralTTbar_ttbar ** 2
    )
    try:
        contamination_ttbar = (
            integralMCothers_ttbar + integralDYJets_ttbar + integralQCD_ttbar
        ) / (integralMCall_ttbar + integralQCD_ttbar)
    except ZeroDivisionError:
        print("ERROR: ZeroDivisionError: integralMCall_ttbar+integralQCD_ttbar is zero; integralMCall_ttbar=", integralMCall_ttbar, "integralQCD_ttbar=", integralQCD_ttbar)
        contamination_ttbar = -1

    # integrals: dyjets
    integralDATA_dyjets = GetIntegralTH1(
        plotObjDYJets.histoDATA, plotObjDYJets.xmin, plotObjDYJets.xmax
    )
    ERRintegralDATA_dyjets = GetErrorIntegralTH1(
        plotObjDYJets.histoDATA, plotObjDYJets.xmin, plotObjDYJets.xmax
    )
    integralMCall_dyjets = GetIntegralTH1(
        plotObjDYJets.histoMCall, plotObjDYJets.xmin, plotObjDYJets.xmax
    )
    ERRintegralMCall_dyjets = GetErrorIntegralTH1(
        plotObjDYJets.histoMCall, plotObjDYJets.xmin, plotObjDYJets.xmax
    )
    integralTTbar_dyjets = GetIntegralTH1(
        plotObjDYJets.histoTTbar, plotObjDYJets.xmin, plotObjDYJets.xmax
    )
    # print "plotObjDYJets.histoTTbar=", plotObjDYJets.histoTTbar
    ERRintegralTTbar_dyjets = GetErrorIntegralTH1(
        plotObjDYJets.histoTTbar, plotObjDYJets.xmin, plotObjDYJets.xmax
    )
    integralDYJets_dyjets = GetIntegralTH1(
        plotObjDYJets.histoZJet, plotObjDYJets.xmin, plotObjDYJets.xmax
    )
    ERRintegralDYJets_dyjets = GetErrorIntegralTH1(
        plotObjDYJets.histoZJet, plotObjDYJets.xmin, plotObjDYJets.xmax
    )
    if doQCD:
        integralQCD_dyjets = GetIntegralTH1(
            plotObjDYJets.histoQCD, plotObjDYJets.xmin, plotObjDYJets.xmax
        )
        ERRintegralQCD_dyjets = GetErrorIntegralTH1(
            plotObjDYJets.histoQCD, plotObjDYJets.xmin, plotObjDYJets.xmax
        )
    else:
        integralQCD_dyjets = 0
        ERRintegralQCD_dyjets = 0
    # contamination from other backgrounds (except WJets and TTBar) in the integral range [QCD is not in MCall]
    integralMCothers_dyjets = (
        integralMCall_dyjets - integralDYJets_dyjets - integralTTbar_dyjets
    )
    ERRintegralMCothers_dyjets = math.sqrt(
        ERRintegralMCall_dyjets ** 2 + ERRintegralDYJets_dyjets ** 2
    )
    contamination_dyjets = (
        integralMCothers_dyjets + integralTTbar_dyjets + integralQCD_dyjets
    ) / (integralMCall_dyjets + integralQCD_dyjets)

    # solve the system of equations
    # (1) --> dyjets
    # (2) --> ttbar
    rTTBar = (
        integralDATA_dyjets * integralDYJets_ttbar
        - (integralMCothers_dyjets + integralQCD_dyjets) * integralDYJets_ttbar
    )
    rTTBar += (
        integralMCothers_ttbar + integralQCD_ttbar - integralDATA_ttbar
    ) * integralDYJets_dyjets
    try:
        rTTBar /= (
            integralTTbar_dyjets * integralDYJets_ttbar
            - integralTTbar_ttbar * integralDYJets_dyjets
        )
    except ZeroDivisionError:
        print("ERROR in rTTBar: ZeroDivisionError: one or more of the integrals is zero:")
        print(integralTTbar_dyjets, integralDYJets_ttbar, "-", integralTTbar_ttbar, integralDYJets_dyjets)
    rDYJets = (
        integralDATA_dyjets * integralTTbar_ttbar
        - (integralMCothers_dyjets + integralQCD_dyjets) * integralTTbar_ttbar
    )
    rDYJets += (
        integralMCothers_ttbar + integralQCD_ttbar - integralDATA_ttbar
    ) * integralTTbar_dyjets
    try:
        rDYJets /= (
            integralTTbar_ttbar * integralDYJets_dyjets
            - integralTTbar_dyjets * integralDYJets_ttbar
        )
    except ZeroDivisionError:
        print("ERROR in rWJets: ZeroDivisionError: one or more of the integrals is zero:")
        print(integralTTbar_ttbar, integralDYJets_dyjets, "-", integralTTbar_dyjets, integralDYJets_ttbar)

    # FIXME
    ##draw histo
    # self.histoMCall.SetFillColor(kBlue)
    # self.histoDATA.SetMarkerStyle(20)

    # self.histoMCall.Draw("HIST")
    # self.histoDATA.Draw("psame")
    # self.histoMCall.GetXaxis().SetRangeUser(self.xminplot,self.xmaxplot)
    # self.histoMCall.GetYaxis().SetRangeUser(self.yminplot,self.ymaxplot)

    # canvas.Update()
    # gPad.RedrawAxis()
    # gPad.Modified()
    ##canvas.SaveAs(self.name + ".eps","eps")
    ##canvas.SaveAs(self.name + ".pdf","pdf")
    # canvas.Print(fileps)
    # canvas.Print(self.name + ".C")
    ## make root file
    # tfile = TFile(self.name+'.root','recreate')
    # tfile.cd()
    # self.histoDATA.Write()
    # self.histoTTbar.Write()
    # self.histoMCall.Write()
    # self.histoQCD.Write()
    # self.histoZJet.Write()
    # self.histoWJet.Write()
    # self.histoSingleTop.Write()
    # self.histoDiboson.Write()
    # tfile.Close()

    # printout
    print()
    print(" TTBar ")
    print("######################################## ")
    print("name:                         " + plotObjTTBar.name)
    print("integral range:               " + str(
        plotObjTTBar.xmin
    ) + " < Mee < " + str(plotObjTTBar.xmax) + " GeV/c2")
    print("integral MC All:              " + str(integralMCall_ttbar) + " +/- " + str(
        ERRintegralMCall_ttbar
    ))
    if doQCD:
        print("integral QCD:                 " + str(integralQCD_ttbar) + " +/- " + str(
            ERRintegralQCD_ttbar
        ))
    print("integral MC TTbar:            " + str(integralTTbar_ttbar) + " +/- " + str(
        ERRintegralTTbar_ttbar
    ))
    print("integral MC DYJets:            " + str(integralDYJets_ttbar) + " +/- " + str(
        ERRintegralDYJets_ttbar
    ))
    print("integral MC other:            " + str(
        integralMCothers_ttbar
    ) + " +/- " + str(ERRintegralMCothers_ttbar))
    print("rescaled integral MC TTbar:   " + str(
        rTTBar * integralTTbar_ttbar
    ) + " +/- " + str(rTTBar * ERRintegralTTbar_ttbar))
    print("rescaled integral MC DYJets:   " + str(
        rDYJets * integralDYJets_ttbar
    ) + " +/- " + str(rDYJets * ERRintegralDYJets_ttbar))
    print("integral DATA:                " + str(integralDATA_ttbar) + " +/- " + str(
        ERRintegralDATA_ttbar
    ))
    print("contribution from other bkgs (except TTbar): " + str(
        contamination_ttbar * 100
    ) + "% [purity " + str(100-contamination_ttbar*100) + "%]")
    # print "integral DATA (corrected for contribution from other bkgs): "  + str( integralDATAcorr_ttbar ) + " +/- " + str( ERRintegralDATAcorr_ttbar )
    print("rescale factor for TTbar background: " + str(rTTBar) + " +/- " + str(
        rTTBarSigma
    ))
    # print "systematical uncertainty of TTbar background modeling: " + str(relERRrescale*100) + "%"
    print("######################################## ")
    print(" DYJets ")
    print("######################################## ")
    print("name:                         " + plotObjDYJets.name)
    print("integral range:               " + str(
        plotObjDYJets.xmin
    ) + " < Mee < " + str(plotObjDYJets.xmax) + " GeV/c2")
    print("integral MC All:              " + str(integralMCall_dyjets) + " +/- " + str(
        ERRintegralMCall_dyjets
    ))
    if doQCD:
        print("integral QCD:                 " + str(integralQCD_dyjets) + " +/- " + str(
            ERRintegralQCD_dyjets
        ))
    print("integral MC TTbar:            " + str(integralTTbar_dyjets) + " +/- " + str(
        ERRintegralTTbar_dyjets
    ))
    print("integral MC DYJets:            " + str(integralDYJets_dyjets) + " +/- " + str(
        ERRintegralDYJets_dyjets
    ))
    print("integral MC other:            " + str(
        integralMCothers_dyjets
    ) + " +/- " + str(ERRintegralMCothers_dyjets))
    print("rescaled integral MC TTbar:   " + str(
        rTTBar * integralTTbar_dyjets
    ) + " +/- " + str(rTTBar * ERRintegralTTbar_dyjets))
    print("rescaled integral MC DYJets:   " + str(
        rDYJets * integralDYJets_dyjets
    ) + " +/- " + str(rDYJets * ERRintegralDYJets_dyjets))
    print("integral DATA:                " + str(integralDATA_dyjets) + " +/- " + str(
        ERRintegralDATA_dyjets
    ))
    print("contribution from other bkgs (except wjets): " + str(
        contamination_dyjets * 100
    ) + "% [purity " + str(100-contamination_dyjets*100) + "%]")
    # print "integral DATA (corrected for contribution from other bkgs): "  + str( integralDATAcorr_wjets ) + " +/- " + str( ERRintegralDATAcorr_wjets )
    print("rescale factor for DYJets background: " + str(rDYJets) + " +/- " + str(
        rDYJetsSigma
    ))
    print("######################################## ")
    print()

    if plotObjTTBar.writeXSecFile:
        # create new cross section file -- ttbar
        originalFileName = plotObjTTBar.fileXsectionNoRescale.split("/")[-1].split(".")[0]
        ttbarFileName = originalFileName + "_" + plotObjTTBar.name + ".txt"
        os.system("rm -f " + ttbarFileName)
        outputFile = open(ttbarFileName, "w")

        for line in open(plotObjTTBar.fileXsectionNoRescale):
            line = line.strip("\n")
            lineNoComments = line.split("#")[
                0
            ]  # strip off anything after any '#' if present
            # ignore empty lines
            if len(lineNoComments) <= 0:
                print(line, file=outputFile)
                continue
            # FIXME this doesn't support keeping comments at the end of lines

            if re.search(plotObjTTBar.datasetName, line):
                list = re.split("\s+", line)
                newline = (
                    str(list[0]) + "    " + str("%.6f" % (float(list[1]) * float(rTTBar)))
                )
                print(newline, file=outputFile)
            else:
                print(line, file=outputFile)

        outputFile.close()
        print("New xsection file (after TTBar rescaling) is: " + ttbarFileName)
        print(" ")

    if plotObjDYJets.writeXSecFile:
        # NB: for the moment, we only support writing the DYJets SF when the TTBar SF has also been written
        originalFileName = ttbarFileName.split("/")[-1].split(".")[0]
        newFileName = originalFileName + "_" + plotObjDYJets.name + ".txt"
        os.system("rm -f " + newFileName)
        outputFile = open(newFileName, "w")

        for line in open(ttbarFileName):
            line = line.strip("\n")
            lineNoComments = line.split("#")[
                0
            ]  # strip off anything after any '#' if present
            # ignore empty lines
            if len(lineNoComments) <= 0:
                print(line, file=outputFile)
                continue

            if re.search(plotObjDYJets.datasetName, line):
                list = re.split("\s+", line)
                newline = (
                    str(list[0]) + "    " + str("%.6f" % (float(list[1]) * float(rDYJets)))
                )
                print(newline, file=outputFile)
            else:
                print(line, file=outputFile)

        outputFile.close()
        print("New xsection file (after DYJets rescaling) is: " + newFileName)
        print(" ")
    return rDYJets, rDYJetsSigma, rTTBar, rTTBarSigma


# ############ DON'T NEED TO MODIFY ANYTHING HERE - END #######################
##############################################################################
zjetSamples = ["ZJet_amcatnlo_ptBinned_IncStitch", "ZJet_madgraphLO_HT", "ZJet_powhegminnlo", "ZJet_amcatnlo_Inc"]


##############################################################################
# ############ USER CODE - BEGIN ##############################################
doQCD = True
useSingleFakeMC = False
digitsForTable = 3
if not doQCD:
    print("INFO: ignoring QCD")
if useSingleFakeMC and doQCD:
    raise RuntimeError("Can't specify QCD from data as well as single fakes from MC")
if len(sys.argv) < 4:
    print("ERROR: did not find MC/data combined plot file or QCD plot file or year")
    print("Usage: python calc_DYJetsAndTTBarRescale_And_xsecFile.py combinedQCDPlotFile.root combinedDataMCPlotFile.root year")
    exit(-1)
if len(sys.argv) > 5:
    print("ERROR: found extra arguments")
    print("Usage: python calc_DYJetsAndTTBarRescale_And_xsecFile.py combinedQCDPlotFile.root combinedDataMCPlotFile.root year")
    exit(-1)

qcdFile = sys.argv[1]
mcFile = sys.argv[2]
year = sys.argv[3]
if len(sys.argv) > 4:
    zjet = sys.argv[4]
    if not any(sample in zjet for sample in zjetSamples):
        raise RuntimeError("zjet sample '{}' found as command-line argument is not known by the script. Possibilities are: {}.".format(zjet, zjetSamples))
else:
    zjet = "ZJet_amcatnlo_ptBinned_IncStitch"
    # zjet = "ZJet_madgraphLO_HT"
    # zjet = "ZJet_powhegminnlo"
    #zjet = "ZJet_amcatnlo_Inc"
    print("INFO: defaulting to zjet sample '{}' as none was given as a command-line argument".format(zjet))
zjetDatasetName = "DYJetsToEE.+" if zjet == "ZJet_powhegminnlo" else "DYJetsToLL.+"

# --- Input files
if doQCD:
    File_QCD_preselection = GetFile(qcdFile)
File_preselection = GetFile(mcFile)

xsectionFile = "/afs/cern.ch/user/s/scooper/work/private/LQNanoAODAttempt/Leptoquarks/analyzer/rootNtupleAnalyzerV2/config/xsection_13TeV_2022.txt"

histNameDefault = "histo1D__SAMPLE__"
hist2DNameDefault = "histo2D__SAMPLE__"
# histNameReleaseMee = "histo1D__SAMPLE__cutHisto_allOtherCuts___________"
histNameAllPrevCuts = "histo1D__SAMPLE__cutHisto_allPreviousCuts________"
#
# samples to use
# wjet = "WJet_amcatnlo_Inc"
# wjet = "WJet_amcatnlo_ptBinned"
wjet = "WJet_amcatnlo_jetBinned"
# wjet = "WJet_amcatnlo_ptBinned"
allBkg = "allBkg"
data = "DATA"
if not doQCD and useSingleFakeMC:
    ttbar = "TTbar_powheg_all"
    diboson = "DIBOSON_nlo_all"
else:
    ttbar = "TTTo2L2Nu"
    diboson = "DIBOSON_nlo"
ttbarDatasetName = "TTT"
# ttbar = "TTbar_powheg_all"
# diboson = "DIBOSON_nlo_all"
# ttbarDatasetName = ttbar
singletop = "SingleTop"
qcd = "QCDFakes_DATA"

# --- Rescaling of DY+jets and ttbar+jets backgrounds
histBaseNames = []
# nominal
histBaseNames.append("Mee_BkgControlRegion")
histBaseNames.append("Mee_BkgControlRegion_BJETBIN1")
histBaseNames.append("Mee_BkgControlRegion_BJETBIN2")
histBaseNames.append("Mee_EBEB_BkgControlRegion")
histBaseNames.append("Mee_EBEE_BkgControlRegion")
histBaseNames.append("Mee_EEEE_BkgControlRegion")
#histBaseNames.append("Mee_70_110_LQ300")
#histBaseNames.append("Mee_70_110_LQ600")
#histBaseNames.append("Mee_70_110_LQ800")
#histBaseNames.append("Mee_70_110_LQ900")
# SF variations
histBaseNames.append("MeeVsNJet_BkgControlRegion_BJETBIN2")
histBaseNames.append("MeeVsST_BkgControlRegion_BJETBIN2")
histBaseNames.append("MeeVsSTjet_BkgControlRegion_BJETBIN2")
histBaseNames.append("MeeVsSTlep_BkgControlRegion_BJETBIN2")
histBaseNames.append("MeeVsMejMin_BkgControlRegion_BJETBIN2")
histBaseNames.append("MeeVsMejMax_BkgControlRegion_BJETBIN2")
histBaseNames.append("MeeVsMeejj_BkgControlRegion_BJETBIN2")
histBaseNames.append("MeeVsMe1j1_BkgControlRegion_BJETBIN2")
histBaseNames.append("MeeVsMe1j2_BkgControlRegion_BJETBIN2")
histBaseNames.append("MeeVsMe2j1_BkgControlRegion_BJETBIN2")
histBaseNames.append("MeeVsMe2j2_BkgControlRegion_BJETBIN2")
histBaseNames.append("MeeVsEle1Pt_BkgControlRegion_BJETBIN2")
histBaseNames.append("MeeVsEle2Pt_BkgControlRegion_BJETBIN2")
histBaseNames.append("MeeVsJet1Pt_BkgControlRegion_BJETBIN2")
histBaseNames.append("MeeVsJet2Pt_BkgControlRegion_BJETBIN2")
histBaseNames.append("MeeVsPFMETType1Pt_BkgControlRegion_BJETBIN2")

STBins = [300, 400, 500, 700, 900, 1250, "Inf"]
STLepBins = [300, 400, 500, 700, "Inf"]
MejBins = [100, 200, 300, 400, 500, "Inf"]
MejMinBins = [100, 200, 300, "Inf"]
MejMaxBins = [100, 200, 300, 400, "Inf"]
# PtBins = [100, 200, 300, 400, 500, "Inf"]
PtBins = [50, 100, 200, "Inf"]
Pt2Bins = [50, 100, "Inf"]
METBins =  [0, 50, 100, 200, "Inf"]
binsForVarDict = {}
#binsForVarDict["NJet"] = [1.5, 2.5, 3.5, 4.5, 5.5, "Inf"]
binsForVarDict["NJet"] = [1.5, 2.5, 3.5, 4.5, "Inf"]
binsForVarDict["ST"] = STBins
binsForVarDict["STjet"] = STBins
binsForVarDict["STlep"] = STLepBins
binsForVarDict["MejMin"] = MejMinBins
binsForVarDict["MejMax"] = MejMaxBins
binsForVarDict["Meejj"] = MejBins
binsForVarDict["Me1j1"] = MejBins
binsForVarDict["Me1j2"] = MejBins
binsForVarDict["Me2j1"] = MejBins
binsForVarDict["Me2j2"] = MejBins
binsForVarDict["Ele1Pt"] = PtBins
binsForVarDict["Ele2Pt"] = Pt2Bins
binsForVarDict["Jet1Pt"] = PtBins
binsForVarDict["Jet2Pt"] = PtBins
binsForVarDict["PFMETType1Pt"] = METBins

varToXTitleDict = {}
varToXTitleDict["NJet"] = "N_{jet}"
varToXTitleDict["ST"] = "S_{T} [GeV]"
varToXTitleDict["STjet"] = "S_{T}(jets) [GeV]"
varToXTitleDict["STlep"] = "S_{T}(electrons) [GeV]"
varToXTitleDict["MejMin"] = "M_{ej}^{min} [GeV]"
varToXTitleDict["MejMax"] = "M_{ej}^{max} [GeV]"
varToXTitleDict["Meejj"] = "M_{eejj} [GeV]"
varToXTitleDict["Me1j1"] = "M_{e1j1} [GeV]"
varToXTitleDict["Me1j2"] = "M_{e1j2} [GeV]"
varToXTitleDict["Me2j1"] = "M_{e2j1} [GeV]"
varToXTitleDict["Me2j2"] = "M_{e2j2} [GeV]"
varToXTitleDict["Ele1Pt"] = "P_{T}(e_{1}) [GeV]"
varToXTitleDict["Ele2Pt"] = "P_{T}(e_{2}) [GeV]"
varToXTitleDict["Jet1Pt"] = "P_{T}(j_{1}) [GeV]"
varToXTitleDict["Jet2Pt"] = "P_{T}(j_{2}) [GeV]"
varToXTitleDict["PFMETType1Pt"] = "PFMET [GeV]"

varToYRangeDict = {}
varToYRangeDict["NJet"] = [0.5, 4]
varToYRangeDict["ST"] = [0.5, 1.5]
varToYRangeDict["STjet"] = [0.5, 1.5]
varToYRangeDict["STlep"] = [0.5, 1.5]
varToYRangeDict["MejMin"] = [0.5, 1.5]
varToYRangeDict["MejMax"] = [0.5, 1.5]
varToYRangeDict["Meejj"] = [0.5, 1.5]
varToYRangeDict["Me1j1"] = [0.5, 1.5]
varToYRangeDict["Me1j2"] = [0.5, 1.5]
varToYRangeDict["Me2j1"] = [0.5, 1.5]
varToYRangeDict["Me2j2"] = [0.5, 1.5]
varToYRangeDict["Ele1Pt"] = [0.5, 1.5]
varToYRangeDict["Ele2Pt"] = [0.5, 1.5]
varToYRangeDict["Jet1Pt"] = [0.5, 1.5]
varToYRangeDict["Jet2Pt"] = [0.5, 1.5]
varToYRangeDict["PFMETType1Pt"] = [0, 8]

legCoordsDict = {}
legCoordsDict["upper-left"] = [0.12,0.65,0.5,0.85]
legCoordsDict["upper-left-thin"] = [0.12,0.65,0.495,0.85]
legCoordsDict["upper-right"] = [0.5,0.65,0.88,0.85]
legCoordsDict["lower-right"] = [0.5,0.15,0.88,0.35]
legCoordsDict["middle-right"] = [0.5,0.25,0.88,0.55]

varToLegPosDict = { "DYJets": {}, "TTBar": {}}
varToLegPosDict["DYJets"]["STlep"] = "upper-right"
varToLegPosDict["DYJets"]["NJet"] = "upper-left-thin"
varToLegPosDict["DYJets"]["Ele1Pt"] = "upper-right"
varToLegPosDict["DYJets"]["Ele2Pt"] = "upper-right"
varToLegPosDict["DYJets"]["MejMin"] = "lower-right"
varToLegPosDict["DYJets"]["MejMax"] = "lower-right"
varToLegPosDict["DYJets"]["PFMETType1Pt"] = "middle-right"
varToLegPosDict["TTBar"]["ST"] = "upper-right"
varToLegPosDict["TTBar"]["STjet"] = "upper-left"
varToLegPosDict["TTBar"]["STlep"] = "upper-left"
varToLegPosDict["TTBar"]["NJet"] = "upper-left-thin"
varToLegPosDict["TTBar"]["Ele1Pt"] = "upper-left"
varToLegPosDict["TTBar"]["Ele2Pt"] = "upper-right"
varToLegPosDict["TTBar"]["Jet1Pt"] = "upper-right"
varToLegPosDict["TTBar"]["Jet2Pt"] = "upper-right"
varToLegPosDict["TTBar"]["MejMin"] = "upper-right"
varToLegPosDict["TTBar"]["MejMax"] = "upper-right"
varToLegPosDict["TTBar"]["Me1j1"] = "upper-right"
varToLegPosDict["TTBar"]["Me1j2"] = "upper-right"
varToLegPosDict["TTBar"]["Me2j1"] = "upper-right"
varToLegPosDict["TTBar"]["Me2j2"] = "upper-right"
varToLegPosDict["TTBar"]["PFMETType1Pt"] = "upper-right"

histNameBaseForXSecFile = "Mee_BkgControlRegion_BJETBIN2"
nominalBkgSFPlotNames = ['Mee_BkgControlRegion_DYJets', 'Mee_BkgControlRegion_gteTwoBtaggedJets_TTbar']

meeMinTTBar = 140
meeMaxTTBar = 220 # TODO check range
meeMinDYJets = 80
meeMaxDYJets = 100
meePlotMaxDYJets = 2000
meePlotMaxTTBar = 2000

plotsTTBar = []
plotsDYJets = []
plotsNotFound = []

for idx, histBaseName in enumerate(histBaseNames):
    histName = histNameDefault
    writeXSecFile = (histBaseName == histNameBaseForXSecFile)

    bins = []
    if "Vs" in histBaseName:
        histName = hist2DNameDefault
        var = histBaseName[histBaseName.find("Vs")+2:histBaseName.find("_", histBaseName.find("Vs"))]
        bins = binsForVarDict[var]

    try:
        plotsDYjets, plotsTTbar = SetupPlots(histName, histBaseName, bins)
        for plot in plotsDYjets:
            plot.writeXSecFile = writeXSecFile
        for plot in plotsTTbar:
            plot.writeXSecFile = writeXSecFile
        plotsTTBar.extend(plotsTTbar)
        plotsDYJets.extend(plotsDYjets)
    except RuntimeError as e:
        print("Caught exception while getting histo: '", e, "'; skipping this one")
        plotsNotFound.append(histBaseName)
        continue

# -----------------------------------------------------------------------------------


# ############ USER CODE - END ################################################
##############################################################################


# --- Generate and print the plots from the list 'plots' defined above

# --- Output files
fileps = "allPlots_calc_dyJetsAndTTBarRescale_And_xsecFile.ps"

# --- Generate and print the plots from the list 'plots' define above
plotErrors = []
nominalScaleFactors = {}
nominalScaleFactorErrs = {}
sfXValsForVar = {}
sfXErrsForVar = {}
sfYValsForVar = {}
sfYErrsForVar = {}
sfXValsForVar["DYJets"] = {}
sfXErrsForVar["DYJets"] = {}
sfYValsForVar["DYJets"] = {}
sfYErrsForVar["DYJets"] = {}
sfXValsForVar["TTBar"] = {}
sfXErrsForVar["TTBar"] = {}
sfYValsForVar["TTBar"] = {}
sfYErrsForVar["TTBar"] = {}
for idx, plot in enumerate(plotsTTBar):
    dyjPlot = plotsDYJets[idx]
    print("INFO: Calculating scale factors for plot={}".format(dyjPlot.name), flush=True)
    try:
        rDYJets, rDYJetsSigma, rTTBar, rTTBarSigma = CalculateRescaleFactor(plot, dyjPlot, fileps)
    except Exception as e:
        plotErrors.append("Had an exception doing CalculateRescaleFactor for plot {}: {}".format(plot.name, e))
        print("WARN: Had an exception doing CalculateRescaleFactor for plot {}: {}".format(plot.name, e), flush=True)
        continue
    if dyjPlot.name == nominalBkgSFPlotNames[0]:
        nominalScaleFactors["DYJets"] = rDYJets
        nominalScaleFactorErrs["DYJets"] = rDYJetsSigma
    if plot.name == nominalBkgSFPlotNames[1]:
        nominalScaleFactors["TTBar"] = rTTBar
        nominalScaleFactorErrs["TTBar"] = rTTBarSigma
    if "Vs" in plot.name or "Vs" in dyjPlot.name:
        var = plot.name[plot.name.find("Vs")+2:plot.name.find("_", plot.name.find("Vs"))]
        xrange = plot.name.split("_")[-2]
        xmin = float(xrange.split("to")[0])
        xmax = xrange.split("to")[1]
        if xmax == "Inf":
            xmax = plot.xmaxplot
        else:
            xmax = float(xmax)
        xwidth = xmax-xmin
        xcenter = (xmax+xmin)/2
        for bkgName in ["DYJets", "TTBar"]:
            if var not in sfXValsForVar[bkgName].keys():
                sfXValsForVar[bkgName][var] = []
                sfXErrsForVar[bkgName][var] = []
                sfYValsForVar[bkgName][var] = []
                sfYErrsForVar[bkgName][var] = []
            sfVal = rDYJets if bkgName == "DYJets" else rTTBar
            sfErr = rDYJetsSigma if bkgName == "DYJets" else rTTBarSigma
            sfXValsForVar[bkgName][var].append(xcenter)
            sfXErrsForVar[bkgName][var].append(xwidth/2)
            sfYValsForVar[bkgName][var].append(sfVal)
            sfYErrsForVar[bkgName][var].append(sfErr)

for bkgName in ["DYJets", "TTBar"]:
    for var in sfXValsForVar[bkgName].keys():
        xPoints = sfXValsForVar[bkgName][var]
        xPointErrs = sfXErrsForVar[bkgName][var]
        yPoints = sfYValsForVar[bkgName][var]
        yPointErrs = sfYErrsForVar[bkgName][var]
        sfNom = nominalScaleFactors[bkgName]
        sfNomErr = nominalScaleFactorErrs[bkgName]
        canvas = TCanvas()
        canvas.cd()
        canvas.SetGridy()
        graph = TGraphErrors(len(xPoints),np.array(xPoints, dtype="f"),np.array(yPoints, dtype="f"),np.array(xPointErrs, dtype="f"),np.array(yPointErrs, dtype="f"))
        graph.SetTitle('')
        graph.SetMarkerColor(kBlue)
        graph.SetLineColor(kBlue)
        graph.Draw('ap0')
        graph.GetXaxis().SetTitle(varToXTitleDict[var])
        graph.GetXaxis().SetTitleSize(0.05)
        graph.GetXaxis().SetTitleOffset(0.8)
        graph.GetYaxis().SetTitle("data/MC")
        graph.GetYaxis().SetTitleSize(0.05)
        graph.GetYaxis().SetTitleOffset(0.8)
        uncertaintyXpoints = [xPoints[0]-xPointErrs[0],xPoints[-1]+xPointErrs[-1]]
        uncertaintyYpoints = [sfNom,sfNom]
        uncertaintyXpointsErrs = [0,0]
        uncertaintyYpointsErrs = [sfNomErr,sfNomErr]
        uncertaintyRegionGraph = TGraphErrors(len(uncertaintyXpoints),np.array(uncertaintyXpoints, dtype="f"),np.array(uncertaintyYpoints, dtype="f"),np.array(uncertaintyXpointsErrs, dtype="f"),np.array(uncertaintyYpointsErrs, dtype="f"))
        uncertaintyRegionGraph.SetFillColor(15)
        uncertaintyRegionGraph.SetFillStyle(3001)
        uncertaintyRegionGraph.SetLineColor(1)
        uncertaintyRegionGraph.Draw('3l')
        legPos = "upper-left"
        if var in varToLegPosDict[bkgName].keys():
            legPos = varToLegPosDict[bkgName][var]
        legCoords = legCoordsDict[legPos]
        leg = TLegend(*legCoords)
        leg.SetTextSize(0.0265)
        graph.GetYaxis().SetRangeUser(varToYRangeDict[var][0], varToYRangeDict[var][1])
        leg.AddEntry(graph,'{} scale factor variations'.format(bkgName),'lp')
        leg.AddEntry(uncertaintyRegionGraph,'Nominal scale factor = '+str(round(sfNom,3))+' #pm '+str(round(sfNomErr,3)),'fl')
        leg.Draw()
        graph.Draw('p0')
        canvas.Modified()
        baseName = bkgName+'_scaleFactorVariation_{}Bins'.format(var)
        # suppress INFO messages when saving canvases
        prevLevel = ROOT.gErrorIgnoreLevel
        ROOT.gErrorIgnoreLevel = ROOT.kWarning
        canvas.Print(baseName+'.png')
        canvas.Print(baseName+'.pdf')
        ROOT.gErrorIgnoreLevel = prevLevel
        if var == "NJet":
            ratioSub = "R_Z" if bkgName == "DYJets" else "R_ttbar"
            ratioSubLtx = "$R_{Z}$" if bkgName == "DYJets" else "$R_{t\\bar{t}}$"
            titlesText = ["year", ratioSub, "jets"]
            titlesLatex = ["year", ratioSubLtx, "jets"]
            tableText = []
            tableLatex = []
            for idx, jetBin in enumerate(xPoints):
                if idx == len(xPoints)-1:
                    jetsText = ">= {}".format(xPoints[idx-1]+1)
                    jetsLatex = "$\geq$~{}".format(xPoints[idx-1]+1)
                else:
                    jetsText = jetBin
                    jetsLatex = jetBin
                tableText.append([year, "{:.6f} +/- {:.6f}".format(yPoints[idx], yPointErrs[idx]), jetsText])
                tableLatex.append([year, "{:.2f} +/- {:.2f}".format(yPoints[idx], yPointErrs[idx]), jetsLatex])
            print("-"*50)
            print(bkgName)
            print(tabulate(tableText, headers=titlesText, tablefmt="fancy_grid"))
            print()
            print(tabulate(tableLatex, headers=titlesLatex, tablefmt="latex_raw"))
            print("-"*50)

print("For latex scale factor table:")
tableLineStr = "{:12} & {}~$\pm$~{} & {}~$\pm$~{}\\\\".format(year, round(nominalScaleFactors["DYJets"], digitsForTable), round(nominalScaleFactorErrs["DYJets"], digitsForTable), round(nominalScaleFactors["TTBar"], digitsForTable), round(nominalScaleFactorErrs["TTBar"], digitsForTable))
print(tableLineStr)
print()

if len(plotsNotFound):
    print("WARN: plots not found:", plotsNotFound, file=sys.stderr)
if len(plotErrors):
    for err in plotErrors:
        print("WARN: {}".format(err), file=sys.stderr)
print("INFO: year = {}".format(year))
print("INFO: using file: " + File_preselection.GetName())
if doQCD:
    print("INFO: using QCD file: " + File_QCD_preselection.GetName())
print("INFO: using samples:")
print("\t DATA ------>", data)
print("\t DY ---------->", zjet, "; datasetname =", zjetDatasetName)
if not doQCD:
    print("\t W ----------->", wjet)
print("\t ttbar ------->", ttbar, "; datasetname =", ttbarDatasetName)
print("\t diboson ----->", diboson)
if doQCD:
    print("\t QCD --------->", qcd)
print("\t SingleTop --->", singletop)
