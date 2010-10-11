#!/usr/bin/env python

##############################################################################
## USER CODE IS TOWARD THE END OF THE FILE
##############################################################################

##############################################################################
############# DON'T NEED TO MODIFY ANYTHING HERE - BEGIN #####################

#---Import
import sys
import string
from optparse import OptionParser
import os.path
from ROOT import *
import re
import ROOT
from array import array

#--- ROOT general options
gROOT.SetBatch(kTRUE);
gStyle.SetOptStat(0)
gStyle.SetPalette(1)
gStyle.SetCanvasBorderMode(0)
gStyle.SetFrameBorderMode(0)
gStyle.SetCanvasColor(kWhite)
gStyle.SetPadTickX(1);
gStyle.SetPadTickY(1);
#--- TODO: WHY IT DOES NOT LOAD THE DEFAULT ROOTLOGON.C ? ---#

def GetFile(filename):
    file = TFile(filename)
    if( not file):
        print "ERROR: file " + filename + " not found"
        print "exiting..."
        sys.exit()
    return file


def GetHisto( histoName , file , scale = 1 ):
    file.cd()
    histo = file.Get( histoName )
    if(scale!=1):
        histo.Scale(scale)
    if( not histo):
        print "ERROR: histo " + histoName + " not found in " + file.GetName()
        print "exiting..."
        sys.exit()
    return histo

def generateHistoList( histoBaseName , samples , variableName, fileName ):
    histolist = []
    for sample in samples:
        hname = (histoBaseName.replace("SAMPLE", sample)).replace("VARIABLE", variableName)
        histolist.append(GetHisto(hname, fileName).Clone())
    return histolist
                                                    
def generateAndAddHistoList( histoBaseName , samples , variableNames, fileName ):
    histolist = []
    for sample in samples:
        iv=0
        for variableName in variableNames:
            hname = (histoBaseName.replace("SAMPLE", sample)).replace("VARIABLE", variableName)
            if (iv==0):
                histo = GetHisto(hname, fileName).Clone()
            else:
                histo.Add(GetHisto(hname, fileName).Clone())
            iv=iv+1
        histolist.append(histo)
    return histolist

def generateHisto( histoBaseName , sample , variableName, fileName ):
    hname = (histoBaseName.replace("SAMPLE", sample)).replace("VARIABLE", variableName)
    histo = GetHisto(hname, fileName).Clone()
    return histo

def generateAndAddHisto( histoBaseName , sample , variableNames, fileName ):
    iv=0
    for variableName in variableNames:
        hname = (histoBaseName.replace("SAMPLE", sample)).replace("VARIABLE", variableName)
        if (iv==0):
            histo = GetHisto(hname, fileName).Clone()
        else:
            histo.Add(GetHisto(hname, fileName).Clone())
        iv=iv+1
    return histo

                                                    
## The Plot class: add members if needed
class Plot:
    histos      = [] # list of histograms to be plotted in this plot
    keys        = [] # list of keys to be put in the legend (1 key per histo)
    histosStack = [] # list of histograms to be plotted in this plot -- stack histo format
    keysStack   = [] # list of keys to be put in the legend (1 key per histo) -- stack histo format
    xtit        = "" # xtitle
    ytit        = "" # ytitle
    xmin        = "" # min x axis range (need to set both min and max. Leave it as is for full range)
    xmax        = "" # max x axis range (need to set both min and max. Leave it as is for full range)
    ymin        = "" # min y axis range (need to set both min and max. Leave it as is for full range)
    ymax        = "" # max y axis range (need to set both min and max. Leave it as is for full range)
    lpos        = "" # legend position (default = top-right, option="bottom-center", "top-left")
    #    xlog        = "" # log scale of X axis (default = no, option="yes") ### IT SEEMS IT DOES NOT WORK
    ylog        = "" # log scale of Y axis (default = no, option="yes")
    rebin       = "" # rebin x axis (default = 1, option = set it to whatever you want )
    name        = "" # name of the final plots
    lint        = "2.9 pb^{-1}" # integrated luminosity of the sample ( example "10 pb^{-1}" )
    addZUncBand = "no" # add an uncertainty band coming from the data-MC Z+jets rescaling (default = "no", option="yes")
    ZUncKey     = "Z/#gamma/Z* + jets unc." # key to be put in the legend for the Z+jets uncertainty band
    ZPlotIndex  = 1 # index of the Z+jets plots in the histosStack list (default = 1)
    ZScaleUnc   = 0.20 # uncertainty of the data-MC Z+jets scale factor
    histodata   = "" # data histogram

    def Draw(self, fileps):

        #-- create canvas
        canvas = TCanvas()
        stack = {}

        #-- log scale
#             xlog may npot work         if (plot.xlog     == "yes"):
#             canvas.SetLogx();
        if (plot.ylog     == "yes"):
            canvas.SetLogy();

        #-- legend
        hsize=0.20
        vsize=0.25
        if (plot.lpos=="bottom-center"):
            xstart=0.35
            ystart=0.25
        elif(plot.lpos=="top-left"):
            xstart=0.12
            ystart=0.63
        else:
            xstart=0.68
            ystart=0.63
        legend = TLegend(xstart, ystart, xstart+hsize, ystart+vsize)
        legend.SetFillColor(kWhite)
        legend.SetMargin(0.2)

        #-- loop over histograms (stacked)
        Nstacked = len(plot.histosStack)
        for iter in range(0, Nstacked):
            #make this stack
            stack[iter] = TH1F()
            Nloop = Nstacked - iter
            for iter1 in range(0,Nloop):
                histo = plot.histosStack[iter1]
                if(iter1==0):
                    stack[iter].SetBins( histo.GetNbinsX(), histo.GetXaxis().GetXmin(), histo.GetXaxis().GetXmax() )
                    #stack[iter].SetName( plot.keysStack[iter] )
                stack[iter].Add(histo)
            #define style
            if(plot.rebin!=""):
                stack[iter].Rebin(plot.rebin)
            stack[iter].SetMarkerStyle(20+2*iter)
            stack[iter].SetMarkerColor(15+10*iter)
            stack[iter].SetLineColor(  15+10*iter)
            stack[iter].SetFillColor(  15+10*iter)
            legend.AddEntry(stack[iter], plot.keysStack[Nstacked - iter - 1],"lf")
            #draw stack
            if iter==0:
                thisMin = stack[iter].GetXaxis().GetXmin()
                thisMax = stack[iter].GetXaxis().GetXmax()
                thisNbins = stack[iter].GetNbinsX()
                newBinning = (thisMax - thisMin) / thisNbins
                stack[iter].SetTitle("")
                stack[iter].GetXaxis().SetTitle(plot.xtit)
                stack[iter].GetYaxis().SetTitle(plot.ytit + " / ( "+ str(newBinning) + " )")
                if (plot.xmin!="" and plot.xmax!=""):
                    stack[iter].GetXaxis().SetRangeUser(plot.xmin,plot.xmax)
                if (plot.ymin!="" and plot.ymax!=""):
                    stack[iter].GetYaxis().SetLimits(plot.ymin,plot.ymax)
                    stack[iter].GetYaxis().SetRangeUser(plot.ymin,plot.ymax)
                #search for maximum of histograms
                #maxHisto = stack[iter].GetMaximum()
                #print maxHisto
                #for hh in plot.histos:
                #    if(plot.rebin!=""):
                #        if(hh.GetMaximum()*plot.rebin > maxHisto):
                #            maxHisto = hh.GetMaximum()*plot.rebin
                #    else:
                #        if(hh.GetMaximum() > maxHisto):
                #            maxHisto = hh.GetMaximum()
                #stack[iter].GetYaxis().SetLimits(0.,maxHisto*1.2)
                #stack[iter].GetYaxis().SetRangeUser(0.001,maxHisto*1.2)
                #draw first histo
                stack[iter].Draw("HIST")
            else:
                stack[iter].Draw("HISTsame")

        #-- Z+jets uncertainty band
        if(plot.addZUncBand == "yes"):
            Zhisto = plot.histosStack[plot.ZPlotIndex].Clone()
            if(plot.rebin!=""):
                Zhisto.Rebin(plot.rebin)
            zUncHisto = stack[0].Clone()
            for bin in range(0,Zhisto.GetNbinsX()):
              zUncHisto.SetBinError(bin+1,plot.ZScaleUnc*Zhisto.GetBinContent(bin+1))
            zUncHisto.SetMarkerStyle(0)
            zUncHisto.SetLineColor(0)
            zUncHisto.SetFillColor(5)
            zUncHisto.SetFillStyle(3154)
            zUncHisto.Draw("E2same")
            legend.AddEntry(zUncHisto, plot.ZUncKey,"f")

        #-- loop over histograms (overlaid)
        ih=0 # index of histo within a plot
        for histo in plot.histos:
            if(plot.rebin!=""):
                histo.Rebin(plot.rebin)
            histo.SetMarkerStyle(20+2*ih)
            histo.SetMarkerColor(2+2*ih)
            histo.SetLineColor(  2+2*ih)
            legend.AddEntry(histo, plot.keys[ih],"l")
            histo.Draw("HISTsame")
            ih=ih+1

        #-- plot data
        if(plot.histodata!=""):
            if(plot.rebin!=""):
                plot.histodata.Rebin(plot.rebin)
            plot.histodata.SetMarkerStyle(20)
            legend.AddEntry(plot.histodata, "data","p")
            plot.histodata.Draw("psame")

        #-- draw label
        l = TLatex()
        l.SetTextAlign(12)
        l.SetTextSize(0.04)
        l.SetTextFont(62)
        l.SetNDC()
#        l.DrawLatex(xstart,ystart-0.05,"CMS Preliminary 2010")
#        l.DrawLatex(xstart,ystart-0.10,"L_{int} = " + plot.lint)
        if (plot.lpos=="bottom-center"):
            l.DrawLatex(0.35,0.20,"CMS Preliminary 2010")
            l.DrawLatex(0.35,0.15,"L_{int} = " + plot.lint)
        if (plot.lpos=="top-left"):
            l.DrawLatex(xstart+hsize+0.02,ystart+vsize-0.03,"CMS Preliminary 2010")
            l.DrawLatex(xstart+hsize+0.02,ystart+vsize-0.08,"L_{int} = " + plot.lint)
        else:
            l.DrawLatex(xstart-hsize-0.10,ystart+vsize-0.03,"CMS Preliminary 2010")
            l.DrawLatex(xstart-hsize-0.10,ystart+vsize-0.08,"L_{int} = " + plot.lint)

        #-- end
        legend.Draw()
        canvas.Update()
        gPad.RedrawAxis()
        gPad.Modified()
        canvas.SaveAs(plot.name + ".eps","eps")
        canvas.SaveAs(plot.name + ".pdf","pdf")
        canvas.Print(fileps)



############# DON'T NEED TO MODIFY ANYTHING HERE - END #######################
##############################################################################


##############################################################################
############# USER CODE - BEGING #############################################

#--- Input root file

File_preselection = GetFile("$LQDATA/collisions/38X_2.9pb-1/output_cutTable_eejjSample_elePt25_jetPt20_withDEtaInEE/analysisClass_eejjSample_plots.root")

File_selection    = File_preselection

#### Common values for plots:
#otherBkgsKey="QCD, single top, VV+jets, W/W*+jets"
otherBkgsKey="Other Bkgs"
zUncBand="no"

pt_xmin=0
pt_xmax=400
pt_ymin=0.001
pt_ymax=300

eta_rebin=10
eta_ymin=0
eta_ymax=30



#--- Final plots are defined here

# Simply add or remove plots and update the list plots = [plot0, plot1, ...] below

histoBaseName = "histo1D__SAMPLE__cutHisto_allPreviousCuts________VARIABLE"

samplesForStackHistos = ["TTbar_Madgraph","ZJetAlpgen","OTHERBKG"]
keysStack =             ["ttbar", "Z/#gamma/Z* + jets", otherBkgsKey]

samplesForHistos = ["LQeejj_M200", "LQeejj_M300"]
keys             = ["LQ eejj M200","LQ eejj M300"]

sampleForDataHisto = "DATA"


#--- Mee_TwoEleOnly ---
variableName = "Mee_TwoEleOnly"

# h_Mee_LQeejj_M100 = GetHisto("histo1D__LQeejj_M100__cutHisto_allPreviousCuts________Mee_TwoEleOnly", File_preselection).Clone()
# h_Mee_LQeejj_M200 = GetHisto("histo1D__LQeejj_M200__cutHisto_allPreviousCuts________Mee_TwoEleOnly", File_preselection).Clone()
# h_Mee_LQeejj_M300 = GetHisto("histo1D__LQeejj_M300__cutHisto_allPreviousCuts________Mee_TwoEleOnly", File_preselection).Clone()
# h_Mee_LQeejj_M400 = GetHisto("histo1D__LQeejj_M400__cutHisto_allPreviousCuts________Mee_TwoEleOnly", File_preselection).Clone()
# h_Mee_LQeejj_M500 = GetHisto("histo1D__LQeejj_M500__cutHisto_allPreviousCuts________Mee_TwoEleOnly", File_preselection).Clone()
# h_Mee_TTbar = GetHisto("histo1D__TTbar_Madgraph__cutHisto_allPreviousCuts________Mee_TwoEleOnly", File_preselection).Clone()
# h_Mee_ZJetAlpgen = GetHisto("histo1D__ZJetAlpgen__cutHisto_allPreviousCuts________Mee_TwoEleOnly", File_preselection).Clone()
# h_Mee_OTHERBKG = GetHisto("histo1D__OTHERBKG__cutHisto_allPreviousCuts________Mee_TwoEleOnly", File_preselection).Clone()
# #h_Mee_QCD_Madgraph = GetHisto("histo1D__QCD_Madgraph__cutHisto_allPreviousCuts________Mee_TwoEleOnly", File_preselection).Clone()
# #h_Mee_QCDPt15 = GetHisto("histo1D__QCDPt15__cutHisto_allPreviousCuts________Mee_TwoEleOnly", File_preselection).Clone()
# h_Mee_SingleTop = GetHisto("histo1D__SingleTop__cutHisto_allPreviousCuts________Mee_TwoEleOnly", File_preselection).Clone()
# h_Mee_VVjets = GetHisto("histo1D__VVjets__cutHisto_allPreviousCuts________Mee_TwoEleOnly", File_preselection).Clone()
# h_Mee_WJetAlpgen = GetHisto("histo1D__WJetAlpgen__cutHisto_allPreviousCuts________Mee_TwoEleOnly", File_preselection).Clone()
h_Mee_DATA = GetHisto("histo1D__DATA__cutHisto_allPreviousCuts________Mee_TwoEleOnly", File_preselection).Clone()

plot0 = Plot()
## inputs for stacked histograms
plot0.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot0.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot0.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot0.keys            = keys
plot0.xtit            = "M(ee) (GeV/c^{2})"
plot0.ytit            = "Number of events"
# plot0.ylog            = "yes"
# plot0.rebin           = 1
# plot0.ymin            = 0.00000001
# plot0.ymax            = 20
plot0.ylog            = "no"
plot0.rebin           = 1
plot0.ymin            = 0
plot0.ymax            = 700
plot0.xmin            = 0
plot0.xmax            = 200
#plot0.lpos = "bottom-center"
plot0.name            = "Mee_allPreviousCuts_ylin"
plot0.addZUncBand     = zUncBand
plot0.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)

plot0_ylog = Plot()
## inputs for stacked histograms
plot0_ylog.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot0_ylog.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot0_ylog.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot0_ylog.keys            = keys
plot0_ylog.xtit            = "M(ee) (GeV/c^{2})"
plot0_ylog.ytit            = "Number of events"
plot0_ylog.ylog            = "yes"
plot0_ylog.rebin           = 1
plot0_ylog.ymin            = 0.001
plot0_ylog.ymax            = 1000
plot0_ylog.xmin            = 0
plot0_ylog.xmax            = 1000
#plot0_ylog.lpos = "bottom-center"
plot0_ylog.name            = "Mee_allPreviousCuts"
plot0_ylog.addZUncBand     = zUncBand
plot0_ylog.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- nEle_PtCut_IDISO_noOvrlp ---
variableName = "nEle_PtCut_IDISO_noOvrlp"

plot1 = Plot()
## inputs for stacked histograms
plot1.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot1.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot1.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot1.keys            = keys
plot1.xtit            = "Number of electrons"
plot1.ytit            = "Number of events"
plot1.ylog            = "yes"
plot1.rebin           = 1
plot1.ymin            = 0.0001
plot1.ymax            = 60000000
#plot1.lpos = "bottom-center"
plot1.name            = "nEle_allPreviousCuts"
plot1.addZUncBand     = zUncBand
plot1.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)




#--- Pt1stEle_PAS  ---
variableName = "Pt1stEle_PAS"

plot2 = Plot()
## inputs for stacked histograms
plot2.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot2.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot2.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot2.keys            = keys
plot2.xtit            = "pT 1st electron (GeV/c)"
plot2.ytit            = "Number of events"
plot2.ylog            = "yes"
plot2.rebin           = 1
plot2.xmin            = pt_xmin
plot2.xmax            = pt_xmax
plot2.ymin            = pt_ymin
plot2.ymax            = pt_ymax
#plot2.lpos = "bottom-center"
plot2.name            = "pT1stEle_allPreviousCuts"
plot2.addZUncBand     = zUncBand
plot2.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- Eta1stEle_PAS  ---
variableName = "Eta1stEle_PAS"

plot3 = Plot()
## inputs for stacked histograms
plot3.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot3.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot3.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot3.keys            = keys
plot3.xtit            = "#eta 1st electron"
plot3.ytit            = "Number of events"
plot3.rebin           = eta_rebin
plot3.ymin            = eta_ymin
plot3.ymax            = eta_ymax
plot3.lpos = "top-left"
plot3.name            = "Eta1stEle_allPreviousCuts"
plot3.addZUncBand     = zUncBand
plot3.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)



#--- Pt2ndEle_PAS  ---
variableName = "Pt2ndEle_PAS"

plot4 = Plot()
## inputs for stacked histograms
plot4.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot4.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot4.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot4.keys            = keys
plot4.xtit            = "pT 2nd electron (GeV/c)"
plot4.ytit            = "Number of events"
plot4.ylog            = "yes"
plot4.rebin           = 1
plot4.xmin            = pt_xmin
plot4.xmax            = pt_xmax
plot4.ymin            = pt_ymin
plot4.ymax            = pt_ymax
#plot4.lpos = "bottom-center"
plot4.name            = "pT2ndEle_allPreviousCuts"
plot4.addZUncBand     = zUncBand
plot4.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- Eta2ndEle_PAS  ---
variableName = "Eta2ndEle_PAS"

plot5 = Plot()
## inputs for stacked histograms
plot5.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot5.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot5.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot5.keys            = keys
plot5.xtit            = "#eta 2nd electron"
plot5.ytit            = "Number of events"
plot5.rebin           = eta_rebin
plot5.ymin            = eta_ymin
plot5.ymax            = eta_ymax
plot5.lpos = "top-left"
plot5.name            = "Eta2ndEle_allPreviousCuts"
plot5.addZUncBand     = zUncBand
plot5.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- nJet_TwoEleOnly_EtaCut ---
variableName = "nJet_TwoEleOnly_EtaCut"

plot6 = Plot()
## inputs for stacked histograms
plot6.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot6.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot6.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot6.keys            = keys
plot6.xtit            = "Number of jets"
plot6.ytit            = "Number of events"
plot6.ylog            = "yes"
plot6.rebin           = 1
plot6.xmin            = 0
plot6.xmax            = 12
plot6.ymin            = 0.01
plot6.ymax            = 2000
#plot6.lpos = "bottom-center"
plot6.name            = "nJet_allPreviousCuts"
plot6.addZUncBand     = zUncBand
plot6.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)



#--- Pt1stJet_PAS ---
variableName = "Pt1stJet_PAS"

plot7 = Plot()
## inputs for stacked histograms
plot7.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot7.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot7.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot7.keys            = keys
plot7.xtit            = "pT 1st jet (GeV/c)"
plot7.ytit            = "Number of events"
plot7.ylog            = "yes"
plot7.rebin           = 1
plot7.xmin            = pt_xmin
plot7.xmax            = pt_xmax
plot7.ymin            = pt_ymin
plot7.ymax            = pt_ymax
#plot7.lpos = "bottom-center"
plot7.name            = "Pt1stJet_allPreviousCuts"
plot7.addZUncBand     = zUncBand
plot7.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- Eta1stJet_PAS ---
variableName = "Eta1stJet_PAS"

plot8 = Plot()
## inputs for stacked histograms
plot8.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot8.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot8.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot8.keys            = keys
plot8.xtit            = "#eta 1st jet"
plot8.ytit            = "Number of events"
plot8.rebin           = eta_rebin
plot8.ymin            = eta_ymin
plot8.ymax            = eta_ymax
plot8.lpos = "top-left"
plot8.name            = "Eta1stJet_allPreviousCuts"
plot8.addZUncBand     = zUncBand
plot8.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- Pt2ndJet_PAS ---
variableName = "Pt2ndJet_PAS"

plot9 = Plot()
## inputs for stacked histograms
plot9.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot9.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot9.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot9.keys            = keys
plot9.xtit            = "pT 2nd jet (GeV/c)"
plot9.ytit            = "Number of events"
plot9.ylog            = "yes"
plot9.xmin            = pt_xmin
plot9.xmax            = pt_xmax
plot9.ymin            = pt_ymin
plot9.ymax            = pt_ymax
#plot9.lpos = "bottom-center"
plot9.name            = "Pt2ndJet_allPreviousCuts"
plot9.addZUncBand     = zUncBand
plot9.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- Eta2ndJet_PAS ---
variableName = "Eta2ndJet_PAS"

plot10 = Plot()
## inputs for stacked histograms
plot10.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot10.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot10.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot10.keys            = keys
plot10.xtit            = "#eta 2nd jet"
plot10.ytit            = "Number of events"
plot10.rebin           = eta_rebin
plot10.ymin            = eta_ymin
plot10.ymax            = eta_ymax
plot10.lpos = "top-left"
plot10.name            = "Eta2ndJet_allPreviousCuts"
plot10.addZUncBand     = zUncBand
plot10.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)

#--- sT ---
variableName = "sT_PAS"

plot11 = Plot()
## inputs for stacked histograms
plot11.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot11.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot11.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot11.keys            = keys
plot11.xtit            = "St (GeV/c)"
plot11.ytit            = "Number of events"
#plot11.xlog            = "yes"
plot11.ylog            = "yes"
plot11.rebin           = 2
plot11.xmin            = 50
plot11.xmax            = 1000
plot11.ymin            = 0.001
plot11.ymax            = 100
#plot11.lpos = "bottom-center"
plot11.name            = "sT_allPreviousCuts"
plot11.addZUncBand     = zUncBand
plot11.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- sTele ---
variableName = "sTele_PAS"

plot11_ele = Plot()
## inputs for stacked histograms
plot11_ele.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot11_ele.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot11_ele.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot11_ele.keys            = keys
plot11_ele.xtit            = "St electrons (GeV/c)"
plot11_ele.ytit            = "Number of events"
#plot11_ele.xlog            = "yes"
plot11_ele.ylog            = "yes"
plot11_ele.rebin           = 2
plot11_ele.xmin            = 50
plot11_ele.xmax            = 1000
plot11_ele.ymin            = 0.001
plot11_ele.ymax            = 100
#plot11_ele.lpos = "bottom-center"
plot11_ele.name            = "sTele_allPreviousCuts"
plot11_ele.addZUncBand     = zUncBand
plot11_ele.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- sTjet ---
variableName = "sTjet_PAS"

plot11_jet = Plot()
## inputs for stacked histograms
plot11_jet.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot11_jet.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot11_jet.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot11_jet.keys            = keys
plot11_jet.xtit            = "St jets (GeV/c)"
plot11_jet.ytit            = "Number of events"
#plot11_jet.xlog            = "yes"
plot11_jet.ylog            = "yes"
plot11_jet.rebin           = 2
plot11_jet.xmin            = 50
plot11_jet.xmax            = 1000
plot11_jet.ymin            = 0.001
plot11_jet.ymax            = 100
#plot11_jet.lpos = "bottom-center"
plot11_jet.name            = "sTjet_allPreviousCuts"
plot11_jet.addZUncBand     = zUncBand
plot11_jet.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)




##--- Mej preselection
variableNames = ["Mej_1stPair_PAS", "Mej_2ndPair_PAS"]

plot12 = Plot()
#plot12.histosStack     = [h_Mej_presel_TTbar, h_Mej_presel_ZJetAlpgen, h_Mej_presel_OTHERBKG]
plot12.histosStack     = generateAndAddHistoList( histoBaseName, samplesForStackHistos, variableNames, File_preselection)
plot12.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot12.histos          = generateAndAddHistoList( histoBaseName, samplesForHistos, variableNames, File_preselection)
plot12.keys            = keys
plot12.xtit            = "Mej (GeV/c^{2})"
plot12.ytit            = "Number of events x 2"
plot12.ylog            = "yes"
plot12.rebin           = 2
plot12.xmin            = 0
plot12.xmax            = 1000
plot12.ymin            = 0.001
plot12.ymax            = 500
#plot12.lpos = "bottom-center"
plot12.name            = "Mej_allPreviousCuts"
plot12.addZUncBand     = zUncBand
#plot12.histodata       = h_Mej_presel_DATA
plot12.histodata       = generateAndAddHisto( histoBaseName, sampleForDataHisto, variableNames, File_preselection)


#--- Mee_PAS (after preselection) ---
variableName = "Mee_PAS"

plot13 = Plot()
plot13.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot13.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot13.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot13.keys            = keys
plot13.xtit            = "M(ee) (GeV/c^{2})"
plot13.ytit            = "Number of events"
# plot13.ylog            = "yes"
# plot13.rebin           = 1
# plot13.ymin            = 0.00000001
# plot13.ymax            = 20
plot13.ylog            = "no"
plot13.rebin           = 1
plot13.ymin            = 0
plot13.ymax            = 40
plot13.xmin            = 0
plot13.xmax            = 200
#plot13.lpos = "bottom-center"
plot13.name            = "Mee_FullPreSel_allPreviousCuts_ylin"
plot13.addZUncBand     = zUncBand
plot13.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)

plot13_ylog = Plot()
## inputs for stacked histograms
plot13_ylog.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot13_ylog.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot13_ylog.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot13_ylog.keys            = keys
plot13_ylog.xtit            = "M(ee) (GeV/c^{2})"
plot13_ylog.ytit            = "Number of events"
plot13_ylog.ylog            = "yes"
plot13_ylog.rebin           = 1 # don't change it (since a rebinning is already applied above on the same histo)
plot13_ylog.ymin            = 0.001
plot13_ylog.ymax            = 100
plot13_ylog.xmin            = 0
plot13_ylog.xmax            = 1000
#plot13_ylog.lpos = "bottom-center"
plot13_ylog.name            = "Mee_FullPreSel_allPreviousCuts"
plot13_ylog.addZUncBand     = zUncBand
plot13_ylog.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)



#--- Mjj_PAS (after preselection) ---
variableName = "Mjj_PAS"

plot14 = Plot()
## inputs for stacked histograms
## it created h_Mjj_FullPreSel_TTbar, h_Mjj_FullPreSel_TTbar+h_Mjj_FullPreSel_ZJetAlpgen , h_Mjj_FullPreSel_TTbar+h_Mjj_FullPreSel_ZJetAlpgen+h_Mjj_FullPreSel_QCD_Madgraph etc..
## and plot them one on top of each other to effectly create a stacked histogram
plot14.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot14.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot14.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot14.keys            = keys
plot14.xtit            = "M(jj) (GeV/c^{2})"
plot14.ytit            = "Number of events"
# plot14.ylog            = "yes"
# plot14.rebin           = 1
# plot14.ymin            = 0.00000001
# plot14.ymax            = 20
plot14.ylog            = "no"
plot14.rebin           = 1
plot14.ymin            = 0
plot14.ymax            = 10
plot14.xmin            = 0
plot14.xmax            = 1000
#plot14.lpos = "bottom-center"
plot14.name            = "Mjj_FullPreSel_allPreviousCuts_ylin"
plot14.addZUncBand     = zUncBand
plot14.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)

plot14_ylog = Plot()
## inputs for stacked histograms
## it created h_Mjj_FullPreSel_TTbar, h_Mjj_FullPreSel_TTbar+h_Mjj_FullPreSel_ZJetAlpgen , h_Mjj_FullPreSel_TTbar+h_Mjj_FullPreSel_ZJetAlpgen+h_Mjj_FullPreSel_QCD_Madgraph etc..
## and plot them one on top of each other to effectly create a stacked histogram
plot14_ylog.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot14_ylog.keysStack       = keysStack

## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot14_ylog.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot14_ylog.keys            = keys
plot14_ylog.xtit            = "M(jj) (GeV/c^{2})"
plot14_ylog.ytit            = "Number of events"
plot14_ylog.ylog            = "yes"
plot14_ylog.rebin           = 1 # don't change it (since a rebinning is already applied above on the same histo)
plot14_ylog.ymin            = 0.001
plot14_ylog.ymax            = 100
plot14_ylog.xmin            = 0
plot14_ylog.xmax            = 1000
#plot14_ylog.lpos = "bottom-center"
plot14_ylog.name            = "Mjj_FullPreSel_allPreviousCuts"
plot14_ylog.addZUncBand     = zUncBand
plot14_ylog.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


##--- Pt Eles AllPreviousCuts ---
variableNames = ["Pt1stEle_PAS","Pt2ndEle_PAS"]

plot2and4 = Plot()
## inputs for stacked histograms
plot2and4.histosStack     = generateAndAddHistoList( histoBaseName, samplesForStackHistos, variableNames, File_preselection)
plot2and4.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot2and4.histos          = generateAndAddHistoList( histoBaseName, samplesForHistos, variableNames, File_preselection)
plot2and4.keys            = keys
plot2and4.xtit            = "pT electrons (GeV/c)"
plot2and4.ytit            = "Number of events x 2"
plot2and4.ylog            = "yes"
plot2and4.rebin           = 1
plot2and4.xmin            = pt_xmin
plot2and4.xmax            = pt_xmax
plot2and4.ymin            = pt_ymin
plot2and4.ymax            = pt_ymax
#plot2and4.lpos = "bottom-center"
plot2and4.name            = "pTEles_allPreviousCuts"
plot2and4.addZUncBand     = zUncBand
plot2and4.histodata       = generateAndAddHisto( histoBaseName, sampleForDataHisto, variableNames, File_preselection)


##--- Eta Eles AllPreviousCuts ---
variableNames = ["Eta1stEle_PAS","Eta2ndEle_PAS"]

plot3and5 = Plot()
## inputs for stacked histograms
plot3and5.histosStack     = generateAndAddHistoList( histoBaseName, samplesForStackHistos, variableNames, File_preselection)
plot3and5.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot3and5.histos          = generateAndAddHistoList( histoBaseName, samplesForHistos, variableNames, File_preselection)
plot3and5.keys            = keys
plot3and5.xtit            = "#eta electrons"
plot3and5.ytit            = "Number of events x 2"
plot3and5.rebin           = eta_rebin/2
plot3and5.ymin            = eta_ymin
plot3and5.ymax            = eta_ymax
plot3and5.lpos            = "top-left"
#plot3and5.lpos = "bottom-center"
plot3and5.name            = "etaEles_allPreviousCuts"
plot3and5.addZUncBand     = zUncBand
plot3and5.histodata       = generateAndAddHisto( histoBaseName, sampleForDataHisto, variableNames, File_preselection)


##--- Pt Jets AllPreviousCuts ---
variableNames = ["Pt1stJet_PAS","Pt2ndJet_PAS"]

plot7and9 = Plot()
## inputs for stacked histograms
plot7and9.histosStack     = generateAndAddHistoList( histoBaseName, samplesForStackHistos, variableNames, File_preselection)
plot7and9.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot7and9.histos          = generateAndAddHistoList( histoBaseName, samplesForHistos, variableNames, File_preselection)
plot7and9.keys            = keys
plot7and9.xtit            = "pT jets (GeV/c)"
plot7and9.ytit            = "Number of events x 2"
plot7and9.ylog            = "yes"
plot7and9.rebin           = 1
plot7and9.xmin            = pt_xmin
plot7and9.xmax            = pt_xmax
plot7and9.ymin            = pt_ymin
plot7and9.ymax            = pt_ymax
#plot7and9.lpos = "bottom-center"
plot7and9.name            = "pTJets_allPreviousCuts"
plot7and9.addZUncBand     = zUncBand
plot7and9.histodata       = generateAndAddHisto( histoBaseName, sampleForDataHisto, variableNames, File_preselection)

##--- Eta Eles AllPreviousCuts ---
variableNames = ["Eta1stJet_PAS","Eta2ndJet_PAS"]

plot8and10 = Plot()
## inputs for stacked histograms
plot8and10.histosStack     = generateAndAddHistoList( histoBaseName, samplesForStackHistos, variableNames, File_preselection)
plot8and10.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot8and10.histos          = generateAndAddHistoList( histoBaseName, samplesForHistos, variableNames, File_preselection)
plot8and10.keys            = keys
plot8and10.xtit            = "#eta jets"
plot8and10.ytit            = "Number of events x 2"
plot8and10.rebin           = eta_rebin/2
plot8and10.ymin            = eta_ymin
plot8and10.ymax            = eta_ymax
plot8and10.lpos            = "top-left"
#plot8and10.lpos = "bottom-center"
plot8and10.name            = "etaJets_allPreviousCuts"
plot8and10.addZUncBand     = zUncBand
plot8and10.histodata       = generateAndAddHisto( histoBaseName, sampleForDataHisto, variableNames, File_preselection)


#--- pfMET ---
variableName = "pfMET_PAS"

plot15 = Plot()
## inputs for stacked histograms
plot15.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot15.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot15.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot15.keys            = keys
plot15.xtit            = "pfMET (GeV/c)"
plot15.ytit            = "Number of events"
#plot15.xlog            = "yes"
plot15.ylog            = "yes"
plot15.rebin           = 1
plot15.xmin            = 0
plot15.xmax            = 300
plot15.ymin            = 0.001
plot15.ymax            = 400
#plot15.lpos = "bottom-center"
plot15.name            = "pfMET_allPreviousCuts"
plot15.addZUncBand     = zUncBand
plot15.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- caloMET ---
variableName = "caloMET_PAS"

plot16 = Plot()
## inputs for stacked histograms
plot16.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot16.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot16.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot16.keys            = keys
plot16.xtit            = "caloMET (GeV/c)"
plot16.ytit            = "Number of events"
#plot16.xlog            = "yes"
plot16.ylog            = "yes"
plot16.rebin           = 1
plot16.xmin            = 0
plot16.xmax            = 300
plot16.ymin            = 0.001
plot16.ymax            = 500
#plot16.lpos = "bottom-center"
plot16.name            = "caloMET_allPreviousCuts"
plot16.addZUncBand     = zUncBand
plot16.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- Pt1stEle_IDISO_NoOvrlp  ---
variableName = "Pt1stEle_IDISO_NoOvrlp"

plot2_nojet = Plot()
## inputs for stacked histograms
plot2_nojet.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot2_nojet.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot2_nojet.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot2_nojet.keys            = keys
plot2_nojet.xtit            = "pT 1st electron (GeV/c)"
plot2_nojet.ytit            = "Number of events"
plot2_nojet.ylog            = "yes"
plot2_nojet.rebin           = 1
plot2_nojet.xmin            = pt_xmin
plot2_nojet.xmax            = pt_xmax
plot2_nojet.ymin            = pt_ymin
plot2_nojet.ymax            = pt_ymax*10
#plot2_nojet.lpos = "bottom-center"
plot2_nojet.name            = "pT1stEle_nojet_allPreviousCuts"
plot2_nojet.addZUncBand     = zUncBand
plot2_nojet.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- Eta1stEle_IDISO_NoOvrlp  ---
variableName = "Eta1stEle_IDISO_NoOvrlp"

plot3_nojet = Plot()
## inputs for stacked histograms
plot3_nojet.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot3_nojet.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot3_nojet.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot3_nojet.keys            = keys
plot3_nojet.xtit            = "#eta 1st electron"
plot3_nojet.ytit            = "Number of events"
plot3_nojet.rebin           = eta_rebin
plot3_nojet.ymin            = eta_ymin
plot3_nojet.ymax            = eta_ymax*10
plot3_nojet.lpos = "top-left"
plot3_nojet.name            = "Eta1stEle_nojet_allPreviousCuts"
plot3_nojet.addZUncBand     = zUncBand
plot3_nojet.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- Pt2ndEle_IDISO_NoOvrlp  ---
variableName = "Pt2ndEle_IDISO_NoOvrlp"

plot4_nojet = Plot()
## inputs for stacked histograms
plot4_nojet.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot4_nojet.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot4_nojet.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot4_nojet.keys            = keys
plot4_nojet.xtit            = "pT 2nd electron (GeV/c)"
plot4_nojet.ytit            = "Number of events"
plot4_nojet.ylog            = "yes"
plot4_nojet.rebin           = 1
plot4_nojet.xmin            = pt_xmin
plot4_nojet.xmax            = pt_xmax
plot4_nojet.ymin            = pt_ymin
plot4_nojet.ymax            = pt_ymax*10
#plot4_nojet.lpos = "bottom-center"
plot4_nojet.name            = "pT2ndEle_nojet_allPreviousCuts"
plot4_nojet.addZUncBand     = zUncBand
plot4_nojet.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)


#--- Eta2ndEle_IDISO_NoOvrlp  ---
variableName = "Eta2ndEle_IDISO_NoOvrlp"

plot5_nojet = Plot()
## inputs for stacked histograms
plot5_nojet.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_preselection)
plot5_nojet.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot5_nojet.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_preselection)
plot5_nojet.keys            = keys
plot5_nojet.xtit            = "#eta 2nd electron"
plot5_nojet.ytit            = "Number of events"
plot5_nojet.rebin           = eta_rebin
plot5_nojet.ymin            = eta_ymin
plot5_nojet.ymax            = eta_ymax*10
plot5_nojet.lpos = "top-left"
plot5_nojet.name            = "Eta2ndEle_nojet_allPreviousCuts"
plot5_nojet.addZUncBand     = zUncBand
plot5_nojet.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_preselection)



# ############################ Plots below to be done after full selection ######################

##--- sT AllOtherCuts ---
variableName = "sT"

plot20 = Plot()
plot20.histosStack     = generateHistoList( histoBaseName, samplesForStackHistos, variableName, File_selection)
plot20.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot20.histos          = generateHistoList( histoBaseName, samplesForHistos, variableName, File_selection)
plot20.keys            = keys
plot20.xtit            = "St (GeV/c)"
plot20.ytit            = "Number of events"
plot20.ylog            = "yes"
plot20.rebin           = 10
plot20.xmin            = 0
plot20.xmax            = 1000
plot20.ymin            = 0.001
plot20.ymax            = 100
#plot20.lpos = "bottom-center"
plot20.name            = "sT_allOtherCuts"
plot20.addZUncBand     = zUncBand
plot20.histodata       = generateHisto( histoBaseName, sampleForDataHisto, variableName, File_selection)


##--- Mej AllOtherCuts ---
variableNames = ["Mej_1stPair","Mej_2ndPair"]

plot21 = Plot()
plot21.histosStack     = generateAndAddHistoList( histoBaseName, samplesForStackHistos, variableNames, File_selection)
plot21.keysStack       = keysStack
## this is the list of histograms that should be simply overlaid on top of the stacked histogram
plot21.histos          = generateAndAddHistoList( histoBaseName, samplesForHistos, variableNames, File_selection)
plot21.keys            = keys
plot21.xtit            = "Mej (GeV/c^{2})"
plot21.ytit            = "Number of events x 2"
plot21.ylog            = "yes"
plot21.rebin           = 1
plot21.xmin            = 0
plot21.xmax            = 600
plot21.ymin            = 0.001
plot21.ymax            = 20
#plot21.lpos = "bottom-center"
plot21.name            = "Mej_allOtherCuts"
plot21.addZUncBand     = zUncBand
plot21.histodata       = generateAndAddHisto( histoBaseName, sampleForDataHisto, variableNames, File_selection)


#-----------------------------------------------------------------------------------


# List of plots to be plotted
plots = [plot0, plot0_ylog, plot1, plot2_nojet, plot2, plot3_nojet, plot3, plot4_nojet, plot4, plot5_nojet, plot5,
         plot2and4, plot3and5, plot6, plot7, plot8,
         plot9, plot10, plot7and9, plot8and10, plot11, plot11_ele, plot11_jet, plot12, plot13, plot13_ylog, plot14, plot14_ylog,
         plot15, plot16,  # produced using preselection root file
         plot20, plot21] # produced using full selection root file



############# USER CODE - END ################################################
##############################################################################


#--- Generate and print the plots from the list 'plots' define above
c = TCanvas()
c.Print("allPlots.ps[")
for plot in plots:
    plot.Draw("allPlots.ps")
c.Print("allPlots.ps]")