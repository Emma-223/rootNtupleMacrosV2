#!/usr/bin/env python3

from ROOT import TGraphErrors,TCanvas,TLegend,TLine,kBlue,TGraph, gROOT
import numpy as numpy
import math
import os
import sys


def GetScaleFactorNominal():
    with open(txtFilename,'r') as thisFile:
        nextRescaleLine = False
        for line in thisFile:
            line = line.strip()
            if nominalBkgPlotName in line:
                nextRescaleLine = True
                continue
            if 'rescale factor' in line and nextRescaleLine:
                line = line.split(':')
                numbers = line[1]
                rescale = float(numbers[0:numbers.find('+/-')])
                rescaleErr = float(numbers[numbers.find('+/-')+4:])
                return rescale, rescaleErr


if len(sys.argv) < 2:
    raise RuntimeError("Did not get any rescale log file. Must specify it as an argument to the script")
txtFilename = sys.argv[1]

# bkgName = 'DYJets'
# nominalBkgPlotName = 'Mee_BkgControlRegion_'+bkgName
bkgName = 'TTbar'
nominalBkgPlotName = 'Mee_BkgControlRegion_gteTwoBtaggedJets_'+bkgName
sfNom, sfNomErr = GetScaleFactorNominal()
varList = ['sT', 'MejMin', 'NJetEq']
#txtFilename = 'ttbar/logAllVariations.log'
#bkgName = 'ttbar'
# use sT and MejMin, one at a time
#varList = ['sT']
#varList = ['MejMin']
#varList = ['LQ']

# for enujj
##bkgName = 'wjet'
#bkgName = 'ttbar'
##varList = ['LQ']
#varList = ['sT']

gROOT.SetBatch(True)

eejjMode = True
# if 'enujj' in txtFilename:
#   eejjMode=False
#   # in enujj mode, default x/y is ttbar

for var in varList:
  xPoints = []
  xPointErrs = []
  yPoints = []
  yPointErrs = []
  if not eejjMode:
    xPointsWJ = []
    xPointErrsWJ = []
    yPointsWJ = []
    yPointErrsWJ = []
  with open(txtFilename,'r') as thisFile:
    isVarPoint = False
    for line in thisFile:
      line = line.strip()
      middle = 0.0
      delta = 0.0
      rescale = 0.0
      rescaleErr = 0.0
      if 'name:' in line and bkgName in line:
        line = line.split(':')
        name = line[1].strip()
        if not var in name:
          continue
        isVarPoint = True
        print('name=', name)
        if not any(specialVar in name for specialVar in ['LQ', 'NJetEq']):
          if not 'To' in name:
            # doesn't quite work for Menu; use separate script for plotting MT variations
            split = name.split('_')
            lowerVarBound = float(split[1])
            upperVarBound = float(split[2])
            middle = float((lowerVarBound+upperVarBound)/2.0)
            delta = upperVarBound-middle
            varBin = name[name.find(var)+len(var):name.rfind('_')]
            print('lowerVarBound=',lowerVarBound,'upperVarBound=', upperVarBound)
          else:
            #print 'varBin=',varBin
            varBin = name[name.find(var)+len(var):name.rfind('_')]
            lowerVarBound = float(varBin[0:varBin.find('To')])
            upperVarBound = float(varBin[varBin.find('To')+2:varBin.rfind('_')])
            # above is for bins; below is for threshold plots
            #lowerVarBound = float(varBin)-5
            #upperVarBound = float(varBin)+5
            if upperVarBound > 2000:
              upperVarBound = 2000
            print('lowerVarBound=',lowerVarBound,'upperVarBound=',upperVarBound)
            middle = float((lowerVarBound+upperVarBound)/2.0)
            delta = upperVarBound-middle
        else:  # other vars: LQ, NJetEq
          #print 'name=',name
          #varBin = name[name.find(var)+len(var):name.rfind('_')]
          varBin = name[name.find(var)+len(var):name.find('_', name.find(var)+len(var))]
          if len(varBin)<=0:
            varBin = name[name.find(var)+len(var):]
          # print('for var={}, varBin={}'.format(var, varBin))
          middle = float(varBin)
          delta = 0
          if var == 'NJetEq':
              delta = 0.5
              middle += 0.5
        if eejjMode or not 'WJet' in name:
          xPoints.append(middle)
          xPointErrs.append(delta)
        elif not eejjMode and 'WJet' in name:
          xPointsWJ.append(middle)
          xPointErrsWJ.append(delta)
      if 'rescale factor' in line and isVarPoint:
        line = line.split(':')
        numbers = line[1]
        #print 'numbers="'+numbers+'"'
        rescale = float(numbers[0:numbers.find('+/-')])
        rescaleErr = float(numbers[numbers.find('+/-')+4:])
        if eejjMode or not 'WJet' in name:
          print('yPoints append: rescale =',rescale,'rescaleErr=',rescaleErr)
          yPoints.append(rescale)
          if rescaleErr < 0:
            rescaleErr = math.fabs(rescaleErr)
          yPointErrs.append(rescaleErr)
        elif not eejjMode and 'WJet' in name:
          print('yPointsWJ append: rescale =',rescale,'rescaleErr=',rescaleErr)
          yPointsWJ.append(rescale)
          if rescaleErr < 0:
            rescaleErr = math.fabs(rescaleErr)
          yPointErrsWJ.append(rescaleErr)
        isVarPoint = False
 

  if len(xPoints) <= 0:
    print('No data in logfile found for var:',var)
    continue
  # make graph
  #print numpy.array(xPoints)
  #print numpy.array(xPointErrs)
  #print numpy.array(yPoints)
  #print numpy.array(yPointErrs)
  if bkgName=='wjet' and not eejjMode:
    xPoints = xPointsWJ
    xPointErrs = xPointErrsWJ
    yPoints = yPointsWJ
    yPointErrs = yPointErrsWJ
  canvas = TCanvas()
  canvas.cd()
  graph = TGraphErrors(len(xPoints),numpy.array(xPoints, dtype="f"),numpy.array(yPoints, dtype="f"),numpy.array(xPointErrs, dtype="f"),numpy.array(yPointErrs, dtype="f"))
  graph.SetTitle('')
  graph.SetMarkerColor(kBlue)
  graph.SetLineColor(kBlue)
  graph.Draw('ap')
  if var=='sT':
    graph.GetXaxis().SetTitle('S_{T} [GeV]')
  elif var=='MejMin':
    graph.GetXaxis().SetTitle('M_{ej}^{min} [GeV]')
  elif 'LQ' in var:
    graph.GetXaxis().SetTitle('M_{LQ} [GeV]')
  elif var=='NJetEq':
    graph.GetXaxis().SetTitle('N_{jet}')
  graph.GetXaxis().SetNdivisions(516)
  if eejjMode:
    if bkgName=='ttbar':
      graph.GetYaxis().SetTitle('t#bar{t} scale factor')
      #sfNom = 0.81479
      #sfNomErr = 0.036965
    elif bkgName=='DYJets':
      graph.GetYaxis().SetTitle('Z+jets scale factor')
      #sfNom = 0.94 # MGHT
      #sfNomErr = 0.01
      #sfNom = 1.03 # amc@NLO PtBinned
      #sfNomErr = 0.02
      #sfNom = 1.05 # amc@NLO PtBinned
      #sfNomErr = 0.01
      #sfNom = 0.98 # amc@NLO PtBinned, Sep. 29 PtEE
      #sfNom = 0.97 # amc@NLO PtBinned, Nov. 19, amcAtNLO diboson
      #sfNomErr = 0.01
  else:
    if bkgName=='ttbar':
      graph.GetYaxis().SetTitle('t#bar{t} scale factor')
      graph.GetYaxis().SetRangeUser(0,2)
      #sfNom = 0.95
      #sfNom = 0.92 # no top pt reweight
      #sfNomErr = 0.01
    elif bkgName=='wjet':
      graph.GetYaxis().SetTitle('W+jets scale factor')
      graph.GetYaxis().SetRangeUser(0,2)
      #sfNom = 0.87
      #sfNomErr = 0.012
  uncertaintyXpoints = [xPoints[0]-xPointErrs[0],xPoints[-1]+xPointErrs[-1]]
  uncertaintyYpoints = [sfNom,sfNom]
  uncertaintyXpointsErrs = [0,0]
  uncertaintyYpointsErrs = [sfNomErr,sfNomErr]
  uncertaintyRegionGraph = TGraphErrors(len(uncertaintyXpoints),numpy.array(uncertaintyXpoints, dtype="f"),numpy.array(uncertaintyYpoints, dtype="f"),numpy.array(uncertaintyXpointsErrs, dtype="f"),numpy.array(uncertaintyYpointsErrs, dtype="f"))
  uncertaintyRegionGraph.SetFillColor(15)
  uncertaintyRegionGraph.SetFillStyle(3001)
  uncertaintyRegionGraph.SetLineColor(1)
  uncertaintyRegionGraph.Draw('3l')
  #lineSFNom = TLine(xPoints[0]-xPointErrs[0],sfNom,xPoints[0]-xPointErrs[0],sfNom)
  #lineSFNom.SetLineColor(2)
  #lineSFNom.Draw()
  # draw the polygon for shading
  #uncertaintyXpoints = [xPoints[0]-xPointErrs[0],xPoints[-1]+xPointErrs[-1],xPoints[-1]+xPointErrs[-1],xPoints[0]-xPointErrs[0]]
  #uncertaintyYpoints = [sfNom+sfNomErr,sfNom+sfNomErr,sfNom-sfNomErr,sfNom-sfNomErr]
  #uncertaintyRegionGraph = TGraph(len(uncertaintyXpoints),numpy.array(uncertaintyXpoints),numpy.array(uncertaintyYpoints))
  #uncertaintyRegionGraph.SetFillColor(15)
  #uncertaintyRegionGraph.SetFillStyle(3001)
  #uncertaintyRegionGraph.Draw('f')
  #lineUp = TLine(graph.GetXaxis().GetXmin(),sfNom*1.1,graph.GetXaxis().GetXmax(),sfNom*1.1)
  #lineUp.SetLineStyle(2)
  #lineUp.Draw()
  #lineDown = TLine(graph.GetXaxis().GetXmin(),sfNom*0.9,graph.GetXaxis().GetXmax(),sfNom*0.9)
  #lineDown.SetLineStyle(2)
  #lineDown.Draw()
  if eejjMode:
    if bkgName.lower() == 'ttbar':
      if var=='sT':
        leg = TLegend(0.19,0.68,0.51,0.89)
      elif var=='MejMin':
        leg = TLegend(0.19,0.19,0.51,0.40)
      else:
        leg = TLegend(0.19,0.19,0.51,0.40)
      leg.AddEntry(graph,'t#bar{t} scale factor variations','lp')
    elif bkgName.lower() == 'dyjets':
      if var=='sT':
        leg = TLegend(0.19,0.19,0.51,0.40)
      elif var=='MejMin':
        leg = TLegend(0.6,0.7,0.9,0.89)
      elif var=='NJetEq':
        leg = TLegend(0.19,0.7,0.51,0.89)
      else:
        leg = TLegend(0.19,0.19,0.51,0.40)
      leg.AddEntry(graph,'DYJets scale factor variations','lp')
    else:
        raise RuntimeError("Could not detect background type; expecting TTBar or DYJets")
  else:
    if bkgName.lower() == 'ttbar':
      if var=='sT':
        leg = TLegend(0.19,0.68,0.51,0.89)
      elif var=='MejMin':
        leg = TLegend(0.19,0.19,0.51,0.40)
      else:
        leg = TLegend(0.19,0.19,0.51,0.40)
      leg.AddEntry(graph,'t#bar{t} scale factor variations','lp')
    elif bkgName.lower() == 'wjet':
      if var=='sT':
        #leg = TLegend(0.19,0.79,0.51,0.89)
        leg = TLegend(0.19,0.19,0.51,0.40)
      elif var=='MejMin':
        #leg = TLegend(0.19,0.79,0.51,0.89)
        leg = TLegend(0.19,0.19,0.51,0.40)
      else:
        leg = TLegend(0.19,0.68,0.51,0.89)
      leg.AddEntry(graph,'W+jets scale factor variations','lp')
    else:
        raise RuntimeError("Could not detect background type; expecting TTBar or WJet")
  leg.AddEntry(uncertaintyRegionGraph,'Nominal scale factor = '+str(round(sfNom,3))+' #pm '+str(round(sfNomErr,3)),'fl')
  #leg.AddEntry(lineUp,'Nominal scale factor #pm 10%','l')
  leg.SetBorderSize(0)
  leg.Draw()
  canvas.Modified()

  # print
  # if var=='sT':
  #   baseName = bkgName+'_scaleFactorVariation_sTBins' 
  # elif var=='MejMin':
  #   baseName = bkgName+'_scaleFactorVariation_mejMinBins'
  # elif 'LQ' in var:
  #   baseName = bkgName+'_scaleFactorVariation_LQBins'
  # elif 'MTenu' in var:
  #   baseName = bkgName+'_scaleFactorVariation_MTBins'
  baseName = bkgName+'_scaleFactorVariation_{}Bins'.format(var)
  canvas.Print(baseName+'.png')
  canvas.Print(baseName+'.pdf')

  # ## wait for input to keep the GUI (which lives on a ROOT event dispatcher) alive
  # if __name__ == '__main__':
  #   rep = ''
  #   while not rep in [ 'c', 'C' ]:
  #    rep = input( 'enter "c" to continue: ' )
  #    if 1 < len(rep):
  #      rep = rep[0]
print('Using nominal scale factor = {} +/- {}'.format(round(sfNom, 3), round(sfNomErr, 3)))
