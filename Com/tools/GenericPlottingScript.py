###############################################################################################################
#####
##### This script takes two input root files and overlays the plots of all of the histograms stored inside them. 
#####
##### To run: python GenericPlottingScript.py
#####
################################################################################################################

import ROOT
import sys
import os
import re
import string
from DataFormats.FWLite import Events, Handle

from ROOT import TCanvas, TFile, TProfile, TNtuple, TH1F, TH2F,TF1, TPad, TPaveLabel, TPaveText, TLegend
from ROOT import gROOT, gBenchmark, gRandom, gSystem, Double

gROOT.Reset()

import Style
#import tdrstyle
thestyle = Style.Style()
thestyle.SetStyle()

#######################################
### List the files you want to run over
#######################################

filelist = []
#filelist.append("/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/May6thPDSkim_GOODCOLL_v1_5202010/CaloLooseQaulityCuts/All.root")
#filelist.append("/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/MinBias_pythia8_Spring10-START3X_V1_5202010/CaloLooseQaulityCuts/All.root")

filelist.append("/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/ForPAS_812010/Full_DataSample_15GeV/All.root")
filelist.append("/uscms_data/d2/mike1886/EwkMETCommissioning/FWLiteAnalyzerOuput/ForPAS_712010/W_Jets/All.root")



################################################################################################################################
### For each file, go into the subdirectory /histos and get a complete list of histograms. Store all these histograms in Hist[]
################################################################################################################################

def GetHistos(tfile):
    Hist = []
    dir=tfile.GetDirectory('histos')
    alist=dir.GetListOfKeys()
    print alist
    for j in alist:
        obj = j.ReadObj()
        if obj.IsA().InheritsFrom(ROOT.TH1.Class()):
            print '  --> found TH1: name = '+j.GetName() + ' title1 = '+j.GetTitle()
            Hist.append(obj)

    print '--------------------------------------------------------------------'                                
    #print "length of MC hist = ",len(Hist)
   
    return Hist




################################################################################################################################
### Draw all Historgrams overlayed for Data and MC. These are passed in by Hist[][]
################################################################################################################################

canvas = {}

def plot(Data,MC,prefix,k):

    if MC.GetEntries() > 0:
        canvas[k] = ROOT.TCanvas( str(prefix) + str(Data.GetName()), str(prefix) + str(Data.GetName()), 1000,700)
        canvas[k].cd()
        canvas[k].SetLogy(1)
        MC.Scale(Data.GetEntries()/MC.GetEntries())
        MC.SetFillColor(2)
        MC.SetLineColor(2)
        MC.GetXaxis().SetTitle(str(MC.GetTitle()))
        MC.GetYaxis().SetTitle('Number of Events')
        MC.Draw()
        Data.SetMarkerStyle(23)
        Data.Draw('SAMES:E')
        #Legend(Data ,MC)
    
        leg=ROOT.TLegend(0.656878,0.769717,0.834672,0.918527)
        leg.AddEntry(Data,'DATA',"lp")
        leg.AddEntry(MC,'MC',"lp")
        leg.SetShadowColor(0)
        leg.SetFillColor(0)
        leg.SetLineColor(0)
        leg.Draw()   

        pl = ROOT.TPaveLabel()
        pl.SetBorderSize(0);
        pl.SetLabel('I CAN ADD TEXT');
        pl.SetTextSize(.5);
        pl.SetFillColor(0);
        pl.SetX1NDC (0.62751);
        pl.SetY1NDC (0.697917);
        pl.SetX2NDC (0.893574);
        pl.SetY2NDC (0.770833);
        #pl.Draw();

        canvas[k].SaveAs('PLOTS/PlotsOf'+str(prefix)+'_'+str(Data.GetName())+'.jpg')
    

### THIS IS NOT WORKING WHEN CALLED FROM PLOT(), IT JUSTS DOES NOTHING!!!
def Legend(Data,MC):
   
    leg = ROOT.TLegend(0.656878,0.769717,0.834672,0.918527)
    leg.AddEntry(MC,"Simulation","f")
    leg.AddEntry(DATA,"DATA","lp");
    leg.SetShadowColor(0)
    leg.SetFillColor(0)
    leg.SetLineColor(0)
    leg.Draw()    


#######################################################################################################
### Loop over each file. Get the list of Histograms (GetHistos()) and then plot each Histogram (plot())
#######################################################################################################
    
files = []
prefix = ['DATA','MC']
j = 0
Histos = {}

for i in filelist:

    files.append(ROOT.TFile(i))
    print 'Now running over ',files[j].GetName()
    Histos[j] = GetHistos(files[j])
   
    j = j + 1


### All Histograms are now stored in Histos[][]. Histos[0][i] is all the Data histograms and Histos[1][i] is all the MC historams

for k in range(len(Histos[0])):
    plot(Histos[0][k],Histos[1][k],prefix[1],k)



  
