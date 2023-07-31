import ROOT
from ROOT import RooMinimizer,TPaveText,RooDataSet,RooNLLVar,TArrayD, RooDataHist, RooClassFactory, TFile, RooRealVar, TH1D, TH1F, RooAbsReal, TLegend, TCanvas, RooArgList, RooGenericPdf, RooArgSet, RooWorkspace, gROOT, RooAddPdf, TLatex, TH2D, TF1
import sys
from ROOT import RooFit
from ROOT import RooCBShape, RooPolynomial, RooExponential, RooGaussian
import array, math

from ROOT import gStyle, kBird, kAquamarine, kBeach, kSolar, kCool
from ROOT import gROOT
from ROOT import TStyle
import numpy as np
gStyle.SetPalette(kBeach);
gStyle.SetPaintTextFormat("4.1f");


from Ostap.LHCbStyle import LHCbStyle
LHCbStyle(name='lhcbStyle2', desc='Standard "LHCb-like" plots style', lineWidth=2, font=132, makeNew=True, force=True, scale=1, colz=True)

#Change the directory 
#l = 0
MyFile = ROOT.TFile.Open("/data/DATA/data.lhcb/SMOG2/Data/255427/tuples/data_run_255427_AlignmentV11_Moore_v54r5_hlt2_thor_SMOG2_2022_Full.root", "READ")
MyTree = MyFile.Get("LambdappiLL/DecayFunTuple")
nLeafs = MyTree.GetListOfBranches()
entries = MyTree.GetEntries()
can = []
hist1 = ROOT.TH1F("hist", "Mass for NPV = 1",150,1100,1140)
hist2 = ROOT.TH1F("hist", "Mass for NPV = 2",150,1100,1140)
hist3 = ROOT.TH1F("hist", "Mass for NPV > 2",150,1100,1140)
X = hist1.GetXaxis()
X.SetTitle("Invariant Mass p#pi(#frac{MeV}{c^{2}})")
Y = hist1.GetYaxis()
Y.SetTitle("No. of candidates (#frac{MeV}{c^{2}})")
for i in range(entries):
    MyTree.GetEntry(i)
    if MyTree.nPVs == 1:
        hist1. Fill(MyTree.Lambda_MASS)
    if MyTree.nPVs == 2:
        hist2.Fill(MyTree.Lambda_MASS)
    if MyTree.nPVs > 2:
        hist3.Fill(MyTree.Lambda_MASS)
hist1.scale()
hist2.scale()
hist3.scale()
can1 = TCanvas('can1', 'can1')
hist1.SetLineColor(ROOT.kBlue)
hist2.SetLineColor(ROOT.kRed)
hist3.SetLineColor(ROOT.kGreen)
hist1.SetMarkerColor(ROOT.kBlue)
hist2.SetMarkerColor(ROOT.kRed)
hist3.SetMarkerColor(ROOT.kGreen)
hist1.SetTitle("Histo")
hist2.SetTitle("Histo")
hist3.SetTitle("Histo")
hist1.Draw()
hist2.Draw("same")
hist3.Draw("same")
legend = ROOT.TLegend(0.1, 0.7, 0.48, 0.9)
legend.AddEntry(hist1,"Mass for NPV = 1")
legend.AddEntry(hist2,"Mass for NPV = 2")
legend.AddEntry(hist3,"Mass for NPV > 2")
legend.Draw()
can1.Draw()