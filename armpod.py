import ROOT
from ROOT import RooMinimizer,TPaveText,RooDataSet,RooNLLVar,TArrayD, RooDataHist, RooClassFactory, TFile, RooRealVar, TH1D, TH1F, RooAbsReal, TLegend, TCanvas, RooArgList, RooGenericPdf, RooArgSet, RooWorkspace, gROOT, RooAddPdf, TLatex, TH2D, TF1,TGraph,TGraphErrors
import sys
from ROOT import RooFit
from ROOT import RooCBShape, RooPolynomial, RooExponential, RooGaussian
import array, math
from ROOT import gStyle, kBird, kAquamarine, kBeach, kSolar, kCool
from ROOT import gROOT
from ROOT import TStyle
import numpy as np
import scipy.integrate 
from scipy.integrate import quad
gStyle.SetPalette(kBeach);
gStyle.SetPaintTextFormat("4.1f");
from Ostap.LHCbStyle import LHCbStyle
LHCbStyle(name='lhcbStyle2', desc='Standard "LHCb-like" plots style', lineWidth=2, font=132, makeNew=True, force=True, scale=1, colz=True)
MyFile = ROOT.TFile.Open("/data/DATA/data.lhcb/SMOG2/Data/255623/tuples/data_run_255623_AlignmentV11_Moore_v54r5_hlt2_thor_SMOG2_2022_Full.root")
MyTree = MyFile.Get("LambdappiLL/DecayFunTuple")
entries = MyTree.GetEntries()
nLeafs = MyTree.GetListOfBranches()
can1 = TCanvas('can1', 'can1')
hist1 = TH2D("h2", "h2 title", 800, -1, 1, 80, 4000, 4000)
hist2 = TH2D("h2","h2title",800, -1, 1, 80, 4000, 4000 )
for i in range(entries):
    MyTree.GetEntry(i)
    if MyTree.Lambda_ID > 0:
         MyTree.alpha = (MyTree.pbar_PZ - MyTree.pim_PZ)/(MyTree.pbar_PZ+MyTree.pim_PZ)
         hist1.Fill(MyTree.alpha, MyTree.Lambda_PT)
    if MyTree.Lambda_ID < 0:
         MyTree.alpha = -(MyTree.pbar_PZ - MyTree.pim_PZ)/(MyTree.pbar_PZ+MyTree.pim_PZ)
         hist1.Fill(MyTree.alpha, MyTree.Lambda_PT)
X = hist1.GetXaxis()
X.SetTitle("#alpha")
Y = hist1.GetYaxis()
Y.SetTitle("p_{t}(#frac{Mev}{c})")
hist1.Draw("colz")
hist2.Draw("same")
can1.Draw()
can1.Update()
