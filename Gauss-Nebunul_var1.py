import ROOT
from ROOT import RooMinimizer,TPaveText,RooDataSet,RooNLLVar,TArrayD, RooDataHist, RooClassFactory, TFile, RooRealVar, TH1D, TH1F, RooAbsReal, TLegend, TCanvas, RooArgList, RooGenericPdf, RooArgSet, RooWorkspace, gROOT, RooAddPdf, TLatex, TH2D, TF1
import sys
from ROOT import RooFit
from ROOT import RooCBShape, RooPolynomial, RooExponential, RooGaussian
import array, math
from ROOT import gStyle, kBird, kAquamarine, kBeach, kSolar, kCool
from ROOT import gROOT
from ROOT import TStyle
from ROOT import TLatex
import math
import numpy as np
import scipy.integrate 
from scipy.integrate import quad
gStyle.SetPalette(kBeach);
gStyle.SetPaintTextFormat("4.1f");
from Ostap.LHCbStyle import LHCbStyle
LHCbStyle(name='lhcbStyle2', desc='Standard "LHCb-like" plots style', lineWidth=2, font=132, makeNew=True, force=True, scale=1, colz=True)
MyFile = ROOT.TFile.Open("/data/DATA/data.lhcb/SMOG2/Data/255623/tuples/data_run_255623_AlignmentV11_Moore_v54r5_hlt2_thor_SMOG2_2022_Full.root")
MyTree = MyFile.Get("LambdappiLL/DecayFunTuple")
nLeafs = MyTree.GetListOfBranches()
Lambda_MASSMC = TH1D("Lambda_MASS","Lambda_MASS",50,1090,1140)
entries = MyTree.GetEntries()
can1 = TCanvas('can1', 'can1')
signal = "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[5])^2)"
signal_2 = "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[5])^2)"
background = "([6]+x*[7])"
combined = signal +  "+" + background 
func = ROOT.TF1("func", combined, 1090, 1140)
func.SetParameters(148000,1116,1,42000,2.3,1116,-350000,320)
f1 = ROOT.TF1("f1", signal_2, 1090, 1140)
func.SetParLimits(5,-1000000,1000000)
for i in range(entries):
    MyTree.GetEntry(i)
    Lambda_MASSMC.Fill(MyTree.Lambda_MASS)
Lambda_MASSMC.Sumw2()
fitResult = Lambda_MASSMC.Fit("func","S")
X = Lambda_MASSMC.GetXaxis()
X.SetTitle("Invariant Mass p#pi(#frac{MeV}{c^{2}})")
Y = Lambda_MASSMC.GetYaxis()
Y.SetTitle("No. of candidates (#frac{MeV}{c^{2}})")
fit = Lambda_MASSMC.GetFunction("func")
chi2 = fit.GetChisquare()
can1 = TCanvas('can1', 'can1')
Lambda_MASSMC.SetTitle("Lambda_MASS")
Lambda_MASSMC.Draw()
cov_matrix = fitResult.GetCovarianceMatrix()
params = np.array([fitResult.Parameter(i) for i in range(func.GetNpar())])
errors = np.array([fitResult.ParError(i) for i in range(func.GetNpar())])
ndof = fitResult.Ndf()
latex = TLatex()
latex.SetNDC()
latex.SetTextSize(0.04)
latex.SetTextAlign(13)
latex.SetTextColor(ROOT.kRed)
latex.DrawLatex(0.18, 0.9, "#mu_{1} = %.2f #pm %.2f" %(func.GetParameter(1), func.GetParError(1)))
latex.DrawLatex(0.18, 0.8, "#mu_{2} = %.2f #pm %.2f" %(func.GetParameter(4), func.GetParError(4)))
latex.DrawLatex(0.18, 0.7, "#sigma_{1} = %.2f #pm %.2f" %(func.GetParameter(2), func.GetParError(2)))
latex.DrawLatex(0.18, 0.6, "#sigma_{2} = %.2f #pm %.2f" %(func.GetParameter(5), func.GetParError(5)))
latex.DrawLatex(0.18, 0.5, "#frac{#chi^{2}}{ndof} = %.2f" %(chi2/ndof))
func.Draw("same")
f1.SetParameters(func.GetParameter(0),func.GetParameter(1),func.GetParameter(2),func.GetParameter(3),func.GetParameter(4), func.GetParameter(5))
cov_matrix1 = ROOT.TMatrixDSym(6)
for i in range(6):
    for j in range(6):
        cov_matrix1[i][j] = cov_matrix[i][j]
params1 = np.array([fitResult.Parameter(i) for i in range(6)])
errors1 = np.array([fitResult.ParError(i) for i in range(6)])
r = f1.Integral(1090,1140)
er = f1.IntegralError(1090,1140,params1,cov_matrix1.GetMatrixArray())
latex.DrawLatex(0.16,0.3,"#signals = %.2f#pm%.5f" %(r,er))
s1 = '[0]*exp(-0.5*((x-[1])/[2])^2)'
fc1 = ROOT.TF1("fc1", s1, 1090, 1140)
fc1.SetParameters(func.GetParameter(0),func.GetParameter(1),func.GetParameter(2))
cov_matrix11 = ROOT.TMatrixDSym(3)
for i in range(3):
    for j in range(3):
        cov_matrix11[i][j] = cov_matrix[i][j]
params11 = np.array([fitResult.Parameter(i) for i in range(3)])
errors11 = np.array([fitResult.ParError(i) for i in range(3)])
a1 = float(fc1.Integral(1090,1140))
ea1 = float(fc1.IntegralError(1090,1140, params11,cov_matrix11.GetMatrixArray()))
s2 = '[0]*exp(-0.5*((x-[1])/[2])^2)'
fc2 = ROOT.TF1("fc1", s2, 1090, 1140)
fc2.SetParameters(func.GetParameter(3),func.GetParameter(4),func.GetParameter(5))
cov_matrix22 = ROOT.TMatrixDSym(3)
for i in range(3,6):
    for j in range(3,6):
        cov_matrix22[i-3][j-3] = cov_matrix[i][j]
params22 = np.array([fitResult.Parameter(i) for i in range(3)])
errors22 = np.array([fitResult.ParError(i) for i in range(3)])
a2 = float(fc2.Integral(1090,1140))
ea2 = float(fc2.IntegralError(1090,1140, params22,cov_matrix22.GetMatrixArray()))
sgm = ((float(func.GetParameter(2)))*(float(a1))+ float(func.GetParameter(4))*(float(a2)))/(float(a1)+float(a2))
sg1 = float(func.GetParameter(2))
esg1 = float(func.GetParError(2))
sg2 = float(func.GetParameter(5))
esg2 = float(func.GetParError(5))
esgm = math.sqrt(((a1**2)*(esg1**2)+(a2**2)*(esg2**2))/((a1+a2)**2)+((sg1-sg2)**2)*((a2**2)*(ea1**2)+(a1**2)*(ea2**2))/((a1+a2)**4))
latex.DrawLatex(0.18,0.4, "<#sigma> = %.2f#pm %.4f" %(sgm,esgm))
print(chi2,ndof)
can1.Update()