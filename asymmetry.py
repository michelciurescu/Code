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
import scipy.integrate 
from scipy.integrate import quad
gStyle.SetPalette(kBeach);
gStyle.SetPaintTextFormat("4.1f");
from Ostap.LHCbStyle import LHCbStyle
LHCbStyle(name='lhcbStyle2', desc='Standard "LHCb-like" plots style', lineWidth=2, font=132, makeNew=True, force=True, scale=1, colz=True)
MyFile = ROOT.TFile.Open("/data/DATA/data.lhcb/SMOG2/Data/255623/tuples/data_run_255623_AlignmentV11_Moore_v54r5_hlt2_thor_SMOG2_2022_Full.root")
# davinci_data_run_255623_hlt2_thor_AlignmentV11_SMOG2_Moorev54r5_noMinBiasnoHidCharm_2022_Final.root
MyTree = MyFile.Get("LambdappiLL/DecayFunTuple")
entries = MyTree.GetEntries()
nLeafs = MyTree.GetListOfBranches()git
signal = "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)"
signal_2 = "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1])/[4])^2)"
background = "sqrt([5]+x*[6])"
combined = signal +  "+" + background 
func = ROOT.TF1("func", combined, 1090, 1140)
func.SetParameters(140000,1116,1.3,42000,2, -3000,450) 
f1 = ROOT.TF1("f1", signal_2, 1090, 1140)
func.SetParLimits(5,-1000000,1000000)
hist = TH1D("Lambda_MASS", "Lambda_MASS", 50 ,1090, 1140)
hist.SetTitle("No. of #bar{#Lambda_{0}} Candidates as a Function of Their Mass")
X = hist.GetXaxis()
X.SetTitle("Invariant Mass p#pi(#frac{MeV}{c^{2}})")
Y = hist.GetYaxis()
Y.SetTitle("No. of #bar{#Lambda_{0}} candidates(#frac{MeV}{c^{2}})")
for i in range(entries):
    MyTree.GetEntry(i)
    if MyTree.Lambda_ID < 0 :
        hist.Fill(MyTree.Lambda_MASS)
fitResult = hist.Fit("func","S")
cov_matrix = fitResult.GetCovarianceMatrix()
params = np.array([fitResult.Parameter(i) for i in range(7)])
errors = np.array([fitResult.ParError(i) for i in range(7)])
ndof = fitResult.Ndf()
ndof = fitResult.Ndf()
fit = hist.GetFunction("func")
chi2 = fit.GetChisquare()
func.SetLineColor(ROOT.kGreen)
func.SetLineColor(ROOT.kGreen)
hist.Sumw2()
hist.Fit("func")
hist.Draw()
func.Draw("Same")
can1 = TCanvas('can1', 'can1')
hist.Draw()
latex = TLatex()
latex.SetNDC()
latex.SetTextSize(0.04)
latex.SetTextAlign(13)
latex.SetTextColor(ROOT.kRed)
latex.DrawLatex(0.2, 0.9, "#mu = %.2f #pm %.2f" %(func.GetParameter(1), func.GetParError(1)))
latex.DrawLatex(0.2, 0.8, "#sigma_{1} = %.2f #pm %.2f" %(func.GetParameter(2), func.GetParError(2)))
latex.DrawLatex(0.2, 0.7, "#sigma_{2} = %.2f #pm %.2f" %(func.GetParameter(4), func.GetParError(4)))
latex.DrawLatex(0.2, 0.6, "#frac{#chi^{2}}{ndof} = %.2f" %(chi2/ndof))
f1.SetParameters(func.GetParameter(0),func.GetParameter(1),func.GetParameter(2),func.GetParameter(3),func.GetParameter(4))
cov_matrix1 = ROOT.TMatrixDSym(5)
for i in range(5):
    for j in range(5):
        cov_matrix1[i][j] = cov_matrix[i][j]
params1 = np.array([fitResult.Parameter(i) for i in range(5)])
errors1 = np.array([fitResult.ParError(i) for i in range(5)])
nosgn = f1.Integral(1090, 1140)
nosgn_error = f1.IntegralError(1090, 1140,params1, cov_matrix1.GetMatrixArray())
latex.DrawLatex(0.16,0.4,"#signals = %.2f #pm %.2f" %(nosgn, nosgn_error))
s1 = '[0]*exp(-0.5*((x-[1])/[2])^2)'
fc1 = ROOT.TF1("fc1", s1, 1090, 1140)
fc1.SetParameters(func.GetParameter(0),func.GetParameter(1),func.GetParameter(2))
a1 = float(fc1.Integral(1090,1140))
s2 = '[0]*exp(-0.5*((x-[1])/[2])^2)'
fc2 = ROOT.TF1("fc1", s2, 1090, 1140)
fc2.SetParameters(func.GetParameter(3),func.GetParameter(1),func.GetParameter(4))
a2 = float(fc2.Integral(1090,1140))
sgm = ((float(func.GetParameter(2)))*(float(a1))+ float(func.GetParameter(4))*(float(a2)))/(float(a1)+float(a2))
latex.DrawLatex(0.2,0.5, "<#sigma> = %.2f" %(sgm))
print(chi2,ndof)
func.Draw("same")
can1.Update()
"Make a change"