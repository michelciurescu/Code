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
MyFile = ROOT.TFile.Open("/data/DATA/data.lhcb/SMOG2/Data/255427/tuples/data_run_255427_AlignmentV11_Moore_v54r5_hlt2_thor_SMOG2_2022_Full.root")
MyTree = MyFile.Get("LambdappiLL/DecayFunTuple")
entries = MyTree.GetEntries()
nLeafs = MyTree.GetListOfBranches()
can1 = TCanvas('can1', 'can1')
hist = TH1D("Lambda_Y", "Lambda_Y", 35,1105,1140)
signal = "[0]*exp(-0.5*((x-[1])/[2])^2)"#+[3]*exp(-0.5*((x-[1])/[4])^2)"
background = "[3]*(1+x*[4])"
combined = signal +  "+" + background 
func = ROOT.TF1("func", combined, 1105, 1140)
func.SetParameters(500, 1115.1, 1,-10, -0.001)
X = hist.GetXaxis()
X.SetTitle("Invariant Mass p#pi(#frac{MeV}{c^{2}})")
Y = hist.GetYaxis()
Y.SetTitle("#candidates #bar{#Lambda_{0}} with 3.9#leqy<5(#frac{MeV}{c^{2}})")
for i in range(entries):
    MyTree.GetEntry(i)
    MyTree.Lambda_Y = math.atanh((float(MyTree.Lambda_PZ))/(float(MyTree.Lambda_ENERGY)))
    if 3.9 <= MyTree.Lambda_Y < 5 and MyTree.Lambda_ID < 0:
        hist.Fill(MyTree.Lambda_MASS)
func.SetParLimits(5,-1000,1000)
fitResult = hist.Fit("func","S")
cov_matrix = fitResult.GetCovarianceMatrix()
params = np.array([fitResult.Parameter(i) for i in range(func.GetNpar())])
errors = np.array([fitResult.ParError(i) for i in range(func.GetNpar())])
ndof = fitResult.Ndf()
fit = hist.GetFunction("func")
chi2 = fit.GetChisquare()
func.SetLineColor(ROOT.kGreen)
hist.Sumw2()
hist.Fit("func")
hist.Draw()
func.Draw("Same")
latex = TLatex()
latex.SetNDC()
latex.SetTextSize(0.04)
latex.SetTextAlign(13)
latex.SetTextColor(ROOT.kRed)
latex.DrawLatex(0.5, 0.9, "#mu = %.2f #pm %.2f" %(func.GetParameter(1), func.GetParError(1)))
latex.DrawLatex(0.5, 0.8, "#sigma = %.2f #pm %.2f" %(func.GetParameter(2), func.GetParError(2)))
#latex.DrawLatex(0.6, 0.7, "#sigma_{2} = %.2f #pm %.2f" %(func.GetParameter(6), func.GetParError(6)))
latex.DrawLatex(0.5, 0.7, "#frac{#chi^{2}}{ndof} = %.2f" %(chi2/ndof))
print(chi2/ndof)
integral = func.Integral(1105,1140)
integral_error = func.IntegralError(1105, 1040, params, cov_matrix.GetMatrixArray())
print(integral, integral_error)
s1 =  "[0]*exp(-0.5*((x-[1])/[2])^2)"
f1 = ROOT.TF1("f1", s1, 1105, 1140)
f1.SetParameters(func.GetParameter(0),func.GetParameter(1),func.GetParameter(2))
params1 = np.array([fitResult.Parameter(i) for i in range(3)])
errors1 = np.array([fitResult.ParError(i) for i in range(3)])
cov_matrix1 = ROOT.TMatrixDSym(3)
for i in range(3):
    for j in range(3):
        cov_matrix1[i][j] = cov_matrix[i][j]
p1 = f1.Integral(1105,1140)
dp1 = f1.IntegralError(1105,1140,params1,cov_matrix1.GetMatrixArray())
print(p1,dp1)
latex.DrawLatex(0.5,0.6,"#signals = %.2f#pm%.2f" %(p1,dp1))
#print(func.GetParameter(1),func.GetParameter(2),func.GetParameter(3))
#s2 =  "[0]*exp(-0.5*((x-[1])/[2])^2)"
#f2 = ROOT.TF1("f2", s2, 1090, 1140)
#f2.SetParameters(func.GetParameter(0),func.GetParameter(1),func.GetParameter(2))
#params2 = np.array([fitResult.Parameter(i) for i in [3,1,4]])
#errors2 = np.array([fitResult.ParError(i) for i in [3,1,4]])
#p2 = f2.Integral(1090,1140)
#dp2 = f2.IntegralError(1090,1140,params2,cov_matrix.GetMatrixArray())
#print(p2)
#e1 = (func.GetParameter(2))
#e2 = (func.GetParameter(4))
#sgma = (p1*e1+p2*e2)/(p1+p2)
#latex.DrawLatex(0.6, 0.4, "<#sigma> = %.2f" %(sgma))
can1.Update()