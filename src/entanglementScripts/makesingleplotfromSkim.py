from __future__ import print_function
import ROOT
from ROOT import gStyle, gPad,gDirectory, gROOT
from ROOT import TCanvas, TPad
from ROOT import TFile, TH1F, TLegend, TLine
from ROOT import kBlack, kAzure, kOrange
from ROOT import kFullCircle
ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
import numpy
from numpy import log10
from optparse import OptionParser
import math as _math
import CMSStyle
import Style
import os


### Options stuff --------------------------------------------------------------------------------------------------------------------
parser = OptionParser()
parser.add_option("--Log","--isLog", dest="isLog",   default=False, action="store_true", help="Plot the plots in log ?" )
(options, args) = parser.parse_args()
Log=options.isLog


# Input Skim File and TTree
inFile = TFile.Open("/data/dust/user/ricardo/uhh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/src/entanglementScripts/spincorr_skim_TTToSemiLeptonic_UL18.root")
tree = inFile.Get("entSkim")
# gen distributions must use the GENERATOR weight, never eventweight (which carries reco SFs).
# Falls back to unweighted "1" on skims that predate the genweight branch.
GENW = "genweight" if tree.GetBranch("genweight") else "1"
print("Using GENW = %s" % GENW)

# ==============================
# histograms included in overlay
# ==============================
# hist_reco = "reco_cosTheta1k"
hist_gen  = "gen_cosTheta1k"

hist_reco = "gen_cosTheta1k" # comparing the effect of using eventweight vs genweight for gen_ histograms

# Output directory for histograms
outDirectory = "/data/dust/user/ricardo/uhh2-106X_v2/CMSSW_10_6_28/src/UHH2/ZprimeSemiLeptonic/src/entanglementScripts/skimPlots"
if not os.path.isdir(outDirectory):
    os.makedirs(outDirectory)

# plot file name
# plot_name = hist_reco.replace("reco_", "") + "_genVreco"
plot_name = "gen_cosTheta1k_eventwtVgenwt"
# # Root file name
# outFile = TFile("%s/%s.root" % (outDirectory, plot_name), "RECREATE")

# Define binning for groups of variables, or None to auto-range
def _binning(name):
    """(nbins, lo, hi) for a variable by name pattern, or None to auto-range."""
    base = name.replace("reco_", "").replace("gen_", "")
    if base == "chi2":            return None              # reco-only, range unknown -> auto
    if base == "M_tt":            return (100, 300.,  2500.)
    if base == "pt_hadTop":       return (100,   0.,  1200.)
    if base == "beta":            return ( 50,   0.,     1.)
    if base == "dyreco":          return ( 80,  -4.,     4.)
    if base.startswith("Sigma_phi") or base.startswith("Delta_phi"):
        return (64, -_math.pi, _math.pi)
    if base.endswith("_plus") or base.endswith("_minus"):  # C cross-combos live in [-2,2]
        return (80, -2., 2.)
    return (50, -1., 1.)          # cosTheta*, cHel(_P3n)(+slices), C{ij}, cosThetaStar
 
 
def _valid(name):
    """selection dropping the -10 / 99 sentinels and requiring the right flag."""
    flag = "pass_reco" if name.startswith("reco_") else "gen_semilep"
    return "(%s && %s > -9.9 && %s < 98.)" % (flag, name, name)
 
 
# Histogram Information:
gROOT.SetBatch(True)
thestyle = Style.Style()
thestyle.SetStyle()
print("AddDirectory:", ROOT.TH1.AddDirectoryStatus())   # expect False -> that's the culprit

# Canvas info:
canvas = TCanvas('c1', 'c1', 800, 600)
canvas.SetFillColor(0)
canvas.SetBorderMode(0)
canvas.SetFrameFillStyle(0)
canvas.SetFrameBorderMode(0)
canvas.SetTickx(0)

# Pad info TPad(name, title, xlow, ylow, xup, yup):
pad1 = TPad("pad1","pad1",0,0.3,1,1)
pad1.SetBottomMargin(0.02)
pad1.Draw()
# second pad for ratio plot
pad2 = TPad("pad2","pad2",0,0,1,0.3)
pad2.SetTopMargin(0.02)
pad2.Draw()
# common pad settings
pad_Lmargin = 0.12
pad_Rmargin = 0.03

# Legend
legend = TLegend(0.70, 0.78, 0.89, 0.90)
legend.SetBorderSize(0)
legend.SetFillStyle(0)


# Extract binnning info
binning = _binning(hist_reco)
nb, lo, hi = _binning(hist_reco) if binning is not None else (50, -1., 1.)

# =================================================================
# Auto-create histogram by calling Draw() to fill from TTree branch
# =================================================================
# e.g.: tree.Draw("branch_name>>hist_name", "selection_string", "goff"=graphics off)
# tree.Draw("%s>>h_%s(%d,%g,%g)" % (hist_reco, hist_reco, nb, lo, hi), "eventweight*%s" % _valid(hist_reco), "goff")  # reco -> eventweight
# tree.Draw("%s>>h_%s(%d,%g,%g)" % (hist_gen, hist_gen, nb, lo, hi),   "%s*%s" % (GENW, _valid(hist_gen)),   "goff")  # gen  -> genweight

# # both curves are the SAME branch (gen_cosTheta1k) with DIFFERENT weights -> DISTINCT names
# tree.Draw("%s>>h_eventwt(%d,%g,%g)" % (hist_gen, nb, lo, hi), "eventweight*%s" % _valid(hist_gen), "goff")  # gen x eventweight
# tree.Draw("%s>>h_genwt(%d,%g,%g)"   % (hist_gen, nb, lo, hi), "%s*%s" % (GENW, _valid(hist_gen)),  "goff")  # gen x genweight
# print("  validity of %s = %s" % (hist_gen, _valid(hist_gen)))

# h_reco = gDirectory.Get("h_eventwt")   # eventweighted gen (the incorrect weighting, for comparison)
# h_gen  = gDirectory.Get("h_genwt")     # genweighted gen   (correct)
# print("  h_reco type %s, h_gen type %s" % (type(h_reco), type(h_gen)))

# h_reco.SetDirectory(0)
# h_gen.SetDirectory(0)

tree.Draw("%s>>h_ew(%d,%g,%g)" % (hist_gen, nb, lo, hi), "eventweight*%s" % _valid(hist_gen), "goff")  # gen x eventweight
h_reco = tree.GetHistogram().Clone("h_eventwt"); h_reco.SetDirectory(0)
tree.Draw("%s>>h_gw(%d,%g,%g)" % (hist_gen, nb, lo, hi), "%s*%s" % (GENW, _valid(hist_gen)),  "goff")  # gen x genweight
h_gen  = tree.GetHistogram().Clone("h_genwt");   h_gen.SetDirectory(0)
print("  entries: eventwt=%d  genwt=%d  (type %s)" % (h_reco.GetEntries(), h_gen.GetEntries(), type(h_reco)))
print("  Events for: %s = %d,  %s = %d using Integral()"   % (hist_reco, h_reco.Integral(),   hist_gen, h_gen.Integral()))

# Normalize histograms to unit area
if h_reco.Integral() > 0: h_reco.Scale(1.0 / h_reco.Integral())
if h_gen.Integral() > 0: h_gen.Scale(1.0 / h_gen.Integral())

# Start drawing on Main pad
pad1.cd()

# Line colors: reco(blue), gen(orange)
h_reco.SetLineColor(kAzure + 1)
h_gen.SetLineColor(kOrange + 7)
# Line weights
h_reco.SetLineWidth(2)
h_gen.SetLineWidth(2)

# Common aesthetics: use either histogram
ymax = max(h_reco.GetMaximum(), h_gen.GetMaximum())
if Log:
    pad1.SetLogy()
    h_reco.SetMaximum(10**(1.25*log10(ymax)))
    h_reco.SetMinimum(1)
else:
    h_reco.SetMaximum(ymax * 1.25 if ymax > 0 else 1.0)
    h_reco.SetMinimum(0.0)
h_reco.SetTitle(hist_reco.replace("reco_", "reco vs gen: "))
h_reco.GetYaxis().SetTitle("a.u. (unit area)")
h_reco.GetXaxis().SetLabelSize(0.0)  # hide x-axis labels on top pad
h_reco.Draw("hist")
h_gen.Draw("hist same")

# Legend
# legend.AddEntry(h_reco, "reco", "l")
# legend.AddEntry(h_gen, "gen", "l")
legend.AddEntry(h_reco, "gen #times eventweight", "l")
legend.AddEntry(h_gen,  "gen #times genweight",   "l")
legend.Draw()

### CMS aesthetics
CMSStyle.setCMSEra("UL18")
text = "Work in progress"
dx = -0.11
CMSStyle.extraText = text[0] + "".join("#kern[%g]{%s}"% (dx, c) for c in text[1:])
# CMSStyle.extraText = "Work in progress: #ell +jets"
CMSStyle.extraOverCmsTextSize = 1.0
CMSStyle.setCMSLumiStyle(pad1, 0)
pad1.SetFillColor(0)
pad1.SetBorderMode(0)
pad1.SetFrameFillStyle(0)
pad1.SetTickx(0)
pad1.SetTicky(0)
pad1.SetTopMargin(0.10) # room for the plot title
pad1.SetLeftMargin(pad_Lmargin)
pad1.SetRightMargin(pad_Rmargin) 
pad1.SetBottomMargin(0) # up against the ratio plot
pad1.Update()
canvas.Update()



# Draw ratio plot on second pad
canvas.cd()
pad2.cd()
h_ratio = h_reco.Clone("h_ratio")
h_ratio.Divide(h_gen)

# Adjust zoom for ratio
ratio_max = 1.01
ratio_min = 0.99

# Common aesthetics for ratio plot
h_ratio.SetLineColor(kBlack)
h_ratio.SetMarkerStyle(20)
h_ratio.SetMarkerSize(0.5)
h_ratio.SetMaximum(ratio_max)
h_ratio.SetMinimum(ratio_min)
h_ratio.SetTitle("")
# y-axis
h_ratio.GetYaxis().SetTitle("reco/gen")
h_ratio.GetYaxis().SetNdivisions(505)
h_ratio.GetYaxis().SetTitleSize(0.11)
h_ratio.GetYaxis().SetTitleOffset(0.42)
h_ratio.GetYaxis().SetLabelSize(0.09)
# x-axis
h_ratio.GetXaxis().SetTitle(hist_reco.replace("reco_", ""))
h_ratio.GetXaxis().SetTitleSize(0.12)
h_ratio.GetXaxis().SetLabelSize(0.10)
# draw
h_ratio.Draw("ep")
# draw horizontal line at y=1
line = TLine(h_ratio.GetXaxis().GetXmin(), 1.0, h_ratio.GetXaxis().GetXmax(), 1.0)
line.SetLineStyle(2)
line.Draw()

# pad2 aesthetics
pad2.SetFillColor(0)
pad2.SetBorderMode(0)
pad2.SetFrameFillStyle(0)
pad2.SetTickx(0)
pad2.SetTicky(0)
pad2.SetTopMargin(0) # up against the main plot
pad2.SetLeftMargin(pad_Lmargin)
pad2.SetRightMargin(pad_Rmargin)
pad2.SetBottomMargin(0.3)
pad2.Update()



# Save plot: PDF gets a single .pdf; the canvas + hists also go into the .root
canvas.cd()
canvas.Update()
pdf_path = "%s/%s%s.pdf" % (outDirectory, plot_name, "_log" if Log else "")
canvas.SaveAs(pdf_path)
print("wrote %s" % pdf_path)
 
# # Write histograms and canvas to output root file
# outFile.cd()
# for _h in (h_reco, h_gen, h_ratio):
#     _h.Write()
# canvas.Write("c_" + plot_name)

# Cleanup
# outFile.Close()
inFile.Close()