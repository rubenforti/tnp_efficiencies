"""
"""
import ROOT
from copy import copy
from utilities.CMS_lumi import CMS_lumi
from utilities.plot_utils import style_settings


ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

analysis = "iso_indep_2gev"
folder = "results/iso_indep_2gev"

# -----------------------------------------------------------------------------

file_hres = ROOT.TFile(f"{folder}/res_{analysis}.root", "READ")
file_hres_cmp = ROOT.TFile(f"{folder}/res_cmp_{analysis}.root", "READ")
file_hres_cmp_bmark = ROOT.TFile(f"{folder}/res_cmpBmark_{analysis}.root", "READ")

histos = {}
#[histos.update({key.GetName() : file_hres.Get(key.GetName())}) for key in file_hres.GetListOfKeys()]
[histos.update({key.GetName() : file_hres_cmp.Get(key.GetName())}) for key in file_hres_cmp.GetListOfKeys()]
#[histos.update({key.GetName() : file_hres_cmp_bmark.Get(key.GetName())}) for key in file_hres_cmp_bmark.GetListOfKeys()]

for hist_key in histos.keys():

    style_settings() 

    hist = histos[hist_key]   

    if "2d" in hist_key:
        x_dim = 1600
        y_dim = 1200
        c = ROOT.TCanvas(hist_key, hist_key, x_dim, y_dim)
        c.cd()
        pad_plot = ROOT.TPad("pad_plot", "pad_plot", 0, 0, 1, 0.92)
        pad_title = ROOT.TPad("pad_title", "pad_title", 0, 0.92, 1, 1)
        pad_plot.SetMargin(0.08, 0.12, 0.12, 0.05), pad_plot.Draw()
        pad_title.SetMargin(0.15, 0.15, 0.15, 0.15), pad_title.Draw()
        plot_name = hist.GetTitle()
    else:
        ROOT.gStyle.SetOptStat("emruo")
        ROOT.gStyle.SetStatBorderSize(2)
        if hist_key == "h_efficiency": 
            ROOT.gStyle.SetStatX(0.6)
        else:
            ROOT.gStyle.SetStatX(0.9)
        ROOT.gStyle.SetStatY(0.95)
        x_dim = 1200
        y_dim = 1200
        c = ROOT.TCanvas(hist_key, hist_key, x_dim, y_dim)
        c.cd()
        pad_plot = ROOT.TPad("pad_plot", "pad_plot", 0, 0, 1, 0.95)
        pad_title = ROOT.TPad("pad_title", "pad_title", 0, 0.95, 1, 1)
        pad_plot.SetMargin(0.12, 0.08, 0.12, 0.04), pad_plot.Draw()
        pad_title.SetMargin(0., 0., 0., 0.), pad_title.Draw()
        x_title = copy(hist.GetTitle())
        plot_name = ""

    pad_plot.cd()
    hist.Draw("COLZ")
    hist.GetZaxis().SetTitle("")
    hist.SetTitleSize(0)
    hist.SetTitle("")
    
    if "2d" in hist_key:
        hist.GetXaxis().SetTitle("p_{T}^{#mu} [GeV]")
        hist.GetXaxis().SetTitleOffset(1.2)
        hist.GetXaxis().SetTitleSize(0.04)
        hist.GetYaxis().SetTitle("#eta^{#mu}")
        hist.GetYaxis().SetTitleOffset(0.8)
        hist.GetYaxis().SetTitleSize(0.04)
    else:
        hist.GetXaxis().SetTitle(x_title)
        hist.GetXaxis().SetTitleOffset(1.2)
        hist.GetXaxis().SetTitleSize(0.04)
        hist.GetYaxis().SetTitle("Events")
        hist.GetYaxis().SetTitleOffset(1.25)
        hist.GetYaxis().SetTitleSize(0.04)
        pad_plot.SetLogy()
    
    
    
    CMS_lumi(pad_plot, 5, 0)
    pad_plot.Update()
    
    pad_title.cd()
    titlebox = ROOT.TPaveText(0, 0.05, 1., 0.95, "NDC NB")
    titlebox.SetFillColor(0)
    titlebox.SetTextSize(0.6)
    titlebox.AddText(plot_name)
    titlebox.Draw()
    pad_title.Update()

    c.Update()

    c.SaveAs(f"{folder}/res_plots/{hist_key}.pdf")


