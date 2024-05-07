"""
"""
import ROOT
import math
from utilities.base_lib import bin_dictionary, binning
from utilities.dataset_utils import import_pdf_library

gen_dir = "/scratch/rforti/tnp_efficiencies_results/tracking"

filename = gen_dir + "/legacy_fit_onlyFailSA_allMC/ws_tracking.root"

import_pdf_library("RooCMSShape")

f = ROOT.TFile(filename, "READ")

ws = f.Get("w")

bin_dict = bin_dictionary("pt_tracking", "eta")

file_out = ROOT.TFile("utilities_true/pars_atLim/parsAtLim_tracking.root", "RECREATE")

pars_dict = {
    "total": ROOT.TH2D("total", "At least one par. at limit", len(binning("pt_tracking"))-1, binning("pt_tracking"), len(binning("eta"))-1, binning("eta")),
    "nsig": ROOT.TH2D("nsig", "nsig", len(binning("pt_tracking"))-1, binning("pt_tracking"), len(binning("eta"))-1, binning("eta")),
    "nbkg": ROOT.TH2D("nbkg", "nbkg", len(binning("pt_tracking"))-1, binning("pt_tracking"), len(binning("eta"))-1, binning("eta")),
    "mu": ROOT.TH2D("mu", "#mu conv.", len(binning("pt_tracking"))-1, binning("pt_tracking"), len(binning("eta"))-1, binning("eta")),
    "sigma": ROOT.TH2D("sigma", "#sigma conv.", len(binning("pt_tracking"))-1, binning("pt_tracking"), len(binning("eta"))-1, binning("eta")),
    "alpha": ROOT.TH2D("alpha", "#alpha", len(binning("pt_tracking"))-1, binning("pt_tracking"), len(binning("eta"))-1, binning("eta")),
    "beta": ROOT.TH2D("beta", "#beta", len(binning("pt_tracking"))-1, binning("pt_tracking"), len(binning("eta"))-1, binning("eta")),
    "gamma": ROOT.TH2D("gamma", "#gamma", len(binning("pt_tracking"))-1, binning("pt_tracking"), len(binning("eta"))-1, binning("eta")),
}

absTol = 1e-4
relTol = 1e-5

for b_key, [_, idx_pt, idx_eta] in bin_dict.items():

    res = ws.obj(f"results_fail_{b_key}")
    pars = res.floatParsFinal()

    n_pars_at_lim = 0

    for par_name, par_hist in pars_dict.items():
        
        if par_name=="total": continue

        par = pars.find(par_name+"_fail_"+b_key)

        if math.isclose(par.getVal(), par.getMin(), abs_tol=absTol, rel_tol=relTol):
            par_hist.SetBinContent(idx_pt, idx_eta, -1)
            n_pars_at_lim += 1
        elif math.isclose(par.getVal(), par.getMax(), abs_tol=absTol, rel_tol=relTol):
            par_hist.SetBinContent(idx_pt, idx_eta, 1)
            n_pars_at_lim += 1
        else:
            par_hist.SetBinContent(idx_pt, idx_eta, 0)

    if n_pars_at_lim > 0:
        pars_dict["total"].SetBinContent(idx_pt, idx_eta, 1)
    else:
        pars_dict["total"].SetBinContent(idx_pt, idx_eta, 0)

ROOT.gStyle.SetOptStat(1000000)
ROOT.gStyle.SetStatY(0.95)
ROOT.gStyle.SetStatX(0.9)

for par_hist in pars_dict.values():
    file_out.cd()
    par_hist.Write()

    c = ROOT.TCanvas(par_hist.GetName(), par_hist.GetName(), 800, 600)
    c.cd()
    par_hist.Draw("colz")
    par_hist.SetTitleOffset(0.6, "Y")
    z_axis = par_hist.GetZaxis()
    z_axis.SetRangeUser(-1, 1)
    z_axis.SetNdivisions(3)
    par_hist.GetXaxis().SetTitle("p_{T} [GeV]")
    par_hist.GetXaxis().SetTitleOffset(1.2)
    par_hist.GetYaxis().SetTitle("#eta")
    par_hist.GetYaxis().SetTitleOffset(0.8)

    c.SaveAs(f"utilities_true/pars_atLim/{par_hist.GetName()}_atLim.png")




file_out.Close()

