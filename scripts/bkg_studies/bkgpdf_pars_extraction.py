import ROOT
from utilities.dataset_utils import import_pdf_library
from utilities.base_library import bin_dictionary, binning


bin_dict = bin_dictionary("pt_tracking", "eta")

par_names = ["alpha", "beta", "gamma"]

filename = "/scratch/rforti/tnp_efficiencies_results/tracking/fit_background/ws_background_fit.root"

import_pdf_library("RooCMSShape")

file = ROOT.TFile(filename)
ws = file.Get("w")

file_out = ROOT.TFile("bkgpdf_pars.root", "recreate")

h2d_alpha = ROOT.TH2D("h2d_alpha", "alpha", len(binning("pt_tracking"))-1, binning("pt_tracking"), len(binning("eta"))-1, binning("eta"))
h2d_beta = ROOT.TH2D("h2d_beta", "beta", len(binning("pt_tracking"))-1, binning("pt_tracking"), len(binning("eta"))-1, binning("eta"))
h2d_gamma = ROOT.TH2D("h2d_gamma", "gamma", len(binning("pt_tracking"))-1, binning("pt_tracking"), len(binning("eta"))-1, binning("eta"))



for b_key, [_, pt_idx, eta_idx] in bin_dict.items():

    res = ws.obj("results_fail_"+b_key)

    finalpars = res.floatParsFinal()

    h2d_alpha.SetBinContent(pt_idx, eta_idx, finalpars.find(f"alpha_fail_{b_key}").getVal())
    h2d_beta.SetBinContent(pt_idx, eta_idx, finalpars.find(f"beta_fail_{b_key}").getVal())
    h2d_gamma.SetBinContent(pt_idx, eta_idx, finalpars.find(f"gamma_fail_{b_key}").getVal())

    h2d_alpha.SetBinError(pt_idx, eta_idx, finalpars.find(f"alpha_fail_{b_key}").getError())
    h2d_beta.SetBinError(pt_idx, eta_idx, finalpars.find(f"beta_fail_{b_key}").getError())
    h2d_gamma.SetBinError(pt_idx, eta_idx, finalpars.find(f"gamma_fail_{b_key}").getError())


file_out.cd()
h2d_alpha.Write()
h2d_beta.Write()
h2d_gamma.Write()
file_out.Close()

