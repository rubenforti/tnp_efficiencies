import ROOT
from utilities.base_lib import binning, bin_dictionary, eval_efficiency, sumw2_error
from utilities.res_tools import results_manager
from utilities.dataset_utils import import_pdf_library

import_pdf_library("RooCMSShape")

ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

ch = "plus"

'''
file_nominal = ROOT.TFile("/scratch/rforti/tnp_efficiencies_results/tracking/BBlight_legacySettings/ws_tracking_BBlight.root")
ws_nominal = file_nominal.Get("w")
res_nominal = results_manager("indep", "pt_tracking", "eta", import_ws=ws_nominal)
'''

benchmark_filename = f"/scratch/rforti/tnp_efficiencies_results/egm_tnp_results/tracking{ch}/allEfficiencies.txt"
with open(benchmark_filename, "r") as file_bmark:
    row_list = file_bmark.readlines()
    res_nominal = results_manager("indep", "pt_tracking", "eta", import_txt=row_list)

file_alt = ROOT.TFile(f"/scratch/rforti/tnp_efficiencies_results/tracking{ch}/BBlight_legacySettings/ws_tracking_BBlight.root")
ws_alt = file_alt.Get("w")
res_alt = results_manager("indep", "pt_tracking", "eta", import_ws=ws_alt)

h2d_eff_mc = ROOT.TH2D("EffMC2D", "EffMC2D", len(binning("eta"))-1, binning("eta"), len(binning("pt_tracking"))-1, binning("pt_tracking"))
h2d_eff_nomi = ROOT.TH2D("EffData2D", "EffData2D", len(binning("eta"))-1, binning("eta"), len(binning("pt_tracking"))-1, binning("pt_tracking"))
h2d_eff_alt = ROOT.TH2D("EffDataAltBkg2D", "EffDataAltBkg2D", len(binning("eta"))-1, binning("eta"), len(binning("pt_tracking"))-1, binning("pt_tracking"))
h2d_sf_nomi = ROOT.TH2D("SF2D_nominal", "SF2D_nominal", len(binning("eta"))-1, binning("eta"), len(binning("pt_tracking"))-1, binning("pt_tracking"))
h2d_sf_alt = ROOT.TH2D("SF2D_dataAltBkg", "SF2D_dataAltBkg", len(binning("eta"))-1, binning("eta"), len(binning("pt_tracking"))-1, binning("pt_tracking"))


for bin_key, [gl_idx, pt_idx, eta_idx] in bin_dictionary("pt_tracking", "eta").items():

    eff_nomi, d_eff_nomi = res_nominal.getEff(bin_key)
    eff_alt, d_eff_alt = res_alt.getEff(bin_key)

    eff_mc, d_eff_mc = res_nominal.getEffMC(bin_key)

    h2d_eff_mc.SetBinContent(eta_idx, pt_idx, eff_mc)
    h2d_eff_mc.SetBinError(eta_idx, pt_idx, d_eff_mc)

    h2d_eff_nomi.SetBinContent(eta_idx, pt_idx, eff_nomi)
    h2d_eff_nomi.SetBinError(eta_idx, pt_idx, d_eff_nomi)

    h2d_eff_alt.SetBinContent(eta_idx, pt_idx, eff_alt)
    h2d_eff_alt.SetBinError(eta_idx, pt_idx, d_eff_alt)

    h2d_sf_nomi.SetBinContent(eta_idx, pt_idx, eff_nomi/eff_mc)
    h2d_sf_nomi.SetBinError(eta_idx, pt_idx, d_eff_nomi/eff_mc)

    h2d_sf_alt.SetBinContent(eta_idx, pt_idx, eff_alt/eff_mc)
    h2d_sf_alt.SetBinError(eta_idx, pt_idx, d_eff_alt/eff_mc)


file_out = ROOT.TFile(f"/scratch/rforti/tnp_efficiencies_results/tracking{ch}/BBlight_legacySettings/allEfficiencies_2D.root", "RECREATE")
file_out.cd()
h2d_eff_mc.Write()
h2d_eff_nomi.Write()
h2d_eff_alt.Write()
h2d_sf_nomi.Write()
h2d_sf_alt.Write()
file_out.Close()


                                



