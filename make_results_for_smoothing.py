import ROOT
from utilities.base_library import binning, bin_dictionary, eval_efficiency, sumw2_error
from utilities.results_utils import results_manager
from utilities.dataset_utils import import_pdf_library

import_pdf_library("RooCMSShape")


ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

file_nominal = ROOT.TFile("/scratch/rforti/tnp_efficiencies_results/tracking/BBlight_legacySettings/ws_tracking_BBlight.root")
ws_nominal = file_nominal.Get("w")

file_alt = ROOT.TFile("/scratch/rforti/tnp_efficiencies_results/tracking/BBlight_legacySettings_dRsig05/ws_tracking_BBlight.root")
ws_alt = file_alt.Get("w")


res_nominal = results_manager("indep", "pt_tracking", "eta", import_ws=ws_nominal)
res_alt = results_manager("indep", "pt_tracking", "eta", import_ws=ws_alt)

h2d_eff_mc = ROOT.TH2D("EffMC2D", "EffMC2D", len(binning("eta"))-1, binning("eta"), len(binning("pt_tracking"))-1, binning("pt_tracking"))
h2d_eff_nomi = ROOT.TH2D("EffData2D", "EffData2D", len(binning("eta"))-1, binning("eta"), len(binning("pt_tracking"))-1, binning("pt_tracking"))
h2d_eff_alt = ROOT.TH2D("EffDataAltBkg2D", "EffDataAltBkg2D", len(binning("eta"))-1, binning("eta"), len(binning("pt_tracking"))-1, binning("pt_tracking"))
h2d_sf_nomi = ROOT.TH2D("SF2D_nominal", "SF2D_nominal", len(binning("eta"))-1, binning("eta"), len(binning("pt_tracking"))-1, binning("pt_tracking"))
h2d_sf_alt = ROOT.TH2D("SF2D_dataAltBkg", "SF2D_dataAltBkg", len(binning("eta"))-1, binning("eta"), len(binning("pt_tracking"))-1, binning("pt_tracking"))


for bin_key, [gl_idx, pt_idx, eta_idx] in bin_dictionary("pt_tracking", "eta").items():

    eff_nomi, d_eff_nomi = res_nominal.getEff(bin_key)
    eff_alt, d_eff_alt = res_alt.getEff(bin_key)

    eff_mc, d_eff_mc = eval_efficiency(ws_nominal.data(f"Minv_mc_pass_{bin_key}").sumEntries(), 
                                       ws_nominal.data(f"Minv_mc_fail_{bin_key}").sumEntries(), 
                                       sumw2_error(ws_nominal.data(f"Minv_mc_pass_{bin_key}")),
                                       sumw2_error(ws_nominal.data(f"Minv_mc_fail_{bin_key}")))
    
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


file_out = ROOT.TFile("/scratch/rforti/tnp_efficiencies_results/tracking/BBlight_legacySettings_dRsig05/allEfficiencies_2D.root", "RECREATE")
file_out.cd()
h2d_eff_mc.Write()
h2d_eff_nomi.Write()
h2d_eff_alt.Write()
h2d_sf_nomi.Write()
h2d_sf_alt.Write()
file_out.Close()


                                



