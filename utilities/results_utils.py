"""
"""
import sys
import ROOT


def init_pass_fail_histos(histo_name, histo_title, bins_var, bins_pt, bins_eta, flag=""):

    histos = {}

    flag_list = ["pass", "fail"] if flag=="" else [flag]

    for fl in flag_list:
        h1d = ROOT.TH1D(f"{histo_name}_{flag}", f"{histo_title} - {fl}ing probes", 
                        len(bins_var)-1, bins_var)
        h2d = ROOT.TH2D(f"{histo_name}_{flag}_2d", f"{histo_title} - {fl}ing probes", 
                        len(bins_pt)-1, bins_pt, len(bins_eta)-1, bins_eta)
        histos.update({f"{histo_name}_{fl}" : h1d})
        histos.update({f"{histo_name}_{fl}_2d" : h2d})

    return histos
    
###############################################################################


def init_results_histos(histo_name, histo_title, bins_var, bins_pt, bins_eta):

    histos = {}

    h1d = ROOT.TH1D(f"h_{histo_name}", histo_title, len(bins_var)-1, bins_var)
    h2d = ROOT.TH2D(f"h_{histo_name}_2d", f"{histo_title} ",
                    len(bins_pt)-1, bins_pt, len(bins_eta)-1, bins_eta)
    histos.update({f"{histo_name}" : h1d})
    histos.update({f"{histo_name}_2d" : h2d})

    return histos

###############################################################################


def fill_res_histograms(res_bmark, res_new, hist_dict, bin_dict):
    """
    """
    for bin_key, [gl_idx, bin_pt, bin_eta] in bin_dict.items():

        if type(bin_eta) is list or type(bin_eta) is list:
            sys.exit("ERROR: binning not correcty defined, use get_mergedbins_bounds=False in dictionary generation")

        eff_bmark, deff_bmark = res_bmark.getEff(bin_key)
        eff_test, deff_test = res_new.getEff(bin_key)

        if "delta" in hist_dict.keys():
            hist_dict["delta"].Fill(eff_test-eff_bmark)
            hist_dict["delta_2d"].SetBinContent(bin_pt, bin_eta, eff_test-eff_bmark)
        if "delta_error" in hist_dict.keys():
            hist_dict["delta_error"].Fill(deff_test-deff_bmark)
            hist_dict["delta_error_2d"].SetBinContent(bin_pt, bin_eta, deff_test-deff_bmark)
        if "pull" in hist_dict.keys():
            hist_dict["pull"].Fill((eff_test-eff_bmark)/(deff_bmark**2 + deff_test**2)**0.5)
            hist_dict["pull_2d"].SetBinContent(bin_pt, bin_eta, (eff_test-eff_bmark)/((deff_bmark**2 + deff_test**2)**0.5))
        if "pull_ref" in hist_dict.keys():
            hist_dict["pull_ref"].Fill((eff_test-eff_bmark)/deff_bmark)
            hist_dict["pull_ref_2d"].SetBinContent(bin_pt, bin_eta, (eff_test-eff_bmark)/deff_bmark)
        if "rm1" in hist_dict.keys():
            hist_dict["rm1"].Fill((eff_test/eff_bmark)-1)
            hist_dict["rm1_2d"].SetBinContent(bin_pt, bin_eta, (eff_test/eff_bmark)-1)
        if "ratio_error" in hist_dict.keys():
            hist_dict["ratio_error"].Fill((deff_test/deff_bmark)-1)
            hist_dict["ratio_error_2d"].SetBinContent(bin_pt, bin_eta, (deff_test/deff_bmark)-1)
    
        
###############################################################################
###############################################################################
