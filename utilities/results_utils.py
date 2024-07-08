"""
"""
import sys
import ROOT


def create_histograms(names, bins_var=None, bins_pt=None, bins_eta=None):
    """
    Generic function to create histograms, that can be 1D or 2D depending on
    the input binning. Returns a dictionary with the histograms.
    """
    histos = {}

    for name in names:

        if bins_var is not None:
            histos[name] = ROOT.TH1D(name, name, len(bins_var)-1, bins_var)
        if (bins_pt is not None) and (bins_eta is not None):
            histos[f"{name}_2d"] = ROOT.TH2D(f"{name}_2d", name, len(bins_pt)-1, bins_pt, len(bins_eta)-1, bins_eta)

    return histos

###############################################################################


def init_pass_fail_histos(histo_names, bins_var=None, bins_pt=None, bins_eta=None, flag=""):

    histos = {}

    flag_list = ["pass", "fail"] if flag=="" else [flag]

    pf_histos = {}
    
    for fl in flag_list:
        histo_names_flagged = [f"{histo_name}_{fl}" for histo_name in histo_names]
        pf_histos.update( create_histograms(histo_names_flagged, bins_var, bins_pt, bins_eta) )

    return pf_histos
    
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
