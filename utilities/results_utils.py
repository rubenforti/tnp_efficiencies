"""
"""
import sys
import ROOT
from array import array
from utilities.base_lib import binnings


def create_histograms(names, titles, bins_var=None, bins_pt=None, bins_eta=None):
    """
    Generic function to create histograms, that can be 1D or 2D depending on
    the input binning. Returns a dictionary with the histograms.
    """
    histos = {}

    for name, title in zip(names , titles):
        if bins_var is not None:
            if type(bins_var) is array:
                histos[name] = ROOT.TH1D(name, title, len(bins_var)-1, bins_var)
            elif type(bins_var) is dict:
                histos[name] = ROOT.TH1D(name, title, len(bins_var[name])-1, bins_var[name])
            else:
                sys.exit("ERROR: wrong binning definition")

        if (bins_pt is not None) and (bins_eta is not None):
            histos[f"{name}_2d"] = ROOT.TH2D(f"{name}_2d", title, len(bins_pt)-1, bins_pt, len(bins_eta)-1, bins_eta)

    return histos

###############################################################################


def init_pass_fail_histos(histo_names, bins_var=None, bins_pt=None, bins_eta=None, flag=""):

    flag_list = ["pass", "fail"] if flag=="" else [flag]

    pf_histos = {}
    
    for fl in flag_list:
        histo_names_flagged = [f"{histo_name}_{fl}" for histo_name in histo_names]
        pf_histos.update( create_histograms(histo_names_flagged, bins_var, bins_pt, bins_eta) )

    return pf_histos

###############################################################################


def init_histos(histo_names, histo_titles, binning_var=None, binning_pt=None, binning_eta=None, flag="all"):

    histos = {}

    bins_var = binnings[binning_var] if type(binning_var) is str else binning_var

    bins_pt, bins_eta = binnings[binning_pt], binnings[binning_eta]

    flag_list = ["pass", "fail"] if flag=="all" else [flag]
    
    for fl in flag_list:
        fl_postfix = "" if fl=="" else f"_{fl}"

        hists_fl = create_histograms([f"{histo_name}{fl_postfix}" for histo_name in histo_names], 
                                     [f"{histo_title}{fl_postfix}" for histo_title in histo_titles],
                                      bins_var, bins_pt, bins_eta )
        histos.update(hists_fl)

    return histos
    
###############################################################################


def init_results_histos(histo_name, histo_title, bins_var, bins_pt, bins_eta):

    histos = {}

    h1d = ROOT.TH1D(f"h_{histo_name}", histo_title, len(bins_var)-1, bins_var)
    h2d = ROOT.TH2D(f"h_{histo_name}_2d", f"{histo_title} ", len(bins_pt)-1, bins_pt, len(bins_eta)-1, bins_eta)
    histos.update({f"{histo_name}" : h1d})
    histos.update({f"{histo_name}_2d" : h2d})

    return histos

###############################################################################

def fill_res_histograms(res, hist_dict, bin_dict):
    """
    """
    for bin_key, [_, bin_pt, bin_eta] in bin_dict.items():

        eff, deff = res.getEff(bin_key)
        effMC, deffMC = res.getEffMC(bin_key)

        if isinstance(bin_eta,list) or isinstance(bin_eta,list):
            sys.exit("ERROR: binning not correcty defined, use get_mergedbins_bounds=False in dictionary generation")

        if "efficiency" in hist_dict.keys():
            hist_dict["efficiency"].Fill(                             eff )
            hist_dict["efficiency_2d"].SetBinContent(bin_pt, bin_eta, eff )
        if "rel_err_efficiency" in hist_dict.keys():
            hist_dict["rel_err_efficiency"].Fill(                             deff/eff ) 
            hist_dict["rel_err_efficiency_2d"].SetBinContent(bin_pt, bin_eta, deff/eff )
        if "efficiency_MC" in hist_dict.keys():
            hist_dict["efficiency_MC"].Fill(                             effMC )
            hist_dict["efficiency_MC_2d"].SetBinContent(bin_pt, bin_eta, effMC )
        if "rel_err_efficiency_MC" in hist_dict.keys():
            hist_dict["rel_err_efficiency_MC"].Fill(                             deffMC/effMC ) 
            hist_dict["rel_err_efficiency_MC_2d"].SetBinContent(bin_pt, bin_eta, deffMC/effMC )
        if "sf" in hist_dict.keys():
            hist_dict["sf"].Fill(                             eff/effMC )
            hist_dict["sf_2d"].SetBinContent(bin_pt, bin_eta, eff/effMC )

###############################################################################


def fill_resCmp_histograms(res_bmark, res_test, hist_dict, bin_dict, isPseudodata=False):
    """
    """
    for bin_key, [_, bin_pt, bin_eta] in bin_dict.items():

        if isinstance(bin_eta,list) or isinstance(bin_eta,list):
            sys.exit("ERROR: binning not correcty defined, use get_mergedbins_bounds=False in dictionary generation")

        eff_test, deff_test = res_test.getEff(bin_key)

        if isPseudodata is False:
            eff_bmark, deff_bmark = res_bmark.getEff(bin_key)  
        else:
            eff_bmark, deff_bmark = res_test.getEffMC(bin_key)
        
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        delta       =  eff_test - eff_bmark
        delta_error =  deff_test - deff_bmark
        err_ref  =  deff_bmark
        eff_ref  =  eff_bmark
        if isPseudodata:
            delta = -delta
            delta_error = -delta_error
            err_ref = deff_test
            eff_ref = eff_test
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        pull = delta/(deff_bmark**2 + deff_test**2)**0.5
        pull_ref = delta/err_ref
        rm1 = delta/eff_ref
        ratio_error = delta_error/err_ref
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if "delta" in hist_dict.keys():
            hist_dict["delta"].Fill(                             delta )
            hist_dict["delta_2d"].SetBinContent(bin_pt, bin_eta, delta )
        if "delta_err" in hist_dict.keys():
            hist_dict["delta_err"].Fill(                             delta_error ) 
            hist_dict["delta_err_2d"].SetBinContent(bin_pt, bin_eta, delta_error )
        if "pull" in hist_dict.keys():
            hist_dict["pull"].Fill(                              pull )
            hist_dict["pull_2d"].SetBinContent(bin_pt, bin_eta,  pull )
        if "pull_ref" in hist_dict.keys():
            hist_dict["pull_ref"].Fill(                             pull_ref )
            hist_dict["pull_ref_2d"].SetBinContent(bin_pt, bin_eta, pull_ref )
        if "rm1" in hist_dict.keys():
            hist_dict["rm1"].Fill(                             rm1 )
            hist_dict["rm1_2d"].SetBinContent(bin_pt, bin_eta, rm1 )
        if "ratio_err" in hist_dict.keys():
            hist_dict["ratio_err"].Fill(                             ratio_error )
            hist_dict["ratio_err_2d"].SetBinContent(bin_pt, bin_eta, ratio_error )    
    
        
###############################################################################
###############################################################################
