"""
"""
import ROOT
import time
# from copy import copy
# from utilities.base_library import lumi_factors, binning, bin_dictionary
from utilities.dataset_utils import ws_init, gen_import_dictionary
from utilities.bkg_utils import bkg_2d_distrib, bkg_mass_distribution, show_negweighted_bins

t0 = time.time()

ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

gen_res_folder = "/scratch/rforti/tnp_efficiencies_results"


# -----------------------------------------------------------------------------
#  GENERAL SETTINGS
# ------------------

type_eff = "tracking"
type_analysis = "indep"
bkg_categories = ["bkg_WW", "bkg_WZ", "bkg_ZZ", 
                  "bkg_TTSemileptonic", "bkg_TTFullyleptonic", "bkg_Ztautau",
                  "bkg_WplusJets", "bkg_WminusJets",                  
                  "bkg_SameCharge"]


binning_pt = "pt_tracking"
binning_eta = "eta"
binning_mass = "mass_50_130"

study_SS_bkg = True

folder = gen_res_folder+"/tracking/bkg_figs/SS"

ws_filename = folder+"/ws_tracking_bkg_SS.root"

generate_datasets = False       

figpath = folder
filepath = folder

negweights_eval = False
binnings_list = [["pt", "eta"]]
                 

minv_plots = {
    "flag" : True,
    "plot_on_data" : True,
    "plot_fit_bkgpdf" : False,
    "plot_on_sig" : True,
    "compare_bkgfrac" : False,
    "logscale" : "hybrid",
}
plot_bkg_distrib = {
    "flag" : False,
    "norm_on_data" : False,
    "norm_on_sig" : False,
    "norm_tot_bkg" : False,
    "plot_projected" : False,
}

for tp in ["data", "sig"]:
    minv_plots.update({f"plot_on_{tp}" : bool(minv_plots[f"plot_on_{tp}"]*minv_plots["flag"])})
    plot_bkg_distrib.update({f"norm_on_{tp}" : bool(plot_bkg_distrib[f"norm_on_{tp}"]*plot_bkg_distrib["flag"])})


# -----------------------------------------------------------------------------------------------------------
#  DATASET GENERATION
# --------------------

if generate_datasets is True:

    # base_folder = "/scratchnvme/rajarshi/Latest_3D_Steve_Histograms_22_Sep_2023"
    base_folder = "datasets"

    import_categories = ["data", "mc"] + bkg_categories
    
    import_dictionary = gen_import_dictionary(base_folder, type_eff, import_categories,
                                              ch_set=["plus", "minus"], scale_MC=True, 
                                              add_SS_mc=True, add_SS_bkg=study_SS_bkg)
    
    for k, v in import_dictionary.items():
        print(k, v)


    
    ws = ws_init(import_dictionary, type_analysis, binning_pt, binning_eta, binning_mass)
    ws.writeToFile(ws_filename)
    


# -----------------------------------------------------------------------------------------------------------
#  BKG EVALUATION
# ----------------

if negweights_eval: 
    for binning_pt, binning_eta in binnings_list:
        show_negweighted_bins(ws_filename, bkg_categories, binning_pt, binning_eta, filepath=filepath) 

if minv_plots["flag"]: 
    print("Plotting minv distributions")
    bkg_mass_distribution(type_eff, ws_filename, bkg_categories, binning_pt, binning_eta,
                          study_SS_bkg=study_SS_bkg,
                          plot_on_data=minv_plots["plot_on_data"], 
                          plot_fit_bkgpdf=minv_plots["plot_fit_bkgpdf"],
                          plot_on_signal=minv_plots["plot_on_sig"], 
                          compare_bkgfrac=minv_plots["compare_bkgfrac"],
                          logscale=minv_plots["logscale"],
                          figpath=figpath)

if plot_bkg_distrib["flag"]:
    print("Plotting bkg 2D distributions")
    bkg_2d_distrib(ws_filename, bkg_categories, binning_pt, binning_eta,
                       study_SS_bkg=study_SS_bkg,
                       norm_data=plot_bkg_distrib["norm_on_data"],
                       norm_sig=plot_bkg_distrib["norm_on_sig"],
                       norm_tot_bkg=plot_bkg_distrib["norm_tot_bkg"], 
                       plot_projected=plot_bkg_distrib["plot_projected"],
                       filepath=filepath)





