"""
"""
import sys, time
import ROOT
from utilities.dataset_utils import ws_init, gen_import_dictionary
from utilities.bkg_utils import bkg_2d_distrib, bkg_mass_distribution, show_negweighted_bins


t0 = time.time()

ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

gen_res_folder = "/scratch/rforti/tnp_efficiencies_results"


# -----------------------------------------------------------------------------
#  GENERAL SETTINGS
# ------------------

input_folder = "/scratch/rforti/steve_histograms_2016"

type_eff = "reco"
type_analysis = "indep"
bkg_categories = ["bkg_SameCharge",
                  "bkg_WW", "bkg_WZ", "bkg_ZZ", 
                  "bkg_TTSemileptonic", "bkg_TTFullyleptonic", "bkg_Ztautau",
                  "bkg_WplusJets", "bkg_WminusJets", "bkg_Zjets"
                  ]


binning_pt = "pt_singlebin"
binning_eta = "eta_singlebin"
binning_mass = "mass_60_120"

study_SS_bkg = False 

charge_flag = "OS" if not study_SS_bkg else "SS"

folder = "/scratch/rforti/tnp_efficiencies_results/recoplus/prova"
ws_filename = folder+f"/ws_reco_singlebin.root"

# folder = gen_res_folder+f"/tracking/bkg_figs"

#ws_filename = folder+f"/ws_tracking_bkg_{charge_flag}_singlebin.root"
'''

folder = "zjets_tests/prova_bkg_distr/OS"
ws_filename = "zjets_tests/ws_prova.root"
'''

generate_datasets = True       

figpath = folder
filepath = folder

negweights_eval = False
binnings_list = [["pt", "eta"]]
                 

minv_plots = {
    "flag" : True,
    "plot_on_data" : False,
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


if study_SS_bkg:
    bkg_categories = [f"{cat}_SS" if (cat!="bkg_SameCharge") else cat for cat in bkg_categories]


# -----------------------------------------------------------------------------------------------------------
#  DATASET GENERATION
# --------------------

if generate_datasets is True:

    # base_folder = "/scratchnvme/rajarshi/Latest_3D_Steve_Histograms_22_Sep_2023"
    base_folder = f"/scratch/rforti/steve_histograms_2016/{type_eff}"
    
    import_categories = bkg_categories
    if minv_plots["flag"] and minv_plots["plot_on_data"]:
        import_categories += ["data"]
    if minv_plots["flag"] and minv_plots["plot_on_sig"]:
        import_categories += ["mc"]
    
    import_dictionary = gen_import_dictionary(base_folder, type_eff, import_categories,
                                              ch_set=["plus", "minus"],
                                              do_OS_tracking=True,
                                              add_SS_mc=study_SS_bkg)
    
    '''
    for k, v in import_dictionary.items():
        print(k, " ", v)
        print(" ")
    sys.exit()  
    '''

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
                          is_SS_bkg=study_SS_bkg,
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





