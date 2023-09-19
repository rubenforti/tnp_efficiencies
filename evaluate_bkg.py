"""
"""
import ROOT
import time
from utilities.base_library import lumi_factors, binning, bin_dictionary
from utilities.dataset_utils import ws_init
from utilities.bkg_utils import gen_bkg_2d_distrib, bkg_mass_distribution, show_negweighted_bins

t0 = time.time()

ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

# -----------------------------------------------------------------------------
#  GENERAL SETTINGS
# ------------------

type_eff = "iso"
type_analysis = "indep"
bkg_categories = ["WW", "WZ", "ZZ", "TTSemileptonic", "Ztautau", "SameCharge"]
# bkg_categories = ["Ztautau"]

binning_pt = "pt"
binning_eta = "eta"
binning_mass = "mass_60_120"


ws_filename = "root_files/ws_bkg_prova.root"
generate_datasets = False
local_datasets = False

figpath = "bkg_studies/iso/negweights_bkg"
filepath = "bkg_studies/iso"

negweights_eval = False
binnings_list = [["pt", "eta"], ["pt", "eta_24bins"], ["pt", "eta_16bins"], 
                 ["pt", "eta_8bins"], ["pt_12bins", "eta"], ["pt_9bins", "eta"],
                 ["pt_6bins", "eta"], ["pt_12bins", "eta_24bins"], ["pt_9bins", "eta_16bins"],
                 ["pt_12bins", "eta_16bins"], ["pt_9bins", "eta_24bins"]]


minv_plots = {
    "flag" : False,
    "plot_on_data" : False,
    "plot_fit_bkgpdf" : False,
    "plot_on_sig" : True,
    "compare_bkgfrac" : False,
    "logscale" : True,
}
plot_bkg_distrib = {
    "flag" : True,
    "norm_on_data" : False,
    "norm_on_sig" : False,
    "norm_tot_bkg" : False,
}

for tp in ["data", "sig"]:
    minv_plots.update({f"plot_on_{tp}" : minv_plots[f"plot_on_{tp}"]*minv_plots["flag"]})
    plot_bkg_distrib.update({f"norm_on_{tp}" : plot_bkg_distrib[f"norm_on_{tp}"]*plot_bkg_distrib["flag"]})


# -----------------------------------------------------------------------------------------------------------
#  DATASET GENERATION
# --------------------
if generate_datasets is True:

    if local_datasets:
        filename_data = f"root_files/datasets/tnp_{type_eff}_data_vertexWeights1_oscharge1.root"
        filename_mc = f"root_files/datasets/tnp_{type_eff}_mc_vertexWeights1_oscharge1.root"
        dirname_bkg = "root_files/datasets/bkg"
    else:
        filename_data = f"/scratchnvme/wmass/Steve_root_files/Standard_SF_files/tnp_{type_eff}_data_vertexWeights1_oscharge1.root"
        filename_mc = f"/scratchnvme/wmass/Steve_root_files/Standard_SF_files/tnp_{type_eff}_mc_vertexWeights1_oscharge1.root"
        dirname_bkg = "/scratchnvme/rajarshi/Bkg_TNP_3D_Histograms/OS"
    
    dirname_bkg = "root_files/datasets/bkg"

    lumi_scales = lumi_factors(type_eff, bkg_categories)
    sig_lumi_scale = lumi_scales.pop("Zmumu")

    bkg_filenames = {}
    if "SameCharge" in bkg_categories:
        bkg_categories.remove("SameCharge")
        bkg_filenames.update({"SameCharge" : 
            f"{dirname_bkg}/tnp_{type_eff}_data_vertexWeights1_oscharge0.root"})
    [bkg_filenames.update({cat : 
        f"{dirname_bkg}/tnp_{type_eff}_{cat}_vertexWeights1_oscharge1.root"}) for cat in bkg_categories]
    
    
    import_dictionary = {"bkg" : {"filenames":bkg_filenames, "lumi_scales":lumi_scales} }
    if minv_plots["plot_on_data"] or plot_bkg_distrib["norm_on_data"]:
        import_dictionary.update({"data" : filename_data})
    if minv_plots["plot_on_sig"] or plot_bkg_distrib["norm_on_sig"]:
        import_dictionary.update({"mc" : {"filename":filename_mc, "lumi_scale":sig_lumi_scale}})

    bin_dict = bin_dictionary(binning_pt, binning_eta)
    
    ws = ws_init(import_dictionary, type_analysis, binning_pt, binning_eta, binning_mass)
    ws.writeToFile(ws_filename)


# -----------------------------------------------------------------------------------------------------------
#  BKG EVALUATION
# ----------------

if negweights_eval: 
    for binning_pt, binning_eta in binnings_list:
        show_negweighted_bins(ws_filename, bkg_categories, binning_pt, binning_eta, filepath=filepath) 

if minv_plots["flag"]: 
    bkg_mass_distribution(type_eff, ws_filename, bkg_categories, binning_pt, binning_eta,
                          plot_on_data=minv_plots["plot_on_data"], 
                          plot_fit_bkgpdf=minv_plots["plot_fit_bkgpdf"],
                          plot_on_signal=minv_plots["plot_on_sig"], 
                          compare_bkgfrac=minv_plots["compare_bkgfrac"],
                          logscale=minv_plots["logscale"], 
                          figpath=figpath)

if plot_bkg_distrib["flag"]:
    gen_bkg_2d_distrib(ws_filename, bkg_categories, binning_pt, binning_eta, 
                       norm_data=plot_bkg_distrib["norm_on_data"],
                       norm_sig=plot_bkg_distrib["norm_on_sig"],
                       norm_tot_bkg=plot_bkg_distrib["norm_tot_bkg"], 
                       filepath=filepath)








