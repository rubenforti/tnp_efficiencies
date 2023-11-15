"""
"""
import ROOT
import time
import json
import sys
from copy import copy
from multiprocessing import Pool
from utilities.base_library import bin_dictionary, lumi_factors
from utilities.dataset_utils import ws_init
from make_fits import runFits, runParallelFits

t0 = time.time()

ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

ROOT.Math.MinimizerOptions.SetDefaultMinimizer("Minuit2")


# -----------------------------------------------------------------------------
#  GENERAL SETTINGS
# ------------------
type_eff = "trackingplus"
type_analysis = "indep"
name_analysis = ""

main_folder = "root_files"


local_datasets = False
load_McBkg = False

generate_datasets = False
default_fit_settings = False
nRUN = 3

binning_pt, binning_eta, binning_mass = "pt_tracking", "eta", "mass_50_130"
mergedbins_bkg = False
binning_pt_bkg, binning_eta_bkg = "pt_12bins", "eta_16bins"
bkg_categories = ["WW", "WZ", "ZZ", "TTFullyleptonic", "Ztautau", "SameCharge"]

fit_on_pseudodata = False

fit_verb = -1
parallel_fits = False
refit_nobkg = True
useMinos = False

#folder = f"{main_folder}/{type_eff}_{type_analysis}_r628"
#workspace_name = f"{folder}/ws_{type_eff}_{type_analysis}.root"
folder = f"trackingplus_res"
workspace_name = f"trackingplus_res/ws_trackingplus.root"

import_pdfs = True
savefigs = True



figpath = {"good": f"{folder}/fit_plots", 
           "check": f"{folder}/fit_plots/check"} 
'''
figpath = {"good": folder, 
           "check": folder}
'''
# figpath = {"check" : "."}

if mergedbins_bkg and (binning_pt != "pt" or binning_eta != "eta"):
    print("ERROR: Evaluation of background in merged bins for its comparison on data is allowed only wrt reco-bins of pt and eta for data")
    sys.exit()

if parallel_fits is True and type_analysis != "indep":
    print("ERROR: parallelization of fits is allowed only for independent analysis")
    sys.exit()

# -----------------------------------------------------------------------------------------------------------
#  DATASET GENERATION
# --------------------
if generate_datasets:

    if load_McBkg or fit_on_pseudodata:
        lumi_scales = lumi_factors(type_eff, bkg_categories)
        lumi_scale_signal = lumi_scales.pop("Zmumu")
        print("Lumi scale sig:  ", lumi_scale_signal)
    else: 
        lumi_scale_signal = 1

    if local_datasets:
        filename_data = f"root_files/datasets/tnp_{type_eff}_data_vertexWeights1_oscharge1.root"
        filename_mc = f"root_files/datasets/tnp_{type_eff}_mc_vertexWeights1_oscharge1.root"
        dirname_bkg = "root_files/datasets/bkg"
    else:
        filename_data = f"/scratchnvme/rajarshi/Latest_3D_Steve_Histograms_22_Sep_2023/OS/tnp_{type_eff}_data_vertexWeights1_oscharge1.root"
        filename_mc = f"/scratchnvme/rajarshi/Latest_3D_Steve_Histograms_22_Sep_2023/OS/tnp_{type_eff}_mc_vertexWeights1_oscharge1.root"
        dirname_bkg = "/scratchnvme/rajarshi/Latest_3D_Steve_Histograms_22_Sep_2023/OS"

    import_dictionary = {
        "data" : filename_data, 
        "mc" : {"filename": filename_mc, "lumi_scale" : lumi_scale_signal}}
   
    ws = ws_init(import_dictionary, type_analysis, binning_pt, binning_eta, binning_mass)
    ws.writeToFile(workspace_name)

    if load_McBkg or fit_on_pseudodata:
        bkg_filenames = {}

        load_bkgCat = copy(bkg_categories)
        if "SameCharge" in bkg_categories:
            load_bkgCat.remove("SameCharge")
            bkg_filenames.update({"SameCharge" : 
            f"{dirname_bkg}/../SS/tnp_{type_eff}_data_vertexWeights1_oscharge0.root"})
        [bkg_filenames.update({cat : 
            f"{dirname_bkg}/tnp_{type_eff}_{cat}_vertexWeights1_oscharge1.root"}) for cat in load_bkgCat]

        import_dict_bkg = {"bkg" : {"filenames" : bkg_filenames, "lumi_scales" : lumi_scales}}

        if mergedbins_bkg is False:
            ws_bkg = ws_init(import_dict_bkg, type_analysis, binning_pt, binning_eta, 
                             binning_mass, import_existing_ws=True, existing_ws_filename=workspace_name, 
                             altBinning_bkg=False)

        else:
            ws_bkg = ws_init(import_dict_bkg, type_analysis, binning_pt_bkg, binning_eta_bkg, 
                             binning_mass, import_existing_ws=True, existing_ws_filename=workspace_name, 
                             altBinning_bkg=True)

        ws_bkg.writeToFile(workspace_name)


# -----------------------------------------------------------------------------------------------------------
# FIT SETTINGS
# -------------

fit_settings_filename = "default_fit_settings.json" if default_fit_settings else "custom_fit_settings.json"

with open(fit_settings_filename) as file:
    fit_settings_json = json.load(file)

if default_fit_settings:
    if ("idip" in type_eff) or ("trigger" in type_eff) or ("iso" in type_eff):
        fit_settings = fit_settings_json["idip_trig_iso"]
    else:
        for pm_flag in ["plus", "minus"]:
            if pm_flag in type_eff: fit_settings = fit_settings_json[type_eff.replace(pm_flag, "")]

else:
    fit_settings = fit_settings_json["run_1"]
    [fit_settings.update(fit_settings_json[f"run_{run_idx}"]) for run_idx in range(2, nRUN+1)]

fit_settings.update({"type_analysis" : type_analysis})
fit_settings.update({"refit_nobkg" : refit_nobkg})
fit_settings.update({"fit_on_pseudodata" : fit_on_pseudodata})
fit_settings.update({"fit_verb" : fit_verb})
fit_settings.update({"useMinos" : useMinos})

if load_McBkg or fit_on_pseudodata:
    fit_settings.update({"bkg_categories" : bkg_categories})

# -----------------------------------------------------------------------------------------------------------
#  RUNNING FITS
# --------------

if parallel_fits is False:
    runFits(workspace_name, bin_dictionary(binning_pt, binning_eta), fit_settings, 
            import_pdfs=import_pdfs, savefigs=savefigs, figpath=figpath)
else:
    runParallelFits(workspace_name, bin_dictionary(binning_pt, binning_eta), fit_settings,
                    import_pdfs=import_pdfs, savefigs=savefigs, figpath=figpath)


# print(fit_settings)


# -----------------------------------------------------------------------------------------------------------

t1 = time.time()

print(f"TEMPO = {(t1-t0)/60} min")
