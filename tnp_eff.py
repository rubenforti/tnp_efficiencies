"""
"""
import ROOT
import time
import json
import sys
from utilities.base_library import bin_dictionary
from utilities.dataset_utils import ws_init, gen_import_dictionary
from make_fits import runFits, runParallelFits


t0 = time.time()

ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

ROOT.Math.MinimizerOptions.SetDefaultMinimizer("Minuit2")

gen_res_folder = "/scratchnvme/rforti/tnp_efficiencies_results"


# -----------------------------------------------------------------------------
#  GENERAL SETTINGS
# ------------------
type_eff = "tracking"
type_analysis = "indep"
charge_selection = ["plus", "minus"]

folder = gen_res_folder+"/results_tracking/benchmark_tracking"
ws_filename = folder+"/ws_tracking_indep_benchmark.root"

generate_datasets = False

fit_settings = "custom_run3"

binning_pt, binning_eta, binning_mass = "pt_tracking", "eta", "mass_50_130"

mergedbins_bkg, binning_pt_bkg, binning_eta_bkg = False, "pt_12bins", "eta_16bins"

load_bkg_datasets = False
bkg_categories = ["bkg_WW", "bkg_WZ", "bkg_ZZ", "bkg_TTFullyleptonic", "bkg_Ztautau", "bkg_SameCharge"]
import_bkg_samesign = False
import_mc_samesign = False

fit_on_pseudodata = False

fit_verb = -1

parallel_fits = False

refit_nobkg = True

useMinos = False

import_pdfs = True

savefigs = True

figpath = {"good": f"{folder}/fit_plots", 
           "check": f"{folder}/fit_plots/check"} 

if mergedbins_bkg and (binning_pt != "pt" or binning_eta != "eta"):
    sys.exit("ERROR: Evaluation of background in merged bins for its comparison on data is allowed only wrt reco-bins of pt and eta for data")


# -----------------------------------------------------------------------------------------------------------
#  DATASET GENERATION
# --------------------
if generate_datasets:

    datasets_folder = "/scratchnvme/rajarshi/Latest_3D_Steve_Histograms_22_Sep_2023"

    import_categories = ["data", "mc"]
    scale_MC= True if bool(load_bkg_datasets or fit_on_pseudodata) else False

    import_dictionary = gen_import_dictionary(datasets_folder, type_eff, import_categories,
                                              ch_set=charge_selection, scale_MC=scale_MC, 
                                              add_SS_mc=import_mc_samesign)

    ws = ws_init(import_dictionary, type_analysis, binning_pt, binning_eta, binning_mass)
    ws.writeToFile(ws_filename)

    if load_bkg_datasets or fit_on_pseudodata: 

        import_dict_bkg = gen_import_dictionary(datasets_folder, type_eff, bkg_categories,
                                                ch_set=charge_selection, add_SS_bkg=import_bkg_samesign)

        if mergedbins_bkg is False:
            binning_pt_bkg, binning_eta_bkg = binning_pt, binning_eta
        
        ws_bkg = ws_init(import_dict_bkg, type_analysis, binning_pt_bkg, binning_eta_bkg, 
                             binning_mass, import_existing_ws=True, existing_ws_filename=ws_filename, 
                             altBinning_bkg=mergedbins_bkg)
        ws_bkg.writeToFile(ws_filename)

# -----------------------------------------------------------------------------------------------------------
# FIT SETTINGS
# -------------

if fit_settings == "default":
    with open(f"default_fit_settings.json") as file: fit_settings_json = json.load(file)
    if type_eff in ["idip", "trigger", "iso"]:
        fit_settings = fit_settings_json["idip_trig_iso"]
    else:
        fit_settings = fit_settings_json[type_eff]

elif "custom" in fit_settings:
    with open(f"custom_fit_settings.json") as file: fit_settings_json = json.load(file)
    run_key = fit_settings.split("_")[1]
    nRUN = int(run_key.replace("run", ""))
    fit_settings = fit_settings_json["run1"]
    [fit_settings.update(fit_settings_json[run_key]) for run_idx in range(2, nRUN+1)]

else:
    sys.exit("ERROR: wrong fit settings indicated")

fit_settings.update({"type_analysis" : type_analysis})
fit_settings.update({"refit_nobkg" : refit_nobkg})
fit_settings.update({"fit_on_pseudodata" : fit_on_pseudodata})
fit_settings.update({"fit_verb" : fit_verb})
fit_settings.update({"useMinos" : useMinos})

if load_bkg_datasets or fit_on_pseudodata: fit_settings.update({"bkg_categories" : bkg_categories})


# -----------------------------------------------------------------------------------------------------------
#  RUNNING FITS
# --------------

if parallel_fits is False:
    runFits(ws_filename, bin_dictionary(binning_pt, binning_eta), fit_settings, 
            import_pdfs=import_pdfs, savefigs=savefigs, figpath=figpath)
else:
    runParallelFits(ws_filename, bin_dictionary(binning_pt, binning_eta), fit_settings,
                    import_pdfs=import_pdfs, savefigs=savefigs, figpath=figpath)


# print(fit_settings)


# -----------------------------------------------------------------------------------------------------------

t1 = time.time()

print(f"TEMPO = {(t1-t0)/60} min")
