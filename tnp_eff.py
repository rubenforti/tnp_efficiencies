"""
"""
import ROOT
import time
import json
import sys
import os
from utilities.base_library import bin_dictionary
from utilities.dataset_utils import ws_init, gen_import_dictionary
from make_fits import runFits, runParallelFits


t0 = time.time()

ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

ROOT.Math.MinimizerOptions.SetDefaultMinimizer("Minuit2")

gen_res_folder = "/scratch/rforti/tnp_efficiencies_results"


# -----------------------------------------------------------------------------
#  GENERAL SETTINGS
# ------------------
type_eff = "tracking"
type_analysis = "indep"
charge_selection = ["plus", "minus"]

folder = gen_res_folder+f"/tracking/prefitSS_fail_new"
ws_filename = folder+"/ws_tracking_indep_prefitSS_fail.root"
'''

folder = "."
ws_filename = folder+"/ws_prova.root"

'''

generate_datasets = False

use_extended_sig_template_fail = False

fit_settings = "custom_run3"

binning_pt, binning_eta, binning_mass = "pt_tracking", "eta", "mass_50_130"

mergedbins_bkg, binning_pt_bkg, binning_eta_bkg = False, "pt_12bins", "eta_16bins"

load_bkg_datasets = True
bkg_categories = ["bkg_WW", "bkg_WZ", "bkg_ZZ", "bkg_TTFullyleptonic", "bkg_Ztautau",
                  "bkg_WplusJets", "bkg_WminusJets", "bkg_Zjets",
                  "bkg_SameCharge"]
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
           "check": f"{folder}/fit_plots/check",
           "prefit": f"{folder}/fit_plots/prefit_bkg"}


if mergedbins_bkg and (binning_pt != "pt" or binning_eta != "eta"):
    sys.exit("ERROR: Evaluation of background in merged bins for its comparison on data is allowed only wrt reco-bins of pt and eta for data")

if fit_on_pseudodata: load_bkg_datasets = True


# -----------------------------------------------------------------------------------------------------------
#  DATASET GENERATION
# --------------------
if generate_datasets:

    if not os.path.exists(folder): os.makedirs(folder)

    datasets_folder = "../steve_hists_tmp" 
    import_categories = ["data", "mc"]

    # Not used, MC is always scaled
    if load_bkg_datasets: 
        scale_MC = True,
        if not mergedbins_bkg: import_categories += bkg_categories
    else:
        scale_MC = False
    

    import_dictionary = gen_import_dictionary(datasets_folder, type_eff, import_categories,
                                              ch_set=charge_selection, scale_MC=True, 
                                              add_SS_mc=import_mc_samesign, add_SS_bkg=import_bkg_samesign)

    ws = ws_init(import_dictionary, type_analysis, binning_pt, binning_eta, binning_mass, 
                 lightMode_bkg=True, fail_template_with_all_SA=use_extended_sig_template_fail)

    ws.writeToFile(ws_filename)


    if load_bkg_datasets and mergedbins_bkg: 

        import_dict_bkg = gen_import_dictionary(datasets_folder, type_eff, bkg_categories,
                                                ch_set=charge_selection, add_SS_bkg=import_bkg_samesign)

        if mergedbins_bkg is False:
            binning_pt_bkg, binning_eta_bkg = binning_pt, binning_eta

        lightBkg = True if not fit_on_pseudodata else False
        
        ws = ws_init(import_dict_bkg, type_analysis, binning_pt_bkg, binning_eta_bkg, 
                     binning_mass, import_existing_ws=True, existing_ws_filename=ws_filename, 
                     lightMode_bkg=lightBkg, altBinning_bkg=mergedbins_bkg,
                     fail_template_with_all_SA=use_extended_sig_template_fail)
        ws.writeToFile(ws_filename)

# -----------------------------------------------------------------------------------------------------------
# FIT SETTINGS
# -------------

if fit_settings == "legacy":
    with open(f"legacy_fit_settings.json") as file: fit_settings_json = json.load(file)
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

fit_settings["type_analysis"] = type_analysis
fit_settings["refit_nobkg"] = refit_nobkg
fit_settings["fit_on_pseudodata"] = fit_on_pseudodata
fit_settings["import_mc_samesign"] = import_mc_samesign
fit_settings["fit_verb"] = fit_verb
fit_settings["useMinos"] = useMinos

if load_bkg_datasets or fit_on_pseudodata: 
    fit_settings["bkg_categories"] = bkg_categories
    fit_settings["import_bkg_samesign"] = import_bkg_samesign



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
