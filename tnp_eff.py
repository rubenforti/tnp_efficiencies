"""
"""
import ROOT
import time
import json
import sys
import os
import argparse
from copy import deepcopy as dcp
from utilities.binning_utils import bin_dictionary
from utilities.dataset_utils import ws_init, gen_import_dictionary
from make_fits import runFits, runParallelFits


t0 = time.time()

ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

ROOT.Math.MinimizerOptions.SetDefaultMinimizer("Minuit2")


gen_res_folder = "/scratch/rforti/tnp_efficiencies_results"
# datasets_folder = "../steve_hists_tmp" 
datasets_folder = "/scratch/rforti/steve_histograms_2016/tracking"

parser = argparse.ArgumentParser(description="Run the fits for the tracking efficiency")

parser.add_argument("--gen_folder", default=gen_res_folder,
                    help="General folder where the fits will be saved")

parser.add_argument("--process_name",  help="Name that characterizes the fit strategy")

parser.add_argument("--datasets_folder", default=datasets_folder,
                    help="Folder where the datasets are stored")

parser.add_argument("--ws_flags", nargs="+", default=[],
                    help="Flags to be added to the workspace name")

parser.add_argument("-g", "--generate_ws", action="store_true", help="Generate the workspace")

parser.add_argument("--type_fit", default="indep", help="Type of fit to be performed")

parser.add_argument("--type_eff", default="tracking",
                    choices=["reco", "tracking", "idip", "trigger", "iso"],
                    help="Type of efficiency")

parser.add_argument("-ch", "--charge_selection", default="all",
                    choices=["plus", "minus", "all"], help="Charge selection")

args = parser.parse_args()

gen_res_folder = args.gen_folder

if args.charge_selection != "all":
    eff_print = args.type_eff+args.charge_selection 
else:
    eff_print = args.type_eff

folder = gen_res_folder+"/"+eff_print+"/"+args.process_name
datasets_folder = args.datasets_folder


ws_name_str = f"ws_{args.type_eff}"
for postfix in args.ws_flags: ws_name_str += f"_{postfix}"
ws_name_str += ".root"


fit_settings = {

    "ws_name" : folder+"/"+ws_name_str,

    "type_eff" : args.type_eff,

    "type_analysis" : args.type_fit,

    "generate_datasets" : args.generate_ws,

    "charge_selection" : ["plus", "minus"] if args.charge_selection == "all" else [args.charge_selection],

    "binning_pt" : "pt_reco",
    "binning_eta" : "eta",
    "binning_mass" : "mass_60_120",

    "fitOnlyBkg" : False,

    "fitPseudodata" : False,

    "use_extended_sig_template_fail" : False,

    "par_fit_settings_type" : "custom",

    "load_bkg_datasets" : True,

    "bkg_categories" : ["bkg_WW", "bkg_WZ", "bkg_ZZ", 
                        "bkg_TTFullyleptonic", "bkg_TTSemileptonic",
                        "bkg_WplusJets", "bkg_WminusJets", "bkg_Zjets",
                        "bkg_Ztautau",
                        # "bkg_SameCharge"
                        ],
    
    "import_bkg_samesign" : False,
    "import_mc_samesign" : False,

    "mergedbins_bkg" : False,
    "binning_pt_bkg" : "pt_12bins",
    "binning_eta_bkg" : "eta_16bins",

    "fit_verb" : -1,

    "parallel_fits" : False,

    "refit_nobkg" : True,

    "useMinos" : False,

    "import_pdfs" : True,

    "savefigs" : True, 
    "figpath" : {"good": f"{folder}/fit_plots",
                "check": f"{folder}/fit_plots/check",
                "prefit": f"{folder}/fit_plots/prefit_bkg"}

}

if fit_settings["mergedbins_bkg"] and (fit_settings["binning_pt"] != "pt" or fit_settings["binning_eta"] != "eta"):
    sys.exit("ERROR: Evaluation of background in merged bins for its comparison on data is allowed only wrt reco-bins of pt and eta for data")

if fit_settings["fitPseudodata"] or fit_settings["fitOnlyBkg"]: 
    fit_settings["load_bkg_datasets"] = True


# -----------------------------------------------------------------------------------------------------------
#  SUMMARY OF THE SETTINGS
# -------------------------
print("\n")
print("--"+('------------'*4)+"--")
print("|-"+('------------'*4)+"-|")
print("|         ~ Summary of the fit settings ~          |")
print("|-"+('------------'*4)+"-|")
print("--"+('------------'*4)+"--")
for key, value in fit_settings.items(): print("| ", key, " ", value, f"\n| {'------------'*4}")
print("--"+('------------'*4)+"--")
print("\n\n")

confirm = input("Do you confirm the settings and start the fit? (y/n) ")
if confirm != "y": sys.exit("Exiting...")


# -----------------------------------------------------------------------------------------------------------
#  DATASET GENERATION
# --------------------

if fit_settings["generate_datasets"] and os.path.exists(fit_settings["ws_name"]): 
    confirm_gen = input("Workspace already exists, do you REALLY want to generate a new one from scratch? (y/n) ")
    if confirm_gen != "y": fit_settings["generate_datasets"] = False


if fit_settings["generate_datasets"]:

    if not os.path.exists(folder): os.makedirs(folder)

    import_categories = ["data", "mc"]

    if fit_settings["load_bkg_datasets"] and not fit_settings["mergedbins_bkg"]: 
        import_categories += fit_settings["bkg_categories"]

    import_dictionary = gen_import_dictionary(datasets_folder, fit_settings["type_eff"], import_categories,
                                              ch_set=fit_settings["charge_selection"],
                                              do_OS_tracking=False,
                                              add_SS_mc=fit_settings["import_mc_samesign"], 
                                              add_SS_bkg=fit_settings["import_bkg_samesign"])

    
    
    ws = ws_init(import_dictionary, fit_settings["type_analysis"], 
                 fit_settings["binning_pt"], fit_settings["binning_eta"], fit_settings["binning_mass"], 
                 lightMode_bkg=True, # (True if fit_settings["import_bkg_samesign"] is False else False), 
                 fail_template_with_all_SA=fit_settings["use_extended_sig_template_fail"])

    ws.writeToFile(fit_settings["ws_name"])
    ws.Print()
    
    if fit_settings["load_bkg_datasets"] and fit_settings["mergedbins_bkg"]: 

        import_dict_bkg = gen_import_dictionary(datasets_folder, fit_settings["type_eff"], fit_settings["bkg_categories"],
                                                ch_set=fit_settings["charge_selection"], add_SS_bkg=fit_settings["import_bkg_samesign"])

        if fit_settings["mergedbins_bkg"] is False:
            binning_pt_bkg, binning_eta_bkg = fit_settings["binning_pt"], fit_settings["binning_eta"]
        
        ws = ws_init(import_dict_bkg, fit_settings["type_analysis"], binning_pt_bkg, binning_eta_bkg, 
                     fit_settings["binning_mass"], import_existing_ws=True, existing_ws_filename=fit_settings["ws_name"], 
                     lightMode_bkg=True, altBinning_bkg=fit_settings["mergedbins_bkg"],
                     fail_template_with_all_SA=fit_settings["use_extended_sig_template_fail"])
        ws.writeToFile(fit_settings["ws_name"])


sys.exit("Exiting...")


# -----------------------------------------------------------------------------------------------------------
# FIT PARAMETERS SETTINGS
# ------------------------

if fit_settings["par_fit_settings_type"] == "legacy":
    with open(f"legacy_fit_settings.json") as file: fit_settings_json = json.load(file)
    if fit_settings["type_eff"] in ["idip", "trigger", "iso"]:
        par_fit_settings = fit_settings_json["idip_trig_iso"]
    else:
        par_fit_settings = fit_settings_json[fit_settings["type_eff"]]

elif fit_settings["par_fit_settings_type"] == "custom":
    with open(f"custom_fit_settings.json") as file: fit_settings_json = json.load(file)
    #run_key = fit_settings["par_fit_settings_type"].split("_")[1]
    nRUN = input("Insert the run number: ")
    par_fit_settings = fit_settings_json["run1"]
    [par_fit_settings.update(fit_settings_json["run"+str(run_idx)]) for run_idx in range(2, int(nRUN)+1)]
else:
    sys.exit("ERROR: wrong fit settings indicated")

fit_settings.update(par_fit_settings)

for k,v in fit_settings.items():
    print(k, " ", v)

# -----------------------------------------------------------------------------------------------------------
#  RUNNING FITS
# --------------

if fit_settings["parallel_fits"] is False:
    runFits(fit_settings)
else:
    runParallelFits(fit_settings)


# print(fit_settings)


# -----------------------------------------------------------------------------------------------------------

t1 = time.time()

print(f"TEMPO = {(t1-t0)/60} min")
