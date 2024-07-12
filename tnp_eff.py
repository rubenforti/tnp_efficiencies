"""
"""
import os
import sys
import time
import ROOT
import utilities.base_lib as base_lib
from utilities.binning_utils import bin_dictionary
from utilities.bkg_utils import base_parser_bkg
from utilities.fit_utils import base_parser_fit, finalize_fit_parsing, checkImport_custom_pdf, printFitStatus, doSingleFit
from utilities.dataset_utils import ws_init
from fitter import IndepFitter, SimFitter
from multiprocessing import Pool
from itertools import repeat


t0 = time.time()

ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True
ROOT.Math.MinimizerOptions.SetDefaultMinimizer("Minuit2")


parser = base_lib.base_parser()
parser = base_parser_bkg(parser)
parser = base_parser_fit(parser)
args = parser.parse_args()

base_lib.control_parsing(args)
finalize_fit_parsing(args)



if "sim" in args.type_analysis:
    sys.exit("ERROR: Simultaneous fits are implemented but not tested since long time, please check the code before running")


# -------------------------
#  SUMMARY OF THE SETTINGS
# -------------------------
print("\n")
print("--"+('------------'*4)+"--")
print("|-"+('------------'*4)+"-|")
print("|         ~ Summary of the fit settings ~          |")
print("|-"+('------------'*4)+"-|")
print("--"+('------------'*4)+"--")
for key, value in vars(args).items():
    if "fit_settings" in key: continue
    print("| ", key, " ", value, f"\n| {'------------'*4}")
print("--"+('------------'*4)+"--")
print(f"| {'------------'*4}")
for key, value in args.par_fit_settings.items():
    print("| ", key, " ", value, f"\n| {'------------'*4}")
print("--"+('------------'*4)+"--")
print("\n\n")

if not args.auto_run:
    confirm = input("Do you confirm the settings and start the fit? (y/n) ")
    if confirm != "y": sys.exit("Exiting...")


# --------------------
#  DATASET GENERATION
# --------------------
if args.generate_ws is True:
    # Two-step routine: first generate the workspace for data and the signal
    # template; then, if requested, import the background template in the same
    # workspace (with different binning if specified)

    # Step 1
    ws = ws_init(args.input, ["data", "mc"], args.eff, args.binning_pt, args.binning_eta, args.binning_mass,
                 type_analysis=args.type_analysis, ch_set=args.charge_selection, do_OS_tracking=args.OS_tracking,
                 fail_template_with_all_SA=args.all_SA_fail, add_SS_mc=args.import_mc_SS)
    ws.writeToFile(args.ws_filename)

    # Step 2
    if args.importBkg:
        
        using_mergedbins_bkg = (bool(args.mergedbins_bkg[0]) or bool(args.mergedbins_bkg[1]))
        binning_pt_bkg, binning_eta_bkg = args.binning_pt, args.binning_eta if not using_mergedbins_bkg else args.mergedbins_bkg

        ws = ws_init(args.input, args.bkg_categories, args.eff, binning_pt_bkg, binning_eta_bkg, args.binning_mass,
                     type_analysis=args.type_analysis, ch_set=args.charge_selection, do_OS_tracking=args.OS_tracking,
                     add_SS_bkg=args.import_bkg_SS, lightMode_bkg=True, altBinning_bkg=using_mergedbins_bkg,
                     import_existing_ws=True, existing_ws_filename=args.ws_filename)
        ws.writeToFile(args.ws_filename)


# sys.exit("Exiting...")


# --------------
#  RUNNING FITS
# --------------

checkImport_custom_pdf(args.par_fit_settings["bkg_model"])

file_ws = ROOT.TFile(args.ws_filename)
ws = file_ws.Get("w")

if not args.no_saveFigs:
    for path in args.figpath.values(): 
        if not os.path.exists(path): os.makedirs(path)

prob_bins = []

for bin_key in bin_dictionary(args.binning_pt, args.binning_eta).keys():

    # Routine for fits performed in series, the parallel version is in the next block

    if args.parallel: continue

    # if bin_key != "[24.0to35.0][-2.3to-2.2]": continue
    
    if args.type_analysis == "indep":
        dict_flags = ["pass", "fail"]
        fitter = IndepFitter(bin_key, args.fit_settings)            
    elif args.type_analysis in ["sim", "sim_sf"]:
        dict_flags = ["sim"]
        fitter = SimFitter(bin_key, args.fit_settings)
    else:
        sys.exit("ERROR: wrong analysis type indicated")
    
    fitter.manageFit(ws)


    res = {flag : fitter.res_obj[flag] for flag in dict_flags}

    if not args.no_saveFigs and fitter.existingFit is False: 
        fitter.saveFig(ws, args.figpath)
        for flag in ["pass", "fail"]:
            if "prefit" in fitter.settings["bkg_model"][flag]: 
                fitter.saveFig_prefit(flag, args.figpath["prefit"])

    if not args.no_importPdfs: fitter.importFitObjects(ws)

    if fitter.bin_status is False: 
        prob_bins.append(bin_key)
    printFitStatus(args.type_analysis, bin_key, res, fitter.bin_status)



if args.parallel:

    if args.type_analysis == "indep":
        dict_flags = ["pass", "fail"]
        fitters = [IndepFitter(bin_key, args.fit_settings) for bin_key in bin_dictionary(args.binning_pt, args.binning_eta).keys()]
    elif args.type_analysis in ["sim", "sim_sf"]:
        dict_flags = ["sim"]
        fitters = [SimFitter(bin_key, args.fit_settings) for bin_key in bin_dictionary(args.binning_pt, args.binning_eta).keys()]
    else:
        print("ERROR: wrong analysis type indicated")
        sys.exit()


    def doSingleFit(fitter): 
        fitter.manageFit(ws)
        if fitter.bin_status is False: 
            prob_bins.append(fitter.bin_key)

    pool = Pool(processes=8)
    pool.map(doSingleFit, fitters)

    print("FITS DONE")
    
    for fitter_out in fitters:
        print(f"Analyzing bin {fitter_out.bin_key}")
        print(fitter_out.bin_status)
        
        if not args.no_importPdfs: fitter_out.importFitObjects(ws)

        if not args.no_saveFigs and fitter_out.existingFit is False: 
            fitter_out.saveFig(ws)
            for flag in ["pass", "fail"]:
                if "prefit" in fitter_out.settings["bkg_model"][flag]: 
                    fitter_out.saveFig_prefit(flag)

        #if fitter_out.bin_status is False: prob_bins.append(fitter_out.bin_key)
        # printFitStatus(args.type_analysis, fitter_out.bin_key, res, fitter_out.bin_status, prob_bins)


print(f"NUM of problematic bins = {len(prob_bins)}")
print(prob_bins)

ws.writeToFile(args.ws_filename)

# -----------------------------------------------------------------------------------------------------------

t1 = time.time()

print(f"Time spent = {(t1-t0)/60} min")
