"""
"""

import ROOT
import os
import sys
from multiprocessing import Pool
from itertools import repeat
from fitter import IndepFitter, SimFitter
from utilities.dataset_utils import import_pdf_library
from utilities.base_library import bin_dictionary



def printFitStatus(type_analysis, bin_key, res, status):

    if status is True:
        print(f"Bin {bin_key} is OK!\n\n")
    else:
        print(f"Bin {bin_key} has problems!\n")
        if type_analysis == "indep":
            print("****")
            res["pass"].Print()
            res["pass"].correlationMatrix().Print()
            print("****")
            res["fail"].Print()
            res["fail"].correlationMatrix().Print()
            print("****")
            print(res["pass"].status(), res["fail"].status())
            print(res["pass"].covQual(), res["fail"].covQual())
            print(res["pass"].edm(), res["fail"].edm())
            print('\n')
        elif type_analysis == "sim" or type_analysis == "sim_sf":
            print("****")
            res["sim"].Print()
            res["sim"].correlationMatrix().Print()
            print("****")
            print(res["sim"].status())
            print(res["sim"].covQual())
            print(res["sim"].edm())
            print('\n\n')

###############################################################################


def checkImportCustomPdfs(bkg_models):
    """
    """
    dict_classes = {
        "cmsshape" : "RooCMSShape",
        "cmsshape_prefitSS" : "RooCMSShape",
        "cmsshape_new" : "RooCMSShape_mod",
        "CB" : "my_double_CB"
    }
    for flag in ["pass", "fail"]:
        if bkg_models[flag] in dict_classes.keys():
            import_pdf_library(dict_classes[bkg_models[flag]])
            print(f"Imported {dict_classes[bkg_models[flag]]} from pdf_library")

###############################################################################


def doSingleFit(fitter, ws, flags):
    """
    """
    fitter.manageFit(ws)

    res = {}

    for flag in flags: res.update({flag : fitter.res_obj[flag]})

    if fitter.settings["savefigs"] is True and fitter.existingFit is False: 
        fitter.saveFig(ws)
        for flag in ["pass", "fail"]:
            print(fitter.settings["bkg_model"][flag], "prefit" in fitter.settings["bkg_model"][flag])
            if "prefit" in fitter.settings["bkg_model"][flag]: 
                print("\nSaving prefit plots\n")
                fitter.saveFig_prefit(flag)

    return fitter, res, fitter.bin_status

###############################################################################


def runFits(settings):
    """
    """
    path = os.path.dirname(__file__)
    ROOT.gSystem.cd(path)

    checkImportCustomPdfs(settings["bkg_model"])

    file_ws = ROOT.TFile(settings["ws_name"])
    ws = file_ws.Get("w")

    if settings["savefigs"] is True:
        for path in settings["figpath"].values(): 
            if not os.path.exists(path): os.makedirs(path)

    prob_bins = []

    for bin_key in bin_dictionary(settings["binning_pt"], settings["binning_eta"]).keys():
        
        if settings["type_analysis"] == "indep":
            dict_flags = ["pass", "fail"]
            fitter = IndepFitter(bin_key, settings)            
        elif settings["type_analysis"] in ["sim", "sim_sf"]:
            dict_flags = ["sim"]
            fitter = SimFitter(bin_key, settings)
        else:
            print("ERROR: wrong analysis type indicated")
            sys.exit()

        fitter, res, status = doSingleFit(fitter, ws, dict_flags)

        if settings["import_pdfs"]: fitter.importFitObjects(ws)

        if status is False: prob_bins.append(bin_key)
        printFitStatus(settings["type_analysis"], bin_key, res, status)
        # if status is False: prob_bins.append(bin_key)

    print(f"NUM of problematic bins = {len(prob_bins)}")
    print(prob_bins)
    ws.writeToFile(settings["ws_name"])

###############################################################################


def runParallelFits(ws_name, bins_dict, fit_settings, import_pdfs=False,  
                    savefigs=False, figpath={"good":"figs/stuff", "check":"figs/check/stuff"}):
    """
    """

    path = os.path.dirname(__file__)
    ROOT.gSystem.cd(path)

    checkImportCustomPdfs(fit_settings["bkg_model"])

    file_ws = ROOT.TFile(ws_name)
    ws = file_ws.Get("w")

    prob_bins = []
    if fit_settings["type_analysis"] == "indep":
        dict_flags = ["pass", "fail"]
        fitters = [IndepFitter(bin_key, fit_settings) for bin_key in bins_dict.keys()]
    elif fit_settings["type_analysis"] == "sim" or fit_settings["type_analysis"] == "sim_sf":
        dict_flags = ["sim"]
        fitters = [SimFitter(bin_key, fit_settings) for bin_key in bins_dict.keys()]
    else:
            print("ERROR: wrong analysis type indicated")
            sys.exit()
    
    figpath = figpath if savefigs is True else None

    pool = Pool(processes=8)
    fit_outputs = pool.starmap( 
        doSingleFit, zip(fitters, repeat(ws), repeat(dict_flags), repeat(savefigs), repeat(figpath)))

    print("FITS DONE")
    
    for fitter_out, res_fit, status_fit in fit_outputs:
        
        if import_pdfs: fitter_out.importFitObjects(ws)

        print(f"Analyzing bin {fitter_out.bin_key}")
        if status_fit is False: prob_bins.append(fitter_out.bin_key)
        printFitStatus(fit_settings["type_analysis"], fitter_out.bin_key, res_fit, status_fit, prob_bins)
    


    print(f"NUM of problematic bins = {len(prob_bins)}")
    print(prob_bins)
    ws.writeToFile(ws_name)

