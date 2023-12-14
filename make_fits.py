"""
"""

import ROOT
import os
import sys
from multiprocessing import Pool
from itertools import repeat
from fitter import IndepFitter, SimFitter
from utilities.dataset_utils import import_pdf_library



def printFitStatus(type_analysis, bin_key, res, status, prob_bins):

    if status is True:
        print(f"Bin {bin_key} is OK!\n\n")
    else:
        print(f"Bin {bin_key} has problems!\n")
        prob_bins.append(bin_key)
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
        "cmsshape" : "RooCMSShape_old",
        "cmsshape_w_prefitSS" : "RooCMSShape",
        "cmsshape_new" : "RooCMSShape_mod",
        "CB" : "my_double_CB"
    }
    for flag in ["pass", "fail"]:
        if bkg_models[flag] in dict_classes.keys():
            import_pdf_library(dict_classes[bkg_models[flag]])
            print(f"Imported {dict_classes[bkg_models[flag]]} from pdf_library")

###############################################################################


def doSingleFit(fitter, ws, flags, savefigs, figpath=None):
    """
    """
    fitter.manageFit(ws)

    res = {}
    for flag in flags: res.update({flag : fitter.res_obj[flag]})

    if savefigs is True and fitter.existingFit is False: fitter.saveFig(ws, figpath)
    bin_key = fitter.bin_key

    return fitter, res, fitter.bin_status

###############################################################################


def runFits(ws_name, bin_dict, fit_settings, parallelize=True, import_pdfs=False,  
            savefigs=False, figpath={"good":"figs/stuff", "check":"figs/check/stuff"}):
    """
    """
    path = os.path.dirname(__file__)
    ROOT.gSystem.cd(path)

    checkImportCustomPdfs(fit_settings["bkg_model"])

    file_ws = ROOT.TFile(ws_name)
    ws = file_ws.Get("w")

    prob_bins = []


    for bin_key in bin_dict.keys():
        
        # if bin_key != "[45.0to55.0][0.7to0.8]": continue
        # if bin_key != "[35.0to45.0][0.7to0.8]": continue

        if fit_settings["type_analysis"] == "indep":
            dict_flags = ["pass", "fail"]
            fitter = IndepFitter(bin_key, fit_settings)            
        elif fit_settings["type_analysis"] == "sim" or fit_settings["type_analysis"] == "sim_sf":
            dict_flags = ["sim"]
            fitter = SimFitter(bin_key, fit_settings)
        else:
            print("ERROR: wrong analysis type indicated")
            sys.exit()

        fitter, res, status = doSingleFit(fitter, ws, dict_flags, savefigs, figpath=figpath)

        if import_pdfs: fitter.importFitObjects(ws)

        printFitStatus(fit_settings["type_analysis"], bin_key, res, status, prob_bins)
        # if status is False: prob_bins.append(bin_key)

    print(f"NUM of problematic bins = {len(prob_bins)}")
    print(prob_bins)
    ws.writeToFile(ws_name)

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
