"""
"""

import ROOT
import os
import sys
from multiprocessing import Pool
from itertools import repeat
from fitter import IndepFitter, SimFitter


def printFitStatus(type_analysis, bin_key, status, res, prob_bins):

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
        elif type_analysis == "sim":
            print("****")
            res["sim"].Print()
            res["sim"].correlationMatrix().Print()
            print("****")
            print(res["sim"].status())
            print(res["sim"].covQual())
            print(res["sim"].edm())
            print('\n\n')


def runFits(ws_name, bins_dict, fit_settings, parallelize=True, import_pdfs=False,  
            savefigs=False, figpath={"good":"figs/stuff", "check":"figs/check/stuff"}):
    """
    """
    path = os.path.dirname(__file__)
    ROOT.gSystem.cd(path)

    file_ws = ROOT.TFile(ws_name)
    ws = file_ws.Get("w")

    prob_bins = []

    for bin_key in bins_dict.keys():
        
        # if bin_key != "[60.0to65.0][-1.3to-1.2]": continue

        if fit_settings["type_analysis"] == "indep":
            fitter = IndepFitter(bin_key, fit_settings)            
            fitter.manageFit(ws)
            status = bool(fitter.results["pass"]["status"]*fitter.results["fail"]["status"])
            res = {"pass" : fitter.results["pass"]["res_obj"], "fail" : fitter.results["fail"]["res_obj"]}

    
        elif fit_settings["type_analysis"] == "sim" or fit_settings["type_analysis"] == "sim_sf":
            fitter = SimFitter(bin_key, fit_settings)
            fitter.manageFit(ws)
            res = {"sim" : fitter.results["sim"]["res_obj"]}
            status = bool(fitter.results["sim"]["status"])
        
        else:
            print("ERROR: wrong analysis type indicated")
            sys.exit()

        if import_pdfs: fitter.importFitObjects(ws)
        if savefigs: fitter.saveFig(ws, figpath)

        printFitStatus(fit_settings["type_analysis"], bin_key, status, res, prob_bins)


    print(f"NUM of problematic bins = {len(prob_bins)}")
    print(prob_bins)
    ws.writeToFile(ws_name)




def parallelize_fits(fitter, ws, import_pdfs, savefigs, figpath=None):
    """
    """
    fitter.manageFit(ws)

    if fitter.settings["type_analysis"] == "indep":
        status = bool(fitter.results["pass"]["status"]*fitter.results["fail"]["status"])
        res = {"pass" : fitter.results["pass"]["res_obj"], "fail" : fitter.results["fail"]["res_obj"]}
    elif fitter.settings["type_analysis"] == "sim" or fitter.settings["type_analysis"] == "sim_sf":
        res = {"sim" : fitter.results["sim"]["res_obj"]}
        status = bool(fitter.results["sim"]["status"])

    if import_pdfs is True: fitter.importFitObjects(ws)

    if savefigs is True: fitter.saveFig(ws, figpath)
    bin_key = fitter.bin_key

    return (fitter, bin_key, status, res)


def runParallelFits(ws_name, bins_dict, fit_settings, import_pdfs=False,  
                    savefigs=False, figpath={"good":"figs/stuff", "check":"figs/check/stuff"}):
    """
    """

    # ATTENTION: TO NOW IT IS ONLY AVAILABLE FOR INDEPENDENT FITS !!!

    path = os.path.dirname(__file__)
    ROOT.gSystem.cd(path)

    file_ws = ROOT.TFile(ws_name)
    ws = file_ws.Get("w")

    prob_bins = []

    if fit_settings["type_analysis"] == "indep":
        fitters = [IndepFitter(bin_key, fit_settings) for bin_key in bins_dict.keys()]
    elif fit_settings["type_analysis"] == "sim" or fit_settings["type_analysis"] == "sim_sf":
        fitters = [SimFitter(bin_key, fit_settings) for bin_key in bins_dict.keys()]
    
    figpath = figpath if savefigs is True else None
    pool = Pool(processes=16)
    fit_results = pool.starmap(parallelize_fits, 
                               zip(fitters, repeat(ws), repeat(import_pdfs), repeat(savefigs), repeat(figpath)))


    print("FITS DONE")

    
    for fitres in fit_results:
        fitter, bin_key, bin_status, res = fitres
        if import_pdfs: fitter.importFitObjects(ws)

        print(f"Analyzing bin {bin_key}")
        if bin_status is False:
            prob_bins.append(fitres["bin_key"])
        printFitStatus(fit_settings["type_analysis"], fitres[1], bin_status, res, prob_bins)

    

    print(f"NUM of problematic bins = {len(prob_bins)}")
    print(prob_bins)
    ws.writeToFile(ws_name)
