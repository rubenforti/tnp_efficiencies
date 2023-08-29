"""
"""

import ROOT
import os
import sys
import time
import argparse
from utilities.base_library import bin_dictionary, lumi_factors
from independent_fitter import EfficiencyFitter as indep_eff
from simultaneous_fitter import EfficiencyFitter as sim_eff
from utilities.fit_utils import check_existing_fit


def printFitStatus(type_analysis, bin_key, status, res, prob_bins):

    if status is True:
        print(f"\n\nBin {bin_key} is OK!\n")
    else:
        print(f"\n\nBin {bin_key} has problems!\n")
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
            res.Print()
            res.correlationMatrix().Print()
            print("****")
            print(res.status())
            print(res.covQual())
            print(res.edm())
            print('\n')



def runFits(ws_name, type_eff, bins_dict, fit_settings, fit_verb=-1, import_pdfs=False, refit_numbkg=True, 
              savefigs=False, figpath={"good":"figs/stuff", "check":"figs/check/stuff"}):
    """
    """
    path = os.path.dirname(__file__)
    ROOT.gSystem.cd(path)

    file_ws = ROOT.TFile(ws_name)
    ws = file_ws.Get("w")

    prob_bins = []

    for bin_key in bins_dict.keys():
        
        if bin_key!="[24.0to26.0][-2.4to-2.3]":
            sys.exit()

        if check_existing_fit(fit_settings["type_analysis"], ws, bin_key) == 0:

            if fit_settings["type_analysis"] == "indep":
                
                print("AAA")
                fitter = indep_eff(bin_key, fit_settings)
                res, status_dict = {}, {}   

                for flag in ["pass", "fail"]:
                    result, status = fitter.doFit(flag, ws)
                    
                    pars_fitted = result.floatParsFinal()
                    nsig_fitted = pars_fitted.find(f"nsig_{flag}_{bin_key}")
                    nbkg_fitted = pars_fitted.find(f"nbkg_{flag}_{bin_key}")
                
                    low_nbkg = (nbkg_fitted.getVal() < 0.005*nsig_fitted.getVal())

                    if (status is False) and refit_numbkg and low_nbkg:
                        result, status = fitter.refit_noBkg(flag, ws)

                    res.update({flag : result})
                    status_dict.update({flag : status})

                status = bool(status_dict["pass"]*status_dict["fail"])

                if status and import_pdfs:
                    ws.Import(fitter.getFinalPdf("pass")), ws.Import(fitter.getFinalPdf("fail"))
                    ws.Import(res["pass"]), ws.Import(res["fail"])
                
                if savefigs:
                    fitter.saveFig(ws, res, status, figpath)


            elif fit_settings["type_analysis"] == "sim" or fit_settings["type_analysis"] == "sim_sf":

                fitter = sim_eff(bin_key, fit_settings)

                res, status = fitter.doFit(ws)
                
                pars_fitted = res.floatParsFinal()

                refit_flags = []
                if (status is False) and refit_numbkg:
                    for flag in ["pass", "fail"]:
                        nsig_fitted = pars_fitted.find(f"nsig_{flag}_{bin_key}")
                        nbkg_fitted = pars_fitted.find(f"nbkg_{flag}_{bin_key}")
            
                        if (nbkg_fitted.getVal() < 0.005*nsig_fitted.getVal()):
                            refit_flags.append(flag)

                    if len(refit_flags) != 0:
                        res, status = fitter.refit_noBkg(refit_flags, ws)

                if status and import_pdfs:
                    ws.Import(fitter.getFinalPdf("sim")), ws.Import(res)

                if savefigs:
                    fitter.saveFig(ws, res, status, figpath)
            
            else:
                print("ERROR: wrong analysis type indicated")
                sys.exit()
        
        printFitStatus(fit_settings["type_analysis"], bin_key, status, res, prob_bins)
                

    print(f"NUM of problematic bins = {len(prob_bins)}")
    print(prob_bins)

    return ws


if __name__ == '__main__':



    ROOT.gROOT.SetBatch(True)
    ROOT.PyConfig.IgnoreCommandLineOptions = True


    t0 = time.time()
    custom_pdfs = ['RooCBExGaussShape',
                   'RooDoubleCBFast', 'RooCMSShape', 'my_double_CB']
    # import_pdf_library(custom_pdfs[2])

    # -------------------------------------------------------------------------
    # GENERAL SETTINGS
    # ----------------

    type_eff = "iso"

    type_analysis = "indep"
        
    bins = bin_dictionary("pt", "eta")


    # -------------------------------------------------------------------------
    # Dataset generation
    # ------------------


    '''
    filename_data = f"/scratchnvme/wmass/Steve_root_files/Standard_SF_files/tnp_{type_eff}_data_vertexWeights1_oscharge1.root"
    filename_mc = f"/scratchnvme/wmass/Steve_root_files/Standard_SF_files/tnp_{type_eff}_mc_vertexWeights1_oscharge1.root"
    dirname_bkg = "/scratchnvme/rajarshi/Bkg_TNP_3D_Histograms/OS"
    '''

    
    filename_data = f"root_files/datasets/tnp_{type_eff}_data_vertexWeights1_oscharge1.root"
    filename_mc = f"root_files/datasets/tnp_{type_eff}_mc_vertexWeights1_oscharge1.root"
    dirname_bkg = "root_files/datasets"
    
    bkg_types = ["WW", "WZ", "ZZ", "TTSemileptonic", "Ztautau"]

    bkg_filenames = {}
    [bkg_filenames.update({cat : 
        f"{dirname_bkg}/tnp_{type_eff}_{cat}_vertexWeights1_oscharge1.root"}) for cat in bkg_types]
    

    lumi_scales = lumi_factors(type_eff, bkg_types)

    lumi_scale_signal = lumi_scales.pop("Zmumu")
    

    import_dictionary = {
        "data" : filename_data,
        "mc" : {
            "filename": filename_mc,
            "lumi_scale" : lumi_scale_signal
        },
        # "bkg" : {
        #    "filenames" : bkg_filenames,
        #   "lumi_scales" : lumi_scales
        # } 
    }

    # workspace_name = f"root_files/ws/ws_bkg_studies.root"
    # workspace_name = "root_files/ws_{type_eff}_indep_mcbkg_mergedbins.root"

    binning_mass = "mass_60_120"

    workspace_name = f"root_files/iso_workspace.root"
    
    # ws = ws_init(import_dictionary, type_analysis, "pt", "eta", binning_mass)
    # ws.writeToFile(workspace_name)

    # file = ROOT.TFile(workspace_name, "READ")
    # ws = file.Get("w")

    '''
    import_dict_bkg =  {
        "bkg" : {
            "filenames" : bkg_filenames,
            "lumi_scales" : lumi_scales
        } 
    }
    # ws_bkg = ws_init(import_dict_bkg, type_analysis, "pt_12bins", "eta_16bins", binning_mass, 
    #                  import_existing_ws=True, existing_ws_filename=workspace_name, altBinning_bkg=True)
    # ws_bkg.writeToFile(workspace_name)
    '''


    # ------------------------------------------------------------------------
   
    # results_name = f"results/results_{type_eff}_{type_analysis}.pkl"
   
    figpath = {"good": "figs/", "check": "figs/"} 

    ws = make_fits(workspace_name, type_eff, type_analysis, bins, fit_settings_iso, 
                  import_pdfs=False, savefigs=True, figpath=figpath)

    # ws.writeToFile(workspace_name)
 

    #print("RISULTATI SCRITTI SU PICKLE FILE")

    t1 = time.time()

    print(f"TEMPO = {t1-t0}")
