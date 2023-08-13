"""
"""

import ROOT
import os
import sys
import time
import argparse
from utilities.base_library import bin_dictionary, binning, lumi_factors
from indep_efficiency import independent_efficiency
# from indep_eff_mcbkg import independent_efficiency
# from sim_efficiency import simultaneous_efficiency
from utilities.fit_utils import check_existing_fit, fit_quality
from utilities.dataset_utils import ws_init, extend_merged_datasets



settings_dict = {
    "strategy" : "indep",
    "id_name" : "AAAAAAAAAAAAAAAAAA",
    "fit_range" : [60.0, 120.0],
    "bkg_categories" : ["WW", "WZ", "ZZ", "TTSemileptonic", "Ztautau"],
    "fit_on_pseudodata" : True,
    "run" : 1,
    "pass" : {
        "fit_strategy" : 2,
        "Nbins" : 5000,
        "bufFraction" : 0.1, 
        "bkg_shape": "expo",
        "pars" : {
            "mu" : [0.0, -5.0, 5.0],
            "sigma" : [0.5, 0.01, 5.0],
            "tau" : [0.0, -5.0, 5.0],
        },
        "norm" : {
            "nsig" : ["1n", 0.5, "3n"],
            "nbkg" : ["0.05n", 0.5, "1.5n"],
        },
    },
    "fail" : {
        "fit_strategy" : 2,
        "Nbins" : 8000,
        "bufFraction" : 0.5,
        "bkg_shape" : "expo",
        "pars" : {
            "mu" : [0.2, -5.0, 5.0],
            "sigma" : [1.5, 0.2, 5.0],
            "tau" : [-0.5, -2.0, 2.0],
        },
        "norm" : {
            "nsig" : ["1n", 0.5, "3n"],
            "nbkg" : ["0.1n", 0.5, "1.5n"],
        }
    }
}


def make_fits(ws_name, type_eff, type_analysis, bins_dict, fit_settings, fit_verb=-1, import_pdfs=False, 
              savefigs=False, figpath={"good":"figs/stuff", "check":"figs/check/stuff"}):
    """
    """
    path = os.path.dirname(__file__)
    ROOT.gSystem.cd(path)

    file_ws = ROOT.TFile(ws_name)
    ws = file_ws.Get("w")

    prob_bins = []

    # key_prova = "[24.0to26.0][-2.4to-2.3]"
    # bins_dictionary = {key_prova : bins_dictionary[key_prova]}

    if "ws_bkg_merged" in fit_settings.keys():
        file_bkg_merged = ROOT.TFile(fit_settings["ws_bkg_merged"], "READ")
        ws_bkg_merged = file_bkg_merged.Get("w")
        fit_settings.update({"ws_bkg_merged" : ws_bkg_merged})

    for bin_key in bins_dict:

        global_idx, bin_pt, bin_eta = bins_dict[bin_key]

        existingRes = check_existing_fit(type_analysis, ws, bin_key)

        if type_analysis == 'indep':

            if existingRes == 0:

                res_pass, res_fail, status = independent_efficiency(ws, bin_key, fit_settings, 
                                                                    refit_numbkg=True, verb=fit_verb, 
                                                                    import_pdfs=import_pdfs, figs=savefigs,
                                                                    figpath=figpath)
            else:
                res_pass, res_fail = existingRes
                status = True
    

            if status is False:
                print(f"\nBin {bin_key} ({bin_pt}|{bin_eta}) has problems!\n")
                prob_bins.append(bin_key)
                print("****")
                res_pass.Print()
                res_pass.correlationMatrix().Print()
                print("****")
                res_fail.Print()
                res_fail.correlationMatrix().Print()
                print("****")
                print(res_pass.status(), res_fail.status())
                print(res_pass.covQual(), res_fail.covQual())
                print(res_pass.edm(), res_fail.edm())
                print('\n')
            else:
                pass
                #res_object.add_result(ws, bin_pt, bin_eta)
            
        elif type_analysis == 'sim':

            if existingRes == 0:
                results = simultaneous_efficiency(ws, bin_key, bkg_pdf, 
                                                  refit_numbkg=True, test_bkg=False, same_smearing=False, 
                                                  enable_mcfit=False, verb=fit_verbosity, figs=savefigs)
            else: 
                results = existingRes

            status = fit_quality(results, old_checks=True)

            if status is False or 1!=0:
                print(f"\nBin {bin_pt},{bin_eta} has problems!\n")
                Nproblems += 1
                problematic_bins.update({bin_key: bin_dictionary[bin_key]})
                results.Print()
                #results.correlationMatrix().Print()
                print("****")
                print(results.status())
                print(results.covQual())
                print(results.edm())
                print(' ')
            else:
                pass
                #res_object.add_result(ws, bin_pt, bin_eta)
                
        else:
            print("Wrong fit strategy! Closing the program")
            sys.exit()
    
        
    print(f"NUM of problematic bins = {len(prob_bins)}")
    print(prob_bins)

    '''
    print("List of non fitted bins:")
    for probl_bin_key in problematic_bins:
        print(problematic_bins[probl_bin_key][1], problematic_bins[probl_bin_key][2])
    '''

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

    types_efficiencies = ("sa", "global", "ID", "iso", "trigger", "veto")
    type_eff = types_efficiencies[3]

    type_analysis = 'indep'

    bkg_pdf = 'expo'

        
    bins = bin_dictionary("pt", "eta")


    # -------------------------------------------------------------------------
    # Dataset generation
    # ------------------


    '''
    filename_data = "/scratchnvme/wmass/Steve_root_files/Standard_SF_files/tnp_iso_data_vertexWeights1_oscharge1.root"
    filename_mc = "/scratchnvme/wmass/Steve_root_files/Standard_SF_files/tnp_iso_mc_vertexWeights1_oscharge1.root"
    dirname_bkg = "/scratchnvme/rajarshi/Bkg_TNP_3D_Histograms/OS"
    '''
    filename_data = "root_files/datasets/tnp_iso_data_vertexWeights1_oscharge1.root"
    filename_mc = "root_files/datasets/tnp_iso_mc_vertexWeights1_oscharge1.root"
    dirname_bkg = "root_files/datasets"
    
    
    bkg_filenames = {}
    [bkg_filenames.update({cat : 
        f"{dirname_bkg}/tnp_{type_eff}_{cat}_vertexWeights1_oscharge1.root"}) for cat in settings_dict["bkg_categories"]]
    
    print(bkg_filenames)
    

    lumi_scales = lumi_factors(type_eff, settings_dict["bkg_categories"])

    print(lumi_scales)


    lumi_scale_signal = lumi_scales.pop("Zmumu")
    

    import_dictionary = {
        "data" : filename_data,
        "mc" : {
            "filename": filename_mc,
            "lumi_scale" : lumi_scale_signal
        },
        "bkg" : {
           "filenames" : bkg_filenames,
           "lumi_scales" : lumi_scales
        } 
    }

    # workspace_name = f"root_files/ws/ws_bkg_studies.root"
    #  workspace_name = "root_files/ws_iso_indep_mcbkg_mergedbins.root"

    binning_mass = "mass_60_120"

    workspace_name = "root_files/ws_iso_pseudodata.root"
    
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
   
    figpath = {"good": "figs/pseudodata", "check": "figs/check_fits/pseudodata"} 

    ws = make_fits(workspace_name, type_eff, type_analysis, bins, settings_dict, 
                  import_pdfs=True, savefigs=True, figpath=figpath)

    ws.writeToFile(workspace_name)
 

    #print("RISULTATI SCRITTI SU PICKLE FILE")

    t1 = time.time()

    print(f"TEMPO = {t1-t0}")
