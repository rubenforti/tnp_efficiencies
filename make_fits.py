"""
"""

import ROOT
import os
import sys
import time
import argparse
from utilities.base_library import bin_dictionary, binning, bkg_lumi_scales
from indep_efficiency import independent_efficiency
from sim_efficiency import simultaneous_efficiency
from utilities.fit_utils import check_existing_fit, fit_quality
from utilities.dataset_utils import ws_init



def make_fits(ws_name, type_eff, type_analysis, bins_dictionary,
              bkg_pdf='expo', fit_verbosity=-1, savefigs=False):
    """
    """
    path = os.path.dirname(__file__)
    ROOT.gSystem.cd(path)


    file_ws = ROOT.TFile(ws_name)
    ws = file_ws.Get("w")

    Nproblems = 0
    problematic_bins = {}

    for bin_key in bins_dictionary:

        global_idx, bin_pt, bin_eta = bins_dictionary[bin_key]

        existingRes = check_existing_fit(type_analysis, ws, bin_key)

        if type_analysis == 'indep':

            print(type(existingRes[0]), type(existingRes[1]))

            if (existingRes[0] == 0 and existingRes[1] == 0):
                res_pass, res_fail = independent_efficiency(ws, bin_key, bkg_pdf, 
                                                            test_bkg=False, refit_numbkg=True, 
                                                            verb=fit_verbosity, figs=savefigs)
            else:
                res_pass, res_fail = existingRes
            
            status = bool(
                fit_quality(res_pass, old_checks=True)*fit_quality(res_fail, old_checks=True))

            if status is False or 1!=0:
                print(f"\nBin {bin_key} ({bin_pt}|{bin_eta}) has problems!\n")
                Nproblems += 1
                print("****")
                res_pass.Print()
                res_pass.correlationMatrix().Print()
                '''
                pars_pass = res_pass.floatParsFinal()
                nsig_pass = pars_pass.find(f"nsig_pass_{bin_key}")
                print((nsig_pass.getVal()**0.5, nsig_pass.getError()))
                '''
                print("****")
                res_fail.Print()
                res_fail.correlationMatrix().Print()
                '''

                pars_fail = res_fail.floatParsFinal()
                nsig_fail = pars_fail.find(f"nsig_fail_{bin_key}")
                print((nsig_fail.getVal()**0.5, nsig_fail.getError()))
                '''
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
    
        
    print(f"NUM of problematic bins = {Nproblems}\n")

    print("List of non fitted bins:")
    for probl_bin_key in problematic_bins:
        print(problematic_bins[probl_bin_key][1], problematic_bins[probl_bin_key][2])

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

        
    bins = bin_dictionary("pt", "eta_8bins")


    # -------------------------------------------------------------------------
    # Dataset generation
    # ------------------



    filename_data = "/scratchnvme/wmass/Steve_root_files/Standard_SF_files/tnp_iso_data_vertexWeights1_oscharge1.root"
    filename_mc = "/scratchnvme/wmass/Steve_root_files/Standard_SF_files/tnp_iso_mc_vertexWeights1_oscharge1.root"
    dirname_bkg = "/scratchnvme/rajarshi/Bkg_TNP_3D_Histograms/OS"
    bkg_filepaths = {}

    bkg_categories= ["WW", "WZ", "ZZ", "TTSemileptonic", "Ztautau"]
    
    [bkg_filepaths.update({cat : 
        f"{dirname_bkg}/tnp_{type_eff}_{cat}_vertexWeights1_oscharge1.root"}) for cat in bkg_categories]
    
    print(bkg_filepaths)
    

    lumi_scales = bkg_lumi_scales(type_eff, bkg_categories)

    print(lumi_scales)

    

    import_dictionary = {
        "data" : filename_data,
        "mc" : filename_mc,
        "bkg" : {
            "filepaths" : bkg_filepaths,
            "lumi_scales" : lumi_scales
        }
    }

    workspace_name = f"root_files/ws/ws_data_mc_bkg.root"
    ws = ws_init(import_dictionary, type_analysis, bins, binning("mass_60_120"))
    ws.writeToFile(workspace_name)

    
    # ------------------------------------------------------------------------


   
    results_name = f"results/results_{type_eff}_{type_analysis}.pkl"

    ws = make_fits(workspace_name, type_eff, type_analysis, bins, savefigs=True)

    ws.writeToFile(workspace_name)
 

    #print("RISULTATI SCRITTI SU PICKLE FILE")

    t1 = time.time()

    print(f"TEMPO = {t1-t0}")
