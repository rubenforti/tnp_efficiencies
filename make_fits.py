"""
"""

import ROOT
import os
import sys
import time
from results_utils import results_manager
from indep_efficiency import indep_eff_fits
from indep_efficiency import check_existing_fit as fit_exist_indep
from sim_efficiency import simultaneous_efficiency
from utilities import import_pdf_library, fit_quality
from workspace_config import ws_init


def make_fits(ws_name, type_eff, type_estimate, bins, bin_combinations, bkg_pdf='expo', test_bkg=False,
              fit_verbosity=-1, savefigs=False):
    """
    """
    path = os.path.dirname(__file__)
    ROOT.gSystem.cd(path)

    if bin_combinations is True:
        bins_pt = []
        bins_eta = []
        for i in range(len(bins[0])):
            for j in range(len(bins[1])):
                bins_pt.append(bins[0][i])
                bins_eta.append(bins[1][j])
        bins_list = [bins_pt, bins_eta]
    else:
        bins_list = bins
    
    print(bins_list)

    file_ws = ROOT.TFile(ws_name)
    ws = file_ws.Get("w")


    Nproblems = 0
    bins_with_problems = []

    for idx in range(len(bins_list[0])):
        bin_pt, bin_eta = bins_list[0][idx], bins_list[1][idx]

        if type_estimate == 'indep':
            
            isFitted = fit_exist_indep(ws, (bin_pt, bin_eta))
            print(type(isFitted[0]), type(isFitted[1]))

            if (isFitted[0] == 0 and isFitted[1] == 0):
                res_pass, res_fail = indep_eff_fits(type_eff, "indep", ws, (bin_pt, bin_eta), 
                                                    bkg_pdf, test_bkg=False, refit_numbkg=True, 
                                                    verb=fit_verbosity, figs=savefigs)
            else:
                res_pass, res_fail = isFitted
            # results.add_result(res_pass, res_fail, bin_pt, bin_eta, check=True)
            status = bool(
                fit_quality(res_pass, old_checks=True)*fit_quality(res_fail, old_checks=True))

        elif type_estimate == 'sim':
            results = simultaneous_efficiency(type_eff, "indep", ws, (bin_pt, bin_eta), bkg_pdf, 
                                              test_bkg=False, same_smearing=False, enable_mcfit=False, 
                                              verb=fit_verbosity, figs=savefigs)
            status = fit_quality(results, old_checks=True)
        else:
            print("Wrong fit strategy! Closing the program")
            sys.exit()

        if status is False:
            print(f"\nBin {bin_pt},{bin_eta} has problems!\n")
            Nproblems += 1
            bins_with_problems.append(f"{bin_pt},{bin_eta}")
            res_pass.Print()
            res_pass.correlationMatrix().Print()
            print("****")
            res_fail.Print()
            res_fail.correlationMatrix().Print()
            print("****")
            print(res_pass.status(), res_fail.status())
            print(res_pass.covQual(), res_fail.covQual())
            print(res_pass.edm(), res_fail.edm())
            print(' ')

    print(f"NUM of problematic bins = {Nproblems}")
    print(bins_with_problems)
    ws.writeToFile(workspace_name)


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

    '''
    bins_pt = [num for num in range(1, 16)]
    bins_eta = [num for num in range(1, 49)]
    bin_combinations = True
    '''
    # ------------------------------------------------------------------------

    
    bin_keys = ['6,7']

    bin_combinations = False

    bins_pt, bins_eta = [], []

    for key in bin_keys:
        bins_pt.append(key.split(',')[0])
        bins_eta.append(key.split(',')[1])

    print(bins_pt)
    print(bins_eta)
    


    workspace_name = f"root_files/ws/ws_{type_eff}_{type_analysis}.root"

    bins = (bins_pt, bins_eta)
    make_fits(workspace_name, type_eff, type_analysis, bins, bin_combinations, 
              bkg_pdf, fit_verbosity=-1, savefigs=True)

 

    #print("RISULTATI SCRITTI SU PICKLE FILE")

    t1 = time.time()

    print(f"TEMPO = {t1-t0}")
