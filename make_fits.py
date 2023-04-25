"""
"""

import ROOT
import os
import sys
import time
from results_utils import results_manager
from indep_efficiency import fit_on_bin
from sim_efficiency import simultaneous_efficiency
from utilities import import_pdf_library, fit_quality, pearson_chi2_eval, llr_test_bkg


def make_fits(type_eff, type_estimate, bins, bin_combinations=True,
              bkg_pdf='expo', test_bkg=False,
              fit_verbosity=-1, savefigs=False):
    """
    """
    path = os.path.dirname(__file__)
    ROOT.gSystem.cd(path)

    ROOT.Math.MinimizerOptions.SetDefaultTolerance(2e-2)
    ROOT.Math.MinimizerOptions.SetDefaultMaxIterations(100000)
    ROOT.Math.MinimizerOptions.SetDefaultErrorDef(0.5)
    ROOT.Math.MinimizerOptions.SetDefaultStrategy(2)

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

    file_ws = ROOT.TFile(
        f"root_files/ws/{type_eff}_workspace_{type_estimate}.root", "READ")
    ws = file_ws.Get("w")

    bkg_pdf = 'mixed'

    Nproblems = 0

    for idx in range(len(bins_list[0])):
        bin_pt, bin_eta = bins_list[0][idx], bins_list[1][idx]

        if type_estimate == 'indep':
            res_pass = fit_on_bin(type_eff, ws, 'pass', (bin_pt, bin_eta),
                                  bkg_pdf, test_bkg=test_bkg,
                                  verb=fit_verbosity, figs=savefigs)
            res_fail = fit_on_bin(type_eff, ws, 'fail', (bin_pt, bin_eta),
                                  bkg_pdf, test_bkg=test_bkg,
                                  verb=fit_verbosity, figs=savefigs)
            # results.add_result(res_pass, res_fail, bin_pt, bin_eta, check=True)
            status = bool(fit_quality(res_pass)*fit_quality(res_fail))

        elif type_estimate == 'sim':
            results = simultaneous_efficiency(
                type_eff, bin, bkg_pdf, same_smearing=False, test_bkg=False,
                enable_mcfit=False, verb=fit_verbosity, figs=savefigs)
            status = fit_quality(results)

        if status is False:
            print(f"\nBin {bin_pt},{bin_eta} has problems!\n")
            Nproblems += 1

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
    ws.writeToFile(f"root_files/{type_eff}_workspace.root")


if __name__ == '__main__':

    t0 = time.time()
    custom_pdfs = ['RooCBExGaussShape',
                   'RooDoubleCBFast', 'RooCMSShape', 'my_double_CB']
    import_pdf_library(custom_pdfs[2])

    # -------------------------------------------------------------------------
    # GENERAL SETTINGS
    # ----------------

    types_efficiencies = ("sa", "global", "ID", "iso", "trigger", "veto")
    type_eff = types_efficiencies[3]

    type_estimate = 'simult_simple'

    bkg_pdf = 'expo'

    bins_pt = [num for num in range(1, 2)]
    bins_eta = [num for num in range(1, 2)]

    # -------------------------------------------------------------------------

    '''
    bin_keys = results.get_problematic_bins()

    bins_pt, bins_eta = [], []

    for key in bin_keys:
        bins_pt.append(key.split(',')[0])
        bins_eta.append(key.split(',')[1])

    print(bins_pt)
    print(bins_eta)
    '''

    bins = (bins_pt, bins_eta)
    make_fits(type_eff, type_estimate, bins, bkg_pdf,
              direct_mode=True, verbose=0)

    '''
    file_ws = ROOT.TFile(f"root_files/{t}_workspace.root")
    ws = file_ws.Get("w")
    ws.Print()
    '''
    results = results_manager()
    results.open(f"results/{type_eff}_results_{type_estimate}.pkl")
    results.view_efficiencies()

    #print("RISULTATI SCRITTI SU PICKLE FILE")

    t1 = time.time()

    print(f"TEMPO = {t1-t0}")
