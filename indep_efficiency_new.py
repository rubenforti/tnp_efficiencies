"""
"""

import ROOT
import os
from fit_utilities_new import fit_with_bkg
from results_utilities import res_manager_indep
from workspace_config import ws_init_std_pdf
from array import array


def indep_efficiency(type_eff, bins_pt, bins_eta, results,
                     test_bkg=False, verbose=-1):
    """
    """

    path = os.path.dirname(__file__)
    ROOT.gSystem.cd(path)

    for bin_pt in bins_pt:
        for bin_eta in bins_eta:
            model_pass, res_pass = fit_with_bkg(
                type_eff, 'pass', (bin_pt, bin_eta),
                test_bkg=test_bkg, verb=verbose)
            model_fail, res_fail = fit_with_bkg(
                type_eff, 'fail', (bin_pt, bin_eta),
                test_bkg=test_bkg, verb=verbose)

            results.add_result(res_pass, res_fail, bin_pt, bin_eta)
            status = results.check_fit_status(
                bin_pt, bin_eta, conditions='all')

            print(" ")
            if status != 0:
                print(f"Bin {bin_pt},{bin_eta} has {status} problems!")
            if verbose != 0 or status != 0:
                print(res_pass.status(), res_fail.status())
                print(res_pass.covQual(), res_fail.covQual())
                print(res_pass.edm(), res_fail.edm())
                print(' ')


if __name__ == '__main__':

    custom_pdfs = ['RooCBExGaussShape',
                   'RooDoubleCBFast', 'RooCMSShape', 'my_double_CB']
    # import_pdf_library(custom_pdfs[2])

    type_eff = ("sa", "global", "ID", "iso", "trigger", "veto")
    t = type_eff[3]

    results = res_manager_indep()

    bins_pt = array('I', [num for num in range(1, 2)])
    bins_eta = array('I', [num for num in range(1, 2)])

    indep_efficiency(t, bins_pt, bins_eta, results, verbose=-1)

    ''''
    file_ws = ROOT.TFile(f"root_files/{t}_workspace.root")
    ws = file_ws.Get("w")
    # ws_init_std_pdf(ws)
    ws.Print()

    probs = results.get_problematic_bins()
    print(f'NUM PROBLEMI = {len(probs)}')
    '''

    '''
    results.write("indep_eff_results.pkl")
    print("RISULTATI SCRITTI SU PICKLE FILE")
    '''
