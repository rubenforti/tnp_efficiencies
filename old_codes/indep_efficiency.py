"""
"""

import ROOT
import os
import pickle
from utilities import import_Steve_histos, eval_efficiency, add_result
from fit_utilities import fit_with_bkg
from results_utilities import res_manager_indep


def indep_efficiency(workspace, bin_pt, bin_eta, results,
                     test_bkg=False, saveplots=False, verbose=-1):
    """
    """

    h_data, h_mc, n_events, x = import_Steve_histos(type, [bin_pt], [bin_eta])

    path = os.path.dirname(__file__)
    ROOT.gSystem.cd(path)

    bins = (bin_pt, bin_eta)

    tau = ROOT.RooRealVar("tau", "tau", -10, 0)
    expo = ROOT.RooExponential("expo", "expo", x, tau)

    res_pass = fit_with_bkg(x, type, h_data[1], h_mc[1], expo, bins,
                            n_events[0][1], saveplot=saveplots, verb=verbose)

    for par in res_pass.floatParsFinal():
        if par.GetName() == 'nsig':
            Npass = par.getVal()
            sigma_Npass = par.getError()

    res_fail = fit_with_bkg(x, type, h_data[0], h_mc[0], expo, bins,
                            n_events[0][0], saveplot=saveplots, verb=verbose)

    for par in res_fail.floatParsFinal():
        if par.GetName() == 'nsig':
            Nfail = par.getVal()
            sigma_Nfail = par.getError()

    eff, d_eff = eval_efficiency(Npass, Nfail, sigma_Npass, sigma_Nfail)

    print(f'Measured efficiency for {t} is: {eff} +- {d_eff}')
    print(res_pass.status(), res_fail.status())
    print(res_pass.covQual(), res_fail.covQual())
    print(res_pass.edm(), res_fail.edm())
    results.add_result(res_pass, res_fail, (eff, d_eff), bin_pt, bin_eta)

    # return updated_results


if __name__ == '__main__':

    custom_pdfs = ['RooCBExGaussShape',
                   'RooDoubleCBFast', 'RooCMSShape', 'my_double_CB']
    # import_pdf_library(custom_pdfs[2])

    type_eff = ("sa", "global", "ID", "iso", "trigger", "veto")
    t = type_eff[3]

    '''
    idx_cond = 1  # One for "pass", zero for fail
    id_flag = "fail" if idx_cond == 0 else "pass"
    '''

    results = res_manager_indep()

    pt_sel = [1]
    eta_sel = [24]

    for pt in pt_sel:
        for eta in eta_sel:
            indep_efficiency(t, pt, eta, results, saveplots=True, verbose=0)
    '''
    for bin_pt in range(1, 2):
        for bin_eta in range(1, 49):
            indep_efficiency(t, bin_pt, bin_eta, results)
    '''

    probs = results.problematic_bins('all')
    print(f'NUM PROBLEMI = {len(probs)}')
    # results.write("indep_eff_results.pkl")
    # print("RISULTATI SCRITTI SU PICKLE FILE")