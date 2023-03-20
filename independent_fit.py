"""
"""

import ROOT
import os
from utilities import import_Steve_histos, eval_efficiency
from fit_distrib import fit_distribution

if __name__ == '__main__':

    custom_pdfs = ['RooCBExGaussShape',
                   'RooDoubleCBFast', 'RooCMSShape', 'my_double_CB']
    # import_pdf_library(custom_pdfs[2])

    type_eff = ("sa", "global", "ID", "iso", "trigger", "veto")
    t = type_eff[3]

    idx_cond = 1  # One for "pass", zero for fail
    id_flag = "fail" if idx_cond == 0 else "pass"

    NBINS_MASS = 80

    h_data, h_mc, n_events, x = import_Steve_histos(t, [1], [1])

    path = os.path.dirname(__file__)
    ROOT.gSystem.cd(path)

    print(
        f"Num events in data and mc = {n_events[0][idx_cond]}, {n_events[1][idx_cond]}")

    tau = ROOT.RooRealVar("tau", "tau", -10, 0)
    expo = ROOT.RooExponential("expo", "expo", x, tau)

    res_pass = fit_distribution(x, t, h_data[1], h_mc[1], expo, 1, n_events[0][1], saveplot=True)
    for par in res_pass.floatParsFinal():
        if par.GetName() == 'nsig':
            Npass = par.getVal()
            sigma_Npass = par.getError()

    res_fail = fit_distribution(x, t, h_data[0], h_mc[0], expo, 0, n_events[0][0], saveplot=True)
    for par in res_fail.floatParsFinal():
        if par.GetName() == 'nsig':
            Nfail = par.getVal()
            sigma_Nfail = par.getError()

    eff, d_eff = eval_efficiency(Npass, Nfail, sigma_Npass, sigma_Nfail)

    print(f'Measured efficiency for {t} is: {eff} +- {d_eff}')
