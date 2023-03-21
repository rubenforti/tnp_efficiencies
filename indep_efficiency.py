"""
"""

import ROOT
import os
import pickle
from utilities import import_Steve_histos, eval_efficiency, add_result
from fit_distribution import fit_without_bkg, fit_with_bkg


if __name__ == '__main__':

    custom_pdfs = ['RooCBExGaussShape',
                   'RooDoubleCBFast', 'RooCMSShape', 'my_double_CB']
    # import_pdf_library(custom_pdfs[2])

    type_eff = ("sa", "global", "ID", "iso", "trigger", "veto")
    t = type_eff[3]

    idx_cond = 1  # One for "pass", zero for fail
    id_flag = "fail" if idx_cond == 0 else "pass"

    NBINS_MASS = 80

    BIN_PT = 1
    BIN_ETA = 1

    bins = (BIN_PT, BIN_ETA)

    h_data, h_mc, n_events, x = import_Steve_histos(t, [BIN_PT], [BIN_ETA])

    path = os.path.dirname(__file__)
    ROOT.gSystem.cd(path)

    print(
        f"Num events in data and mc = {n_events[0][idx_cond]}, {n_events[1][idx_cond]}")

    tau = ROOT.RooRealVar("tau", "tau", -10, 0)
    expo = ROOT.RooExponential("expo", "expo", x, tau)

    res_pass = fit_with_bkg(
        x, t, h_data[1], h_mc[1], expo, bins, n_events[0][1], saveplot=True)
    for par in res_pass.floatParsFinal():
        if par.GetName() == 'nsig':
            Npass = par.getVal()
            sigma_Npass = par.getError()

    res_fail = fit_with_bkg(
        x, t, h_data[0], h_mc[0], expo, bins, n_events[0][0], saveplot=True)
    for par in res_fail.floatParsFinal():
        if par.GetName() == 'nsig':
            Nfail = par.getVal()
            sigma_Nfail = par.getError()

    eff, d_eff = eval_efficiency(Npass, Nfail, sigma_Npass, sigma_Nfail)

    print(f'Measured efficiency for {t} is: {eff} +- {d_eff}') 

    
    results = {}
    
    results = add_result(results, res_pass, res_fail, (eff, d_eff), BIN_PT, BIN_ETA)


    with open("indep_eff_results.pkl", "wb") as f:
        pickle.dump(results, f)
        f.close()
    print("RISULTATI SCRITTI SU PICKLE FILE")
    

    with open("indep_eff_results.pkl", "rb") as f:
        deserialized_dict = pickle.load(f)
        print(deserialized_dict)
    
