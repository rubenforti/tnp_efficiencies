"""
"""

import ROOT
import os
import pickle
from utilities import import_Steve_histos, eval_efficiency, add_result
from fit_utilities import fit_with_bkg


def sim_efficiency(type_eff, bin_pt, bin_eta, results,
                   test_bkg=False, saveplots=False):
    """
    """

    h_data, h_mc, n_events, x = import_Steve_histos(
        type_eff, [bin_pt], [bin_eta])

    path = os.path.dirname(__file__)
    ROOT.gSystem.cd(path)

    bins = (bin_pt, bin_eta)

    tau1 = ROOT.RooRealVar("tau1", "tau1", -10, 0)
    expo1 = ROOT.RooExponential("expo1", "expo1", x, tau1)

    tau2 = ROOT.RooRealVar("tau2", "tau2", -10, 0)
    expo2 = ROOT.RooExponential("expo2", "expo2", x, tau2)

    pdf_mc_pass = ROOT.RooHistPdf("pdf_mc_pass", "pdf_mc_pass", x, h_mc[1])
    pdf_mc_fail = ROOT.RooHistPdf("pdf_mc_fail", "pdf_mc_fail", x, h_mc[0])

    mean = ROOT.RooRealVar("mean", "mean", 0, -2, 2)
    sigma = ROOT.RooRealVar("sigma", "sigma", 0.5, 0.001, 2)
    smearing = ROOT.RooGaussian("smearing", "smearing", x, mean, sigma)

    print(type(h_data))
    print(type(h_data[0]), type(h_data[1]))

    x.setBins(1000, "cache")
    conv_pass = ROOT.RooFFTConvPdf("conv", "conv", x, pdf_mc_pass, smearing, 3)
    conv_pass.setBufferFraction(0.1)
    conv_fail = ROOT.RooFFTConvPdf("conv", "conv", x, pdf_mc_fail, smearing, 3)
    conv_fail.setBufferFraction(0.1)

    Nsig_pass = ROOT.RooRealVar(
        "nsig_p", "#signal events pass", 250, 0, n_events[0][1])
    Nbkg_pass = ROOT.RooRealVar(
        "nbkg_p", "#background events pass", 0, n_events[0][1])
    f_pass = ROOT.RooRealVar("f_pass", "f_pass", 0, 1)

    Nsig_fail = ROOT.RooRealVar(
        "nsig_f", "#signal events fail", 0, n_events[0][0])
    Nbkg_fail = ROOT.RooRealVar(
        "nbkg_f", "#background events fail", 0, n_events[0][0])
    f_fail = ROOT.RooRealVar("f_fail", "f_fail", 0, 1)

    sum_pass = ROOT.RooAddPdf("sum_pass", "sum_pass", [
                              conv_pass, expo1], [Nsig_pass, Nbkg_pass])
    model_pass = ROOT.RooAddPdf(sum_pass)

    eff = ROOT.RooRealVar("eff", "eff", 0, 1)
    minus_eff = ROOT.RooPolyVar("minus_eff", "minus_eff", eff, [0, -1.])
    prod = ROOT.RooProduct("prod", "prod", [minus_eff, Nsig_pass])
    Nsig_fail_new = ROOT.RooAddition(
        "new_Nsig_fail", "new_Nsig_fail", [Nsig_pass, prod])
    sum_fail = ROOT.RooAddPdf("sum_fail", "sum_fail", [conv_fail, expo2], [
                              Nsig_fail_new, Nbkg_fail])
    model_fail = ROOT.RooAddPdf(sum_fail)

    sample = ROOT.RooCategory("sample", "sample")
    sample.defineType("pass")
    sample.defineType("fail")

    toy_data_pass = model_pass.generate({x}, 10000)
    toy_data_fail = model_fail.generate({x}, 10000)

    comb_dataset = ROOT.RooDataHist(
        "combData",
        "combined datasets",
        ROOT.RooArgSet(x),
        Index=sample, Import={"pass": h_data[1], "fail": h_data[0]})

    print("Data combination OK")

    simPdf = ROOT.RooSimultaneous("simPdf", "simultaneous pdf", sample)
    simPdf.addPdf(model_pass, "pass")
    simPdf.addPdf(model_fail, "fail")

    fitResult = simPdf.fitTo(comb_dataset, Save=True,
                             PrintLevel=-1, Extended=True)
    fitResult.Print()

    '''
    eff, d_eff = eval_efficiency(Npass, Nfail, sigma_Npass, sigma_Nfail)

    print(f'Measured efficiency for {type_eff} is: {eff} +- {d_eff}')

    updated_results = add_result(results, res_pass, res_fail,
                                 (eff, d_eff), bin_pt, bin_eta)

    return updated_results
    '''


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

    results = {}

    # for bin_pt in range(15):
    # for bin_eta in range(48):
    sim_efficiency(t, 1, 1, results)

    '''
    with open("simult_eff_results.pkl", "wb") as f:
        pickle.dump(results, f)
        f.close()
    '''
    print("RISULTATI SCRITTI SU PICKLE FILE")
