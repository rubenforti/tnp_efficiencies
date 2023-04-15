"""
"""

import ROOT
import os
import pickle
from utilities import import_Steve_histos
from results_utilities import res_manager_sim


def fit_on_bin(type_eff, workspace, cond, bin, test_bkg=False, verb=-1):
    """
    """


def sim_efficiency(type_eff, bin_pt, bin_eta, results,
                   same_smearing=True,
                   test_bkg=False,
                   enable_mcfit=False,
                   saveplots=False):
    """
    """

    h_data, h_mc, n_events, x = import_Steve_histos(
        type_eff, [bin_pt], [bin_eta])

    path = os.path.dirname(__file__)
    ROOT.gSystem.cd(path)

    # Backgrounds
    # -----------
    tau_p = ROOT.RooRealVar("tau_p", "tau_p", -10, 0)
    expo_p = ROOT.RooExponential("expo_p", "expo_p", x, tau_p)
    tau_f = ROOT.RooRealVar("tau_f", "tau_f", -10, 0)
    expo_f = ROOT.RooExponential("expo_f", "expo_f", x, tau_f)

    # PDFs from MC datasets
    # ---------------------
    pdf_mc_pass = ROOT.RooHistPdf("pdf_mc_pass", "pdf_mc_pass", x, h_mc[1])
    pdf_mc_fail = ROOT.RooHistPdf("pdf_mc_fail", "pdf_mc_fail", x, h_mc[0])

    # Smearing functions
    # ------------------
    mean = ROOT.RooRealVar("mean", "mean", 0, -2, 2)
    sigma = ROOT.RooRealVar("sigma", "sigma", 0.5, 0.001, 2)
    smear = ROOT.RooGaussian("smear", "smear", x, mean, sigma)
    if same_smearing is not True:
        mean_f = ROOT.RooRealVar("mean_f", "mean_f", 0, -2, 2)
        sigma_f = ROOT.RooRealVar("sigma_f", "sigma_f", 0.5, 0.001, 2)
        smear_fail = ROOT.RooGaussian(
            "smear_fail", "smear_fail", x, mean_f, sigma_f)
    else:
        smear_fail = smear

    # Convolutions with FFT
    # ---------------------
    x.setBins(1000, "cache")
    conv_pass = ROOT.RooFFTConvPdf("conv", "conv", x, pdf_mc_pass, smear, 3)
    conv_pass.setBufferFraction(0.1)
    conv_fail = ROOT.RooFFTConvPdf(
        "conv", "conv", x, pdf_mc_fail, smear_fail, 3)
    conv_fail.setBufferFraction(0.1)

    exp_ntot = n_events[0][0]+n_events[1][1]
    Nsig_tot = ROOT.RooRealVar("Nsig_tot", "#signal events total", exp_ntot,
                               0, exp_ntot+3*ROOT.TMath.Sqrt(exp_ntot))
    efficiency = ROOT.RooRealVar("efficiency", "efficiency", 0, 1)

    Nsig_pass = ROOT.RooProduct(
        "Nsig_pass", "Nsig_pass", [efficiency, Nsig_tot])
    Nbkg_pass = ROOT.RooRealVar(
            "nbkg_p", "#background events pass", 0, n_events[0][1])

    sum_pass = ROOT.RooAddPdf("sum_pass", "sum_pass",
                              [conv_pass, expo_p], [Nsig_pass, Nbkg_pass])
    model_pass = ROOT.RooAddPdf(sum_pass)

    one_minus_eff = ROOT.RooPolyVar(
        "one_minus_eff", "one_minus_eff", efficiency, [1, -1.])
    Nsig_fail = ROOT.RooProduct("prod", "prod", [one_minus_eff, Nsig_tot])

    Nbkg_fail = ROOT.RooRealVar(
        "nbkg_f", "#background events fail", 0, n_events[0][0])

    sum_fail = ROOT.RooAddPdf("sum_fail", "sum_fail",
                              [conv_fail, expo_f], [Nsig_fail, Nbkg_fail])
    model_fail = ROOT.RooAddPdf(sum_fail)

    sample = ROOT.RooCategory("sample", "sample")
    sample.defineType("pass")
    sample.defineType("fail")

    comb_dataset = ROOT.RooDataHist(
        "combData", "combined datasets", ROOT.RooArgSet(x), Index=sample,
        Import={"pass": h_data[1], "fail": h_data[0]})

    simPdf = ROOT.RooSimultaneous("simPdf", "simultaneous pdf", sample)
    simPdf.addPdf(model_pass, "pass")
    simPdf.addPdf(model_fail, "fail")

    fitResult = simPdf.fitTo(comb_dataset, Save=True,
                             PrintLevel=-1, Extended=True)

    results.add_result(fitResult, bin_pt, bin_eta)


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

    results = res_manager_sim()

    for bin_pt in range(15):
        for bin_eta in range(48):
            sim_efficiency(t, bin_pt+1, bin_eta+1,
                           results, same_smearing=False)

    results.write("simult_eff_results.pkl")

    print("RISULTATI SCRITTI SU PICKLE FILE")
