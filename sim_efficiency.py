"""
"""

import ROOT
import time
import os
import sys
from results_utilities import res_manager_sim, fit_quality
from plot_functions import makeAndSavePlot
from stat_functions import pearson_chi2_eval, llr_test_bkg
from utilities import import_pdf_library


def sim_efficiency(type_eff, workspace, bin, bkg_pdf,
                   same_smearing=True,
                   test_bkg=False,
                   enable_mcfit=False,
                   saveplots=False):
    """
    """

    if type(workspace[f'PDF_pass_({bin[0]},{bin[1]})']) is ROOT.RooAddPdf:
        # DA IMPLEMENTARE CORRETTAMENTE
        sys.exit()

    elif type(workspace[f'PDF_pass_({bin[0]},{bin[1]})']) is ROOT.TObject and type(workspace[f'PDF_fail_({bin[0]},{bin[1]})']) is ROOT.TObject:

        histo_data_pass = workspace[f"Minv_data_pass_({bin[0]},{bin[1]})"]
        # histo_mc_pass = workspace[f"Minv_mc_pass_({bin[0]},{bin[1]})"]

        histo_data_fail = workspace[f"Minv_data_fail_({bin[0]},{bin[1]})"]
        # histo_mc_fail = workspace[f"Minv_mc_fail_({bin[0]},{bin[1]})"]

        axis = workspace[f"x_sim_({bin[0]},{bin[1]})"]

        # Backgrounds
        # -----------
        tau_p = ROOT.RooRealVar("tau_p", "tau_p", -10, 0)
        expo_p = ROOT.RooExponential("expo_p", "expo_p", axis, tau_p)
        tau_f = ROOT.RooRealVar("tau_f", "tau_f", -10, 0)
        expo_f = ROOT.RooExponential("expo_f", "expo_f", axis, tau_f)

        # PDFs from MC datasets
        # ---------------------
        pdf_mc_pass = ROOT.RooHistPdf(
            "pdf_mc_pass", "pdf_mc_pass", axis,
            workspace[f"Minv_mc_pass_({bin[0]},{bin[1]})"])
        pdf_mc_fail = ROOT.RooHistPdf(
            "pdf_mc_fail", "pdf_mc_fail", axis,
            workspace[f"Minv_mc_fail_({bin[0]},{bin[1]})"])

        # Smearing functions and FFT convolutions
        # ---------------------------------------
        mean = ROOT.RooRealVar("mean", "mean", 0, -2, 2)
        sigma = ROOT.RooRealVar("sigma", "sigma", 0.5, 0.001, 2)
        smear = ROOT.RooGaussian("smear", "smear", axis, mean, sigma)

        axis.setBins(1000, "cache")

        conv_pass = ROOT.RooFFTConvPdf(
            "conv", "conv", axis, pdf_mc_pass, smear, 3)

        if same_smearing is not True:
            mean_f = ROOT.RooRealVar("mean_f", "mean_f", 0, -2, 2)
            sigma_f = ROOT.RooRealVar("sigma_f", "sigma_f", 0.5, 0.001, 2)

            smear_fail = ROOT.RooGaussian(
                "smear_fail", "smear_fail", axis, mean_f, sigma_f)

            conv_fail = ROOT.RooFFTConvPdf(
                "conv", "conv", axis, pdf_mc_fail, smear_fail, 3)
            conv_fail.setBufferFraction(0.1)
        else:
            conv_fail = ROOT.RooFFTConvPdf(
                "conv", "conv", axis, pdf_mc_fail, smear, 3)

        conv_pass.setBufferFraction(0.1)
        conv_fail.setBufferFraction(0.1)

        expected_ntot = histo_data_pass.sumEntries()+histo_data_fail.sumEntries()
        Nsig_tot = ROOT.RooRealVar(
            "Nsig_tot", "#signal events total",
            expected_ntot, 0, expected_ntot+3*ROOT.TMath.Sqrt(expected_ntot))

        efficiency = ROOT.RooRealVar("efficiency", "efficiency", 0, 1)

        Nsig_pass = ROOT.RooProduct(
            "Nsig_pass", "Nsig_pass", [efficiency, Nsig_tot])
        Nbkg_pass = ROOT.RooRealVar(
                "nbkg_p", "#background events pass", 0, histo_data_pass.sumEntries())

        sum_pass = ROOT.RooAddPdf("sum_pass", "sum_pass",
                                  [conv_pass, expo_p], [Nsig_pass, Nbkg_pass])
        model_pass = ROOT.RooAddPdf(sum_pass)

        one_minus_eff = ROOT.RooPolyVar(
            "one_minus_eff", "one_minus_eff", efficiency, [1, -1.])
        Nsig_fail = ROOT.RooProduct("prod", "prod", [one_minus_eff, Nsig_tot])

        Nbkg_fail = ROOT.RooRealVar(
            "nbkg_f", "#background events fail", 0, histo_data_fail.sumEntries())

        sum_fail = ROOT.RooAddPdf("sum_fail", "sum_fail",
                                  [conv_fail, expo_f], [Nsig_fail, Nbkg_fail])
        model_fail = ROOT.RooAddPdf(sum_fail)

        sample = ROOT.RooCategory("sample", "sample")
        sample.defineType("pass")
        sample.defineType("fail")

        comb_dataset = ROOT.RooDataHist(
            "combData", "combined datasets", ROOT.RooArgSet(axis), Index=sample,
            Import={"pass": histo_data_pass, "fail": histo_data_fail})

        simPdf = ROOT.RooSimultaneous("simPdf", "simultaneous pdf", sample)
        simPdf.addPdf(model_pass, "pass")
        simPdf.addPdf(model_fail, "fail")

        fitResult = simPdf.fitTo(comb_dataset, Save=True,
                                 PrintLevel=-1, Extended=True)

        print(fitResult.floatParsFinal())


if __name__ == '__main__':

    path = os.path.dirname(__file__)
    ROOT.gSystem.cd(path)

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

    for bin_pt in range(1):
        for bin_eta in range(1):
            sim_efficiency(t, bin_pt+1, bin_eta+1, same_smearing=False)

    results.write("simult_eff_results.pkl")

    print("RISULTATI SCRITTI SU PICKLE FILE")
