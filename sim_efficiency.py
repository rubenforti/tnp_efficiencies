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


def simultaneous_efficiency(type_eff, bin, bkg_pdf,
                            same_smearing=True,
                            test_bkg=False,
                            enable_mcfit=False,
                            saveplots=False):
    """
    """

    path = os.path.dirname(__file__)
    ROOT.gSystem.cd(path)

    file_ws = ROOT.TFile(f"root_files/{type_eff}_workspace.root")
    workspace = file_ws.Get("w")

    if type(workspace[f'PDF_pass_({bin[0]},{bin[1]})']) is ROOT.RooAddPdf:
        # DA IMPLEMENTARE CORRETTAMENTE
        sys.exit()

    elif type(workspace[f'PDF_pass_({bin[0]},{bin[1]})']) is ROOT.TObject and type(workspace[f'PDF_fail_({bin[0]},{bin[1]})']) is ROOT.TObject:

        axis, h_data, pdf_mc = [0, 0], [0, 0], [0, 0]

        smearing, conv_pdf = [0, 0], [0, 0]

        if same_smearing is True:
            mean = ROOT.RooRealVar(
                f"mean_({bin[0]},{bin[1]})", "mean", 0, -2, 2)
            sigma = ROOT.RooRealVar(
                f"sigma_({bin[0]},{bin[1]})", "sigma", 0.5, 0.001, 2)
        else:
            mean, sigma = [0, 0], [0, 0]

        for cond in ["pass", "fail"]:

            idx = 1 if cond == "pass" else 0
            print(f'idx = {idx}')

            axis[idx] = workspace[f"x_{cond}_({bin[0]},{bin[1]})"]

            h_data[idx] = workspace[f"Minv_data_{cond}_({bin[0]},{bin[1]})"]
            print(type(h_data[idx]))

            pdf_mc[idx] = ROOT.RooHistPdf(
                f"pdf_mc_{cond}_({bin[0]},{bin[1]})", f"pdf_mc_{cond}",
                axis[idx], workspace[f"Minv_mc_{cond}_({bin[0]},{bin[1]})"])
            print(type(pdf_mc[idx]))

            if same_smearing is True:
                smearing[idx] = ROOT.RooGaussian(
                    f"smearing_({bin[0]},{bin[1]})", "Gaussian smearing", axis[idx], mean, sigma)
            else:
                mean[idx] = ROOT.RooRealVar(f"mean_{cond}_({bin[0]},{bin[1]})",
                                            "mean", 0, -2, 2)
                sigma[idx] = ROOT.RooRealVar(f"sigma_{cond}_({bin[0]},{bin[1]})",
                                             "sigma", 0.5, 0.001, 2)
                smearing[idx] = ROOT.RooGaussian(
                    f"smearing_{cond}_({bin[0]},{bin[1]})", "Gaussian smearing",
                    axis[idx], mean[idx], sigma[idx])

            axis[idx].setBins(3000, "cache")
            conv_pdf[idx] = ROOT.RooFFTConvPdf(
                f"conv_{cond}_({bin[0]}_{bin[1]})", f"Convolution {cond}",
                axis[idx], pdf_mc[idx], smearing[idx], 3)
            conv_pdf[idx].setBufferFraction(0.1)

        # Backgrounds
        # -----------
        tau_p = ROOT.RooRealVar("tau_p", "tau_p", -10, 0)
        expo_p = ROOT.RooExponential("expo_p", "expo_p", axis[1], tau_p)
        tau_f = ROOT.RooRealVar("tau_f", "tau_f", -10, 0)
        expo_f = ROOT.RooExponential("expo_f", "expo_f", axis[0], tau_f)

        expected_ntot = h_data[1].sumEntries()+h_data[0].sumEntries()
        Nsig_tot = ROOT.RooRealVar(
            "Nsig_tot", "#signal events total",
            expected_ntot, 0, expected_ntot+3*ROOT.TMath.Sqrt(expected_ntot))

        efficiency = ROOT.RooRealVar("efficiency", "efficiency", 0, 1)

        Nsig_pass = ROOT.RooProduct(
            "Nsig_pass", "Nsig_pass", [efficiency, Nsig_tot])
        Nbkg_pass = ROOT.RooRealVar(
                "nbkg_p", "#background events pass", 0, h_data[1].sumEntries())

        sum_pass = ROOT.RooAddPdf("sum_pass", "sum_pass",
                                  [conv_pdf[1], expo_p], [Nsig_pass, Nbkg_pass])
        model_pass = ROOT.RooAddPdf(sum_pass)

        one_minus_eff = ROOT.RooPolyVar(
            "one_minus_eff", "one_minus_eff", efficiency, [1, -1.])
        Nsig_fail = ROOT.RooProduct("prod", "prod", [one_minus_eff, Nsig_tot])

        Nbkg_fail = ROOT.RooRealVar(
            "nbkg_f", "#background events fail", 0, h_data[0].sumEntries())

        sum_fail = ROOT.RooAddPdf("sum_fail", "sum_fail",
                                  [conv_pdf[0], expo_f], [Nsig_fail, Nbkg_fail])
        model_fail = ROOT.RooAddPdf(sum_fail)

        sample = ROOT.RooCategory("sample", "sample")
        sample.defineType("pass")
        sample.defineType("fail")

        comb_dataset = ROOT.RooDataHist(
            "combData", "combined datasets", ROOT.RooArgSet(axis[0], axis[1]),
            Index=sample, Import={"pass": h_data[1], "fail": h_data[0]})

        simPdf = ROOT.RooSimultaneous("simPdf", "simultaneous pdf", sample)
        simPdf.addPdf(model_pass, "pass")
        simPdf.addPdf(model_fail, "fail")

        res = simPdf.fitTo(comb_dataset, Save=True,
                           PrintLevel=0, Extended=True)

        print(res.status())
        print(res.covQual())
        print(res.edm())

        '''
        if fit_quality(res) is True:
            workspace.Import(simPdf)
            workspace.Import(res)
        '''


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

    bin = (1, 1)
    simultaneous_efficiency(t, bin, 'expo', same_smearing=True)

    results.write("simult_eff_results.pkl")

    print("RISULTATI SCRITTI SU PICKLE FILE")
