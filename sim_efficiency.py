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

        axis = workspace[f"x_sim_({bin[0]},{bin[1]})"]
        axis.setBins(3000, "cache")

        h_data, pdf_mc = [0, 0], [0, 0]

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

            h_data[idx] = workspace[f"Minv_data_{cond}_({bin[0]},{bin[1]})"]
            print(type(h_data[idx]))

            pdf_mc[idx] = ROOT.RooHistPdf(
                f"pdf_mc_{cond}_({bin[0]},{bin[1]})", f"pdf_mc_{cond}",
                axis, workspace[f"Minv_mc_{cond}_({bin[0]},{bin[1]})"])
            print(type(pdf_mc[idx]))

            if same_smearing is True:
                smearing[idx] = ROOT.RooGaussian(
                    f"smearing_({bin[0]},{bin[1]})", "Gaussian smearing", axis, mean, sigma)
            else:
                mean[idx] = ROOT.RooRealVar(f"mean_{cond}_({bin[0]},{bin[1]})",
                                            "mean", 0, -2, 2)
                sigma[idx] = ROOT.RooRealVar(f"sigma_{cond}_({bin[0]},{bin[1]})",
                                             "sigma", 0.5, 0.001, 2)
                smearing[idx] = ROOT.RooGaussian(
                    f"smearing_{cond}_({bin[0]},{bin[1]})", "Gaussian smearing",
                    axis, mean[idx], sigma[idx])

            conv_pdf[idx] = ROOT.RooFFTConvPdf(
                f"conv_{cond}_({bin[0]}_{bin[1]})", f"Convolution {cond}",
                axis, pdf_mc[idx], smearing[idx], 3)
            conv_pdf[idx].setBufferFraction(0.1)

        # Backgrounds
        # -----------
        if bkg_pdf == 'expo':
            tau_pass = ROOT.RooRealVar(
                f"tau_pass_({bin[0]},{bin[1]})", "tau", -10, 0)
            expo_pass = ROOT.RooExponential(
                f"expo_bkg_pass_({bin[0]},{bin[1]})",
                "Exponential background", axis, tau_pass)
            tau_fail = ROOT.RooRealVar(
                f"tau_fail_({bin[0]},{bin[1]})", "tau", -10, 0)
            expo_fail = ROOT.RooExponential(
                f"expo_bkg_fail_({bin[0]},{bin[1]})",
                "Exponential background", axis, tau_fail)

        expected_ntot = h_data[1].sumEntries()+h_data[0].sumEntries()
        Nsig_tot = ROOT.RooRealVar(
            f"Nsig_tot_({bin[0]},{bin[1]})", "Nsig total",
            expected_ntot, 0, expected_ntot+3*ROOT.TMath.Sqrt(expected_ntot))

        efficiency = ROOT.RooRealVar(
            f"efficiency_({bin[0]},{bin[1]})", "Efficiency", 0, 1)

        Nsig_pass = ROOT.RooProduct(f"Nsig_pass_({bin[0]},{bin[1]})",
                                    "Nsig pass", [efficiency, Nsig_tot])
        Nbkg_pass = ROOT.RooRealVar(f"Nbkg_pass_({bin[0]},{bin[1]})",
                                    "Nbkg pass", 0, h_data[1].sumEntries())

        sum_pass = ROOT.RooAddPdf(f"sum_pass_({bin[0]},{bin[1]})", "Signal+Bkg",
                                  [conv_pdf[1], expo_pass], [Nsig_pass, Nbkg_pass])
        model_pass = ROOT.RooAddPdf(sum_pass)

        one_minus_eff = ROOT.RooPolyVar(f"one_minus_eff_({bin[0]},{bin[1]})",
                                        "1 - efficiency", efficiency, [1, -1.])

        Nsig_fail = ROOT.RooProduct(f"Nsig_fail_({bin[0]},{bin[1]})",
                                    "Nsig fail", [one_minus_eff, Nsig_tot])
        Nbkg_fail = ROOT.RooRealVar(f"Nbkg_fail_({bin[0]},{bin[1]})",
                                    "Nsig fail",  0, h_data[0].sumEntries())

        sum_fail = ROOT.RooAddPdf(f"sum_fail_({bin[0]},{bin[1]})", "Signal+Bkg",
                                  [conv_pdf[0], expo_fail], [Nsig_fail, Nbkg_fail])
        model_fail = ROOT.RooAddPdf(sum_fail)

        sample = ROOT.RooCategory(
            f"sample_({bin[0]},{bin[1]})", "Composite sample")
        sample.defineType("pass")
        sample.defineType("fail")

        comb_dataset = ROOT.RooDataHist(
            f"combData_({bin[0]},{bin[1]})", "Combined datasets", ROOT.RooArgSet(
                axis),
            Index=sample, Import={"pass": h_data[1], "fail": h_data[0]})

        simPdf = ROOT.RooSimultaneous(
            f"simPdf_({bin[0]},{bin[1]})", "Simultaneous pdf", sample)
        simPdf.addPdf(model_pass, "pass")
        simPdf.addPdf(model_fail, "fail")

        res = simPdf.fitTo(comb_dataset, Save=True,
                           PrintLevel=0, Extended=True)

        print(res.status())
        print(res.covQual())
        print(res.edm())

        if fit_quality(res) is True:
            workspace.Import(simPdf)
            workspace.Import(res)


if __name__ == '__main__':

    path = os.path.dirname(__file__)
    ROOT.gSystem.cd(path)

    custom_pdfs = ['RooCBExGaussShape',
                   'RooDoubleCBFast', 'RooCMSShape', 'my_double_CB']
    import_pdf_library(custom_pdfs[2])

    type_eff = ("sa", "global", "ID", "iso", "trigger", "veto")
    t = type_eff[3]

    results = res_manager_sim()

    bin = (1, 1)
    simultaneous_efficiency(t, bin, 'expo', same_smearing=True)

    results.write("simult_eff_results.pkl")

    print("RISULTATI SCRITTI SU PICKLE FILE")
