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
                            verb=-1, figs=False):
    """
    """

    path = os.path.dirname(__file__)
    ROOT.gSystem.cd(path)

    file_ws = ROOT.TFile(f"root_files/{type_eff}_workspace.root")
    workspace = file_ws.Get("w")

    '''
    if type(workspace[f'PDF_pass_({bin[0]},{bin[1]})']) is ROOT.RooAddPdf:
        # DA IMPLEMENTARE CORRETTAMENTE
        sys.exit()

    elif type(workspace[f'PDF_pass_({bin[0]}|{bin[1]})']) is ROOT.TObject and type(workspace[f'PDF_fail_({bin[0]}|{bin[1]})']) is ROOT.TObject:
    '''

    # SI ASSUME CHE LE PDF NON SIANO GIÀ STATE DEFINITE !!!!!!!!!!!!

    axis = workspace[f"x_sim_({bin[0]}|{bin[1]})"]

    axis.setRange("fitRange", 50, 130)

    NBINS = 8000
    binning = ROOT.RooUniformBinning(axis.getRange(
        "fitRange")[0], axis.getRange("fitRange")[1], NBINS)
    axis.setBinning(binning, "cache")

    h_data, pdf_mc = [0, 0], [0, 0]

    smearing, conv_pdf = [0, 0], [0, 0]

    if same_smearing is True:
        mean = ROOT.RooRealVar(
            f"mean_({bin[0]},{bin[1]})", "MEAN", 0, -2, 2)
        sigma = ROOT.RooRealVar(
            f"sigma_({bin[0]},{bin[1]})", "SIGMA", 0.5, 0.001, 2)
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
                                        f"MEAN {cond}", 0, -2, 2)
            sigma[idx] = ROOT.RooRealVar(f"sigma_{cond}_({bin[0]},{bin[1]})",
                                         f"SIGMA {cond}", 0.5, 0.001, 2)
            smearing[idx] = ROOT.RooGaussian(
                f"smearing_{cond}_({bin[0]},{bin[1]})", "Gaussian smearing",
                axis, mean[idx], sigma[idx])

    conv_pdf_pass = ROOT.RooFFTConvPdf(
        f"conv_pass_({bin[0]}|{bin[1]})", f"Convolution pass",
        axis, pdf_mc[1], smearing[1], 3)
    conv_pdf_pass.setBufferFraction(0.1)
    conv_pdf_pass.setBufferStrategy(2)

    conv_pdf_fail = ROOT.RooFFTConvPdf(
        f"conv_pass_({bin[0]}|{bin[1]})", "Convolution fail",
        axis, pdf_mc[0], smearing[0], 3)
    conv_pdf_fail.setBufferFraction(0.1)
    conv_pdf_fail.setBufferStrategy(2)

    # -----------------
    #  Background pass
    # -----------------
    if (bkg_pdf == 'expo') or (bkg_pdf == 'mixed'):
        background_pass = ROOT.RooExponential(
            f"expo_bkg_pass_({bin[0]},{bin[1]})", "Exponential bkg",
            axis, ROOT.RooRealVar(f"tau_pass_({bin[0]},{bin[1]})",
                                  "TAU PASS", -0.5, -2.0, 0.001))
    elif bkg_pdf == 'cmsshape':
        background_pass = ROOT.RooCMSShape(
            f"cmsshape_bkg_pass_({bin[0]}|{bin[1]})", "CMSShape bkg", axis,
            ROOT.RooRealVar(
                f"alpha_pass_({bin[0]}|{bin[1]})", "ALPHA", 60.0, 40.0, 130.0),
            ROOT.RooRealVar(
                f"beta_pass_({bin[0]}|{bin[1]})", "BETA", 2.5, 0.01, 15.0),
            ROOT.RooRealVar(
                f"gamma_pass_({bin[0]}|{bin[1]})", "GAMMA", 0.1, 0.0001, 0.2),
            ROOT.RooRealVar(
                f"peak_pass_({bin[0]}|{bin[1]})", "PEAK", 90.0))
    else:
        print(
            "BKG shape given is not implemented! Retry with 'expo', 'cmsshape', or 'mixed'")
        sys.exit()

    # -----------------
    #  Background fail
    # -----------------
    if bkg_pdf == 'expo':
        background_fail = ROOT.RooExponential(
            f"expo_bkg_fail_({bin[0]},{bin[1]})", "Exponential bkg",
            axis, ROOT.RooRealVar(f"tau_fail_({bin[0]},{bin[1]})",
                                  "TAU FAIL", -0.5, -2.0, 0.001))
    elif (bkg_pdf == 'cmsshape') or (bkg_pdf == 'mixed'):
        background_fail = ROOT.RooCMSShape(
            f"cmsshape_bkg_fail_({bin[0]}|{bin[1]})", "CMSShape bkg", axis,
            ROOT.RooRealVar(
                f"alpha_fail_({bin[0]}|{bin[1]})", "ALPHA", 60.0, 40.0, 130.0),
            ROOT.RooRealVar(
                f"beta_fail_({bin[0]}|{bin[1]})", "BETA", 2.5, 0.01, 15.0),
            ROOT.RooRealVar(
                f"gamma_fail_({bin[0]}|{bin[1]})", "GAMMA", 0.1, 0.0001, 0.2),
            ROOT.RooRealVar(
                f"peak_fail_({bin[0]}|{bin[1]})", "PEAK", 90.0))

    # -----------------------
    #  Model for PASS events
    # -----------------------
    expected_ntot = h_data[1].sumEntries() + h_data[0].sumEntries()
    Nsig_tot = ROOT.RooRealVar(
        f"Nsig_tot_({bin[0]}|{bin[1]})", "Nsig Toral",
        expected_ntot, 0, expected_ntot+4*ROOT.TMath.Sqrt(expected_ntot))

    efficiency = ROOT.RooRealVar(
        f"efficiency_({bin[0]}|{bin[1]})", "Efficiency", 0, 1)

    Nsig_pass = ROOT.RooProduct(
            f"Nsig_pass_({bin[0]},{bin[1]})", "Nsig pass", [efficiency, Nsig_tot])
    Nbkg_pass = ROOT.RooRealVar(
            f"Nbkg_pass_({bin[0]},{bin[1]})", "Nbkg pass",
            0.001*h_data[1].sumEntries(), 0, 0.5*h_data[1].sumEntries())

    sum_pass = ROOT.RooAddPdf(
            f"sum_pass_({bin[0]},{bin[1]})", "Signal+Bkg pass",
            [conv_pdf[1], background_pass], [Nsig_pass, Nbkg_pass])
    model_pass = ROOT.RooAddPdf(sum_pass)

    # -----------------------
    #  Model for FAIL events
    # -----------------------
    one_minus_eff = ROOT.RooPolyVar(f"one_minus_eff_({bin[0]}|{bin[1]})",
                                    "1 - efficiency", efficiency, [1, -1.])

    Nsig_fail = ROOT.RooProduct(
            f"Nsig_fail_({bin[0]}|{bin[1]})", "Nsig fail", [one_minus_eff, Nsig_tot])
    Nbkg_fail = ROOT.RooRealVar(
            f"Nbkg_fail_({bin[0]}|{bin[1]})", "Nbkg fail",
            0.001*h_data[0].sumEntries(), 0, 0.5*h_data[0].sumEntries())

    sum_fail = ROOT.RooAddPdf(
            f"sum_fail_({bin[0]},{bin[1]})", "Signal+Bkg fail",
            [conv_pdf[0], background_fail], [Nsig_fail, Nbkg_fail])
    model_fail = ROOT.RooAddPdf(sum_fail)

    # --------------------
    #  Categories and fit
    # --------------------
    sample = ROOT.RooCategory(
            f"sample_({bin[0]}{bin[1]})", "Composite sample")
    sample.defineType("pass")
    sample.defineType("fail")

    comb_dataset = ROOT.RooDataHist(
            f"combData_({bin[0]}|{bin[1]})", "Combined datasets", ROOT.RooArgSet(
                axis),
            Index=sample, Import={"pass": h_data[1], "fail": h_data[0]})

    simPdf = ROOT.RooSimultaneous(
            f"simPdf_({bin[0]},{bin[1]})", "Simultaneous pdf", sample)
    simPdf.addPdf(model_pass, "pass")
    simPdf.addPdf(model_fail, "fail")

    res = simPdf.fitTo(comb_dataset,
                       Extended=True,
                       Range='fitRange',
                       # ExternalConstraints=ROOT.RooArgSet("nexp"),
                       # BatchMode=False,
                       # Offset=True,
                       Minimizer=("Minuit2", "migrad"),
                       Strategy=2,
                       MaxCalls=100000,
                       Save=True,
                       PrintLevel=verb)

    res.SetName(f"results_{cond}_({bin[0]}|{bin[1]})")

    # pearson_chi2_eval(histo_data, model, histo_data.numEntries(), res)

    '''
    if fit_quality(res) is True:
        workspace.Import(model)
        workspace.Import(res)
    '''

    '''
    else:
        print("******\nERROR in PDF types\n*******")
        sys.exit()
    '''

    '''
    DA ADATTARE !!!!!!
    if figs is True:
        makeAndSavePlot(axis, histo_data, model,
                        bkg_name=background.GetName(), pull=False,
                        name=f"{cond}_{bin[0]}_{bin[1]}.png")
    '''

    return res


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
