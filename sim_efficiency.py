"""
"""

import ROOT
import time
import os
import sys
from results_utils import results_manager
from make_plots import plot_distr_with_fit
from utilities import import_pdf_library, pearson_chi2_eval, llr_test_bkg


def simultaneous_efficiency(type_eff, type_analysis, ws, bin, bkg_pdf,
                            same_smearing=True,
                            test_bkg=False,
                            enable_mcfit=False,
                            verb=-1, figs=False):
    """
    """

    path = os.path.dirname(__file__)
    ROOT.gSystem.cd(path)

    # ----------------------------------
    #  Initial parameters - MODIFY HERE
    # ----------------------------------
    NBINS = [2000, 10000]
    bufFractions = [0.5, 0.05]

    fit_strategy = [2, 2]

    # Lower limit for NSIG and upper limit for NBKG, in terms of total events in the histogram, for pass and fail respectively
    lims_num = [[0.5, 0.2], [0.5, 0.2]]


    axis = ws.var(f"x_sim_({bin[0]}|{bin[1]})")
    axis.setRange("fitRange", 50, 130)
    
    axis.setBinning(
        ROOT.RooUniformBinning(axis[idx].getRange("fitRange")[0],
                               axis[idx].getRange("fitRange")[1], NBINS[idx]), "cache")

    if same_smearing is True:
        mean = ROOT.RooRealVar(
            f"mean_({bin[0]}|{bin[1]})", "mean", 0, -2, 2)
        sigma = ROOT.RooRealVar(
            f"sigma_({bin[0]}|{bin[1]})", "sigma", 0.5, 0.01, 2)
    else:
        mean_p = ROOT.RooRealVar(
            f"mean_pass_({bin[0]}|{bin[1]})", "mean", 0, -2, 2)
        sigma_p = ROOT.RooRealVar(
            f"sigma_pass_({bin[0]}|{bin[1]})", "sigma", 0.5, 0.01, 2)
        mean_f = ROOT.RooRealVar(
            f"mean__fail_({bin[0]}|{bin[1]})", "mean", 0, -2, 2)
        sigma_f = ROOT.RooRealVar(
            f"sigma_fail_({bin[0]}|{bin[1]})", "sigma", 0.5, 0.01, 2)
        
        mean, sigma = [mean_f, mean_p], [sigma_f, sigma_p]


    tau_p = ROOT.RooRealVar(f"tau_pass_({bin[0]}|{bin[1]})", "tau", 0.0, -5, 5)
    tau_f = ROOT.RooRealVar(f"tau_fail_({bin[0]}|{bin[1]})", "tau", 0.0, -5, 5)
    tau = [tau_f, tau_p]
    
    alpha_p = ROOT.RooRealVar(f"alpha_pass_({bin[0]}|{bin[1]})", "alpha", 60.0, 40.0, 130.0)
    beta_p  = ROOT.RooRealVar(f"beta_pass_({bin[0]}|{bin[1]})", "beta",   2.5,  0.1,  40.0)
    gamma_p = ROOT.RooRealVar(f"gamma_pass_({bin[0]}|{bin[1]})", "gamma", 0.1,  0.0,  1.0)
    peak_p  = ROOT.RooRealVar(f"peak_pass_({bin[0]}|{bin[1]})", "peak",   90.0)

    alpha_f = ROOT.RooRealVar(f"alpha_fail_({bin[0]}|{bin[1]})", "alpha", 60.0, 40.0, 130.0)
    beta_f  = ROOT.RooRealVar(f"beta_fail_({bin[0]}|{bin[1]})", "beta",   2.5,  0.1,  40.0)
    gamma_f = ROOT.RooRealVar(f"gamma_fail_({bin[0]}|{bin[1]})", "gamma", 0.1,  0.0,  1.0)
    peak_f  = ROOT.RooRealVar(f"peak_fail_({bin[0]}|{bin[1]})", "peak",   90.0)

    alpha, beta, gamma, peak = [alpha_f, alpha_p], [beta_f, beta_p], [gamma_f, gamma_p], [peak_f, peak_p]


    # ------------------------------
    #  Histograms and PDFs building
    # ------------------------------
    mean, sigma = [mean_f, mean_p], [sigma_f, sigma_p]
    h_data, h_mc, pdf_mc = [0, 0], [0, 0], [0, 0]
    smearing, conv_pdf, background = [0, 0], [0, 0], [0, 0]

    nexp, Nsig, Nbkg = [0, 0], [0, 0], [0, 0]
    sum_func, model, results, fit_status = [0, 0], [0, 0], [0, 0], [0, 0]

    for cond in ["pass", "fail"]:  # PASS has index 1, FAIL has index 0

        idx = 1 if cond == "pass" else 0
        print(f'idx = {idx}')

        h_data[idx] = ws.data(f"Minv_data_{cond}_({bin[0]}|{bin[1]})")
        print(type(h_data[idx]))

        h_mc[idx] = ws.data(f"Minv_mc_{cond}_({bin[0]}|{bin[1]})")

        pdf_mc[idx] = ROOT.RooHistPdf(
            f"pdf_mc_{cond}_({bin[0]}|{bin[1]})", f"pdf_mc_{cond}",
            axis, h_mc[idx])
        print(type(pdf_mc[idx]))

        if same_smearing and idx == 1:  # idx==1 needed not to create two times the same pdfs
            smearing = ROOT.RooGaussian(
                f"smearing_({bin[0]}|{bin[1]})", "Gaussian smearing", axis, mean, sigma)
            conv_pdf[0] = ROOT.RooFFTConvPdf(f"conv_fail_({bin[0]}|{bin[1]})", f"Convolution pdf",
                                               axis, pdf_mc[0], smearing, 3)
            conv_pdf[1] = ROOT.RooFFTConvPdf(f"conv_pass_({bin[0]}|{bin[1]})", f"Convolution pdf",
                                               axis, pdf_mc[1], smearing, 3)
        else:
            smearing[idx] = ROOT.RooGaussian(
                f"smearing_{cond}_({bin[0]}|{bin[1]})", "Gaussian smearing",
                axis, mean[idx], sigma[idx])
            conv_pdf[idx] = ROOT.RooFFTConvPdf(f"conv_{cond}_({bin[0]}|{bin[1]})", f"Convolution pdf",
                                               axis, pdf_mc[idx], smearing[idx], 3)

        conv_pdf[idx].setBufferFraction(bufFractions[idx])
        # conv_pdf[idx].setBufferStrategy(2)
        # print(type(conv_pdf_pass))

        if bkg_pdf == 'expo':
            background[idx] = ROOT.RooExponential(
                f"expo_bkg_{cond}_({bin[0]}|{bin[1]})", "Exponential bkg", axis, tau[idx])
        elif bkg_pdf == 'mixed' and idx == 1:
            background[1] = ROOT.RooExponential(
                f"expo_bkg_pass_({bin[0]}|{bin[1]})", "Exponential bkg", axis, tau[1])
            background[0] = ROOT.RooCMSShape(f"cmsshape_bkg_pass_({bin[0]}|{bin[1]})", "CMSShape bkg", 
                                             axis, alpha[0], beta[0], gamma[0], peak[0])
        elif bkg_pdf == 'cmsshape':
            background[idx] = ROOT.RooCMSShape(f"cmsshape_bkg_{cond}_({bin[0]}|{bin[1]})", "CMSShape bkg", 
                                               axis, alpha[idx], beta[idx], gamma[idx], peak[idx])
        else:
            print("BKG shape given is not implemented! Retry with 'expo', 'cmsshape', or 'mixed'")
            sys.exit()

    # ------------------------------
    #  Model for extended fit on MC
    # ------------------------------
    if enable_mcfit:
        expected_ntot_mc = h_mc[1].sumEntries() + h_mc[0].sumEntries()  # CHECK !!!!!
        Nmc_tot = ROOT.RooRealVar(
            f"Nmc_tot_({bin[0]}|{bin[1]})", "N_mc Total",
            expected_ntot_mc, 0, expected_ntot_mc+5*ROOT.TMath.Sqrt(expected_ntot_mc))
        
        eff_mc = ROOT.RooRealVar(f"efficiency_mc_({bin[0]}|{bin[1]})", "Efficiency MC", 0, 1)
        Nmc_pass = ROOT.RooProduct(f"Nmc_pass_({bin[0]}|{bin[1]})", "Nmc pass", [eff_mc, Nmc_tot])

        one_minus_eff_mc = ROOT.RooPolyVar("one_minus_eff_mc_({bin[0]}|{bin[1]})", 
                                           "1 - efficiency (mc)", eff_mc, [1, -1.])
        Nmc_fail = ROOT.RooProduct(
            f"Nmc_fail_({bin[0]}|{bin[1]})", "Nmc fail", [one_minus_eff_mc, Nmc_tot])

        model_mc_pass = ROOT.RooExtendPdf(f"mc_pass_ext_({bin[0]}|{bin[1]})",
                                          "MC extended pass", pdf_mc[1], Nmc_pass)
        model_mc_fail = ROOT.RooExtendPdf(f"mc_fail_ext_({bin[0]}|{bin[1]})",
                                          "MC extended fail", pdf_mc[0], Nmc_fail)

    # -----------------------
    #  Model for PASS events
    # -----------------------
    expected_ntot = h_data[1].sumEntries() + h_data[0].sumEntries()
    Nsig_tot = ROOT.RooRealVar(
        f"Nsig_tot_({bin[0]}|{bin[1]})", "Nsig Total",
        expected_ntot, 0, expected_ntot+4*ROOT.TMath.Sqrt(expected_ntot))

    if enable_mcfit:
        scale_factor = ROOT.RooRealVar(
            f"scale_factor_({bin[0]}|{bin[1]})", "Scale factor", 0, 2)

        efficiency = ROOT.RooProduct(f"efficiency_({bin[0]}|{bin[1]})",
                                     "Efficiency", [scale_factor, eff_mc])
    else:
        efficiency = ROOT.RooRealVar(f"efficiency_({bin[0]}|{bin[1]})",
                                     "Efficiency", 0, 1)

    Nsig_pass = ROOT.RooProduct(
            f"Nsig_pass_({bin[0]}|{bin[1]})", "Nsig pass", [efficiency, Nsig_tot])
    Nbkg_pass = ROOT.RooRealVar(
            f"Nbkg_pass_({bin[0]}|{bin[1]})", "Nbkg pass",
            0.001*h_data[1].sumEntries(), 0, 0.5*h_data[1].sumEntries())

    sum_pass = ROOT.RooAddPdf(
            f"sum_pass_({bin[0]}|{bin[1]})", "Signal+Bkg pass",
            [conv_pdf_pass, background_pass], [Nsig_pass, Nbkg_pass])
    model_pass = ROOT.RooAddPdf(sum_pass)
    print(type(model_pass))

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
            f"sum_fail_({bin[0]}|{bin[1]})", "Signal+Bkg fail",
            [conv_pdf_fail, background_fail], [Nsig_fail, Nbkg_fail])
    model_fail = ROOT.RooAddPdf(sum_fail)
    print(type(model_fail))

    # --------------------
    #  Categories and fit
    # --------------------
    sample = ROOT.RooCategory(
            f"sample_({bin[0]}{bin[1]})", "Composite sample")
    sample.defineType("pass")
    sample.defineType("fail")
    if enable_mcfit:
        sample.defineType("mc_pass")
        sample.defineType("mc_fail")
        comb_dataset = ROOT.RooDataHist(
                f"combData_({bin[0]}|{bin[1]})", "Combined datasets",
                ROOT.RooArgSet(axis), Index=sample,
                Import={"pass": h_data[1], "fail": h_data[0], "mc_pass": h_mc[1], "mc_fail": h_mc[0]})
    else:
        comb_dataset = ROOT.RooDataHist(
                f"combData_({bin[0]}|{bin[1]})", "Combined datasets",
                ROOT.RooArgSet(axis), Index=sample,
                Import={"pass": h_data[1], "fail": h_data[0]})

    simPdf = ROOT.RooSimultaneous(
            f"simPdf_({bin[0]}|{bin[1]})", "Simultaneous pdf", sample)
    simPdf.addPdf(model_pass, "pass")
    simPdf.addPdf(model_fail, "fail")
    if enable_mcfit:
        simPdf.addPdf(model_mc_pass, "mc_pass")
        simPdf.addPdf(model_mc_fail, "mc_fail")

    res = simPdf.fitTo(comb_dataset,
                       ROOT.RooFit.Range("fitRange"),
                       ROOT.RooFit.Minimizer("Minuit2"),
                       ROOT.RooFit.Strategy(fit_strategy[idx]),
                       # ROOT.RooFit.MaxCalls(100000),
                       ROOT.RooFit.Save(1),
                       ROOT.RooFit.PrintLevel(verb))

    res.SetName(f"results_{cond}_({bin[0]}|{bin[1]})")

    # pearson_chi2_eval(histo_data, model, histo_data.numEntries(), res)

    '''
    if fit_quality(res) is True:
        ws.Import(model)
        ws.Import(res)
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

    results = results_manager("sim")

    bin = (1, 1)
    simultaneous_efficiency(t, bin, same_smearing=True)

    results.write("simult_eff_results.pkl")

    print("RISULTATI SCRITTI SU PICKLE FILE")
