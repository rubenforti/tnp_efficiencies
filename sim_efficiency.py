"""
"""

import ROOT
import time
import os
import sys
from results_utils import results_manager
from make_plots import plot_pass_and_fail
from utilities import import_pdf_library, fit_quality


def check_existing_fit(ws, bin):

    if type(ws.obj(f'simPDF_({bin[0]}|{bin[1]})')) is ROOT.RooSimultaneous:
        print("Not possible to refit an existing PDF! \nReturning the results obtained previously")
        res = ROOT.RooFitResult(ws.obj(f'results_({bin[0]}|{bin[1]})'))
        return res
    else:
        return 0


def simultaneous_efficiency(type_eff, type_analysis, ws, bin, bkg_pdf,
                            same_smearing=False,
                            refit_numbkg=False,
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
    NBINS = 5000
    bufFractions = [0.2, 0.2]  # [Fail, Pass]

    fit_strategy = 2

    # Lower limit for NSIG and upper limit for NBKG, in terms of total events in the histogram, for pass and fail respectively
    lims_num = [[0.5, 0.2], [0.5, 0.2]]

    axis = ws.var(f"x_sim_({bin[0]}|{bin[1]})")
    axis.setRange("fitRange", 60, 120)
    
    axis.setBinning(ROOT.RooUniformBinning(
            axis.getRange("fitRange")[0], axis.getRange("fitRange")[1], NBINS), "cache")

    if same_smearing is True:
        mean = ROOT.RooRealVar(
            f"mean_({bin[0]}|{bin[1]})", "mean", 0, -5.0, 5.0)
        sigma = ROOT.RooRealVar(
            f"sigma_({bin[0]}|{bin[1]})", "sigma", 0.5, 0.1, 5.0)
    else:
        mean_p = ROOT.RooRealVar(
            f"mean_pass_({bin[0]}|{bin[1]})", "mean", 0, -5.0, 5.0)
        sigma_p = ROOT.RooRealVar(
            f"sigma_pass_({bin[0]}|{bin[1]})", "sigma", 0.5, 0.1, 5.0)
        mean_f = ROOT.RooRealVar(
            f"mean_fail_({bin[0]}|{bin[1]})", "mean", 0, -5.0, 5.0)
        sigma_f = ROOT.RooRealVar(
            f"sigma_fail_({bin[0]}|{bin[1]})", "sigma", 0.5, 0.1, 5.0)
        
        mean, sigma = [mean_f, mean_p], [sigma_f, sigma_p]
        


    tau_p = ROOT.RooRealVar(f"tau_pass_({bin[0]}|{bin[1]})", "tau", 0.0, -5.0, 5.0)
    tau_f = ROOT.RooRealVar(f"tau_fail_({bin[0]}|{bin[1]})", "tau", 0.0, -5.0, 5.0)
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
    histo_data, histo_mc, pdf_mc = [0, 0], [0, 0], [0, 0]
    conv_pdf, background = [0, 0], [0, 0]

    if same_smearing:
        smearing = ROOT.RooGaussian(f"smearing_({bin[0]}|{bin[1]})", 
                                    "Gaussian smearing", axis, mean, sigma)
    else:
        smearing = [0, 0]

    nexp, Nsig, Nbkg = [0, 0], [0, 0], [0, 0]
    sum_func, model, results, fit_status = [0, 0], [0, 0], [0, 0], [0, 0]

    for cond in ["pass", "fail"]:  # PASS has index 1, FAIL has index 0

        idx = 1 if cond == "pass" else 0
        print(f'idx = {idx}')

        histo_data[idx] = ws.data(f"Minv_data_{cond}_({bin[0]}|{bin[1]})")
        print(type(histo_data[idx]))

        histo_mc[idx] = ws.data(f"Minv_mc_{cond}_({bin[0]}|{bin[1]})")

        pdf_mc[idx] = ROOT.RooHistPdf(f"pdf_mc_{cond}_({bin[0]}|{bin[1]})", f"pdf_mc_{cond}",
                                      axis, histo_mc[idx], 3)
        print(type(pdf_mc[idx]))

        if same_smearing is True:
            conv_pdf[idx] = ROOT.RooFFTConvPdf(f"conv_{cond}_({bin[0]}|{bin[1]})", f"Convolution pdf",
                                               axis, pdf_mc[idx], smearing)
        else:
            smearing[idx] = ROOT.RooGaussian(
                f"smearing_{cond}_({bin[0]}|{bin[1]})", "Gaussian smearing",
                axis, mean[idx], sigma[idx])
            conv_pdf[idx] = ROOT.RooFFTConvPdf(f"conv_{cond}_({bin[0]}|{bin[1]})", f"Convolution pdf",
                                               axis, pdf_mc[idx], smearing[idx])

        conv_pdf[idx].setBufferFraction(bufFractions[idx])
        # conv_pdf[idx].setBufferStrategy(2)
        # print(type(conv_pdf_pass))

        if bkg_pdf == 'expo':
            background[idx] = ROOT.RooExponential(
                f"expo_bkg_{cond}_({bin[0]}|{bin[1]})", "Exponential bkg", axis, tau[idx])
        elif bkg_pdf == 'mixed' and idx == 0:
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
        expected_ntot_mc = histo_mc[1].sumEntries() + histo_mc[0].sumEntries()  # CHECK !!!!!
        Ntot_mc = ROOT.RooRealVar(f"Nsig_tot_mc_({bin[0]}|{bin[1]})", "N_mc Total", expected_ntot_mc,
                                  expected_ntot_mc-5*ROOT.TMath.Sqrt(expected_ntot_mc), 
                                  expected_ntot_mc+5*ROOT.TMath.Sqrt(expected_ntot_mc))
        
        eff_mc = ROOT.RooRealVar(f"efficiency_mc_({bin[0]}|{bin[1]})", "Efficiency MC", 0.9, 0.0, 1.0)
        Npass_mc = ROOT.RooProduct(f"npass_mc_({bin[0]}|{bin[1]})", "Nmc pass", [eff_mc, Ntot_mc])

        one_minus_eff_mc = ROOT.RooPolyVar("one_minus_eff_mc_({bin[0]}|{bin[1]})", 
                                           "1 - efficiency (mc)", eff_mc, [1, -1.])
        Nfail_mc = ROOT.RooProduct(
            f"Nfail_mc_({bin[0]}|{bin[1]})", "Nmc fail", [one_minus_eff_mc, Ntot_mc])

        model_mc_pass = ROOT.RooExtendPdf(f"mc_pass_ext_({bin[0]}|{bin[1]})",
                                          "MC extended pass", pdf_mc[1], Npass_mc)
        model_mc_fail = ROOT.RooExtendPdf(f"mc_fail_ext_({bin[0]}|{bin[1]})",
                                          "MC extended fail", pdf_mc[0], Nfail_mc)

        scale_factor = ROOT.RooRealVar(
            f"scale_factor_({bin[0]}|{bin[1]})", "Scale factor", 1.0, 0.0, 2.0)

        efficiency = ROOT.RooProduct(f"efficiency_({bin[0]}|{bin[1]})",
                                     "Efficiency", [scale_factor, eff_mc])
    else:
        efficiency = ROOT.RooRealVar(f"efficiency_({bin[0]}|{bin[1]})",
                                     "Efficiency", 0.9, 0, 1)

    # -----------------------
    #  Model for PASS events
    # -----------------------
    expected_ntot = histo_data[1].sumEntries() + histo_data[0].sumEntries()
    Nsig_tot = ROOT.RooRealVar(f"ntot_({bin[0]}|{bin[1]})", "Nsig Total", expected_ntot, 
                           expected_ntot-5*ROOT.TMath.Sqrt(expected_ntot), 
                           expected_ntot+5*ROOT.TMath.Sqrt(expected_ntot))

    Nsig_pass = ROOT.RooProduct(
            f"Nsig_pass_({bin[0]}|{bin[1]})", "Nsig pass", ROOT.RooArgList(efficiency, Nsig_tot))
    Nbkg_pass = ROOT.RooRealVar(
            f"Nbkg_pass_({bin[0]}|{bin[1]})", "Nbkg pass",
            0.1*histo_data[1].sumEntries(), 0.5, 1.5*histo_data[1].sumEntries())

    sum_pass = ROOT.RooAddPdf(
            f"sum_pass_({bin[0]}|{bin[1]})", "Signal+Bkg pass",
            ROOT.RooArgList(conv_pdf[1], background[1]), ROOT.RooArgList(Nsig_pass, Nbkg_pass))
    model_pass = ROOT.RooAddPdf(sum_pass)
    print(type(model_pass))

    # -----------------------
    #  Model for FAIL events
    # -----------------------
    roo_one = ROOT.RooRealVar("one", "one", 1.0)
    roo_minus_one = ROOT.RooRealVar("minus_one", "minus one", -1.0)

    one_minus_eff = ROOT.RooPolyVar(f"one_minus_eff_({bin[0]}|{bin[1]})",
                                    "1 - efficiency", efficiency, ROOT.RooArgList(roo_one, roo_minus_one))
    Nsig_fail = ROOT.RooProduct(
            f"Nsig_fail_({bin[0]}|{bin[1]})", "Nsig fail", ROOT.RooArgList(one_minus_eff, Nsig_tot))
    Nbkg_fail = ROOT.RooRealVar(
            f"Nbkg_fail_({bin[0]}|{bin[1]})", "Nbkg fail",
            0.1*histo_data[0].sumEntries(), 0.5, 1.5*histo_data[0].sumEntries())

    sum_fail = ROOT.RooAddPdf(
            f"sum_fail_({bin[0]}|{bin[1]})", "Signal+Bkg fail",
            ROOT.RooArgList(conv_pdf[0], background[0]), ROOT.RooArgList(Nsig_fail, Nbkg_fail))
    model_fail = ROOT.RooAddPdf(sum_fail)
    print(type(model_fail))

    # --------------------
    #  Categories and fit
    # --------------------
    sample = ROOT.RooCategory(f"sample_({bin[0]}{bin[1]})", "Composite sample")
    sample.defineType("pass")
    sample.defineType("fail")
    if enable_mcfit:
        sample.defineType("mc_pass")
        sample.defineType("mc_fail")
        comb_dataset = ROOT.RooDataHist(
                f"combData_({bin[0]}|{bin[1]})", "Combined datasets",
                ROOT.RooArgSet(axis), ROOT.RooFit.Index(sample),
                ROOT.RooFit.Import("pass", histo_data[1]),
                ROOT.RooFit.Import("fail", histo_data[0]),
                ROOT.RooFit.Import("mc_pass", histo_mc[1]), 
                ROOT.RooFit.Import("mc_fail", histo_mc[0]))
    else:
        comb_dataset = ROOT.RooDataHist(
                f"combData_({bin[0]}|{bin[1]})", "Combined datasets",
                ROOT.RooArgSet(axis), ROOT.RooFit.Index(sample),
                ROOT.RooFit.Import("pass", histo_data[1]),
                ROOT.RooFit.Import("fail", histo_data[0]))

    simPdf = ROOT.RooSimultaneous(
            f"simPDF_({bin[0]}|{bin[1]})", "Simultaneous pdf", sample)
    simPdf.addPdf(model_pass, "pass")
    simPdf.addPdf(model_fail, "fail")
    if enable_mcfit:
        simPdf.addPdf(model_mc_pass, "mc_pass")
        simPdf.addPdf(model_mc_fail, "mc_fail")

    res = simPdf.fitTo(comb_dataset,
                       ROOT.RooFit.Range("fitRange"),
                       ROOT.RooFit.Minimizer("Minuit2"),
                       ROOT.RooFit.Strategy(fit_strategy),
                       # ROOT.RooFit.MaxCalls(100000),
                       ROOT.RooFit.Save(1),
                       ROOT.RooFit.PrintLevel(verb))

    res.SetName(f"results_({bin[0]}|{bin[1]})")

    # pearson_chi2_eval(histo_data, model, histo_data.numEntries(), res)

    
    if fit_quality(res, old_checks=True) is True:
        ws.Import(simPdf)
        ws.Import(res)
        if figs:
            bkg_names =[background[0].GetName(), background[1].GetName()]
            plot_pass_and_fail((axis, axis), histo_data, (model_fail, model_pass), bkg_names, 
                               name=f"figs/fit_iso/bin_{bin[0]},{bin[1]}.pdf")
    else:
        bkg_names =[background[0].GetName(), background[1].GetName()]
        plot_pass_and_fail((axis, axis), histo_data, (model_fail, model_pass), bkg_names, 
                           name=f"figs/check_fits/bin_{bin[0]},{bin[1]}.pdf")

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
