"""
"""

import ROOT
import time
import os
import sys
from results_utils import results_manager
from make_plots import plot_distr_with_fit
from utilities import import_pdf_library, fit_quality, fit_quality_old, pearson_chi2_eval


def check_existing_fit(ws, bin):
    isFittedPass = type(
        ws.obj(f'PDF_pass_({bin[0]}|{bin[1]})')) is ROOT.RooAddPdf
    isFittedFail = type(
        ws.obj(f'PDF_fail_({bin[0]}|{bin[1]})')) is ROOT.RooAddPdf
    if isFittedPass and isFittedFail:
        print("Not possible to refit an existing PDF! \nReturning the results obtained previously")
        res_pass = ROOT.RooFitResult(
            ws.obj(f'results_pass_({bin[0]}|{bin[1]})'))
        res_fail = ROOT.RooFitResult(
            ws.obj(f'results_fail_({bin[0]}|{bin[1]})'))
        return res_pass, res_fail
    else:
        return 0


def indep_eff_fits(type_eff, type_analysis, ws, bin, bkg_pdf, refit_numbkg=False, test_bkg=False,
                   verb=-1, figs=False):
  
    path = os.path.dirname(__file__)
    ROOT.gSystem.cd(path)

    # ----------------------------------
    #  Initial parameters - MODIFY HERE
    # ----------------------------------
    NBINS = [2000, 2000]
    bufFractions = [0.5, 0.5]

    fit_strategy = [2, 2]

    # Lower limit for NSIG and upper limit for NBKG, in terms of total events in the histogram, for pass and fail respectively
    lims_num = [[0.5, 0.2], [0.5, 0.2]]

    mean_p = ROOT.RooRealVar(
        f"mean_pass_({bin[0]}|{bin[1]})", "mean", 0, -5.0, 5.0)
    sigma_p = ROOT.RooRealVar(
        f"sigma_pass_({bin[0]}|{bin[1]})", "sigma", 2, 0.5, 6.0)
    mean_f = ROOT.RooRealVar(
        f"mean_fail_({bin[0]}|{bin[1]})", "mean", 0, -5.0, 5.0)
    sigma_f = ROOT.RooRealVar(
        f"sigma_fail_({bin[0]}|{bin[1]})", "sigma", 2, 0.5, 6.0)

    if bkg_pdf == "expo":
        tau_p = ROOT.RooRealVar(
            f"tau_pass_({bin[0]}|{bin[1]})", "tau", 0.0, -1.5, 1.5)
        tau_f = ROOT.RooRealVar(
            f"tau_fail_({bin[0]}|{bin[1]})", "tau", 0.0, -1.5, 1.5)
    elif bkg_pdf == "mixed":
        pass
    elif bkg_pdf == "cmsshape":
        pass
    else:
        print("REQUESTED BACKGROUND SHAPE IS NOT SUPPORTED")
        sys.exit()

    # ------------------------------
    #  Histograms and PDFs building
    # ------------------------------
    mean, sigma, tau = [mean_f, mean_p], [sigma_f, sigma_p], [tau_f, tau_p]
    axis, histo_data, pdf_mc = [0, 0], [0, 0], [0, 0]
    smearing, conv_pdf, background = [0, 0], [0, 0], [0, 0]

    nexp, Nsig, Nbkg = [0, 0], [0, 0], [0, 0]
    sum_func, model, results = [0, 0], [0, 0], [0, 0]

    for cond in ["pass", "fail"]:  # PASS has index 1, FAIL has index 0

        idx = 1 if cond == "pass" else 0

        axis[idx] = ws.var(f"x_{cond}_({bin[0]}|{bin[1]})")
        axis[idx].setRange("fitRange", 50, 130)

        # binning = ROOT.RooUniformBinning(
        #     axis.getRange("fitRange")[0], axis.getRange("fitRange")[1], NBINS[idx])
        axis[idx].setBinning(
            ROOT.RooUniformBinning(axis[idx].getRange("fitRange")[0],
                                   axis[idx].getRange("fitRange")[1], NBINS[idx]), "cache")

        histo_data[idx] = ws.data(f"Minv_data_{cond}_({bin[0]}|{bin[1]})")

        pdf_mc[idx] = ROOT.RooHistPdf(
            f"pdf_mc_{cond}_({bin[0]}|{bin[1]})", "pdf MC", axis[idx],
            ws.data(f"Minv_mc_{cond}_({bin[0]}|{bin[1]})"))

        smearing[idx] = ROOT.RooGaussian(f"smearing_{cond}_({bin[0]}|{bin[1]})",
                                         "Gaussian smearing", axis[idx], mean[idx], sigma[idx])

        conv_pdf[idx] = ROOT.RooFFTConvPdf(f"conv_{cond}_({bin[0]}|{bin[1]})", f"Convolution pdf",
                                           axis[idx], pdf_mc[idx], smearing[idx], 3)
        conv_pdf[idx].setBufferFraction(bufFractions[idx])
        # conv_pdf[idx].setBufferStrategy(0)
        # conv_pdf[idx].setOperMode(3)

        if bkg_pdf == "expo":
            background[idx] = ROOT.RooExponential(f"expo_bkg_{cond}_({bin[0]},{bin[1]})",
                                                  "Exponential background", axis[idx], tau[idx])
        elif bkg_pdf == "mixed":
            pass
        elif bkg_pdf == "cmsshape":
            pass
        else:
            print("REQUESTED BACKGROUND SHAPE IS NOT SUPPORTED")
            sys.exit()

        # -----------------------
        #  Final models and fits
        # -----------------------

        n_events = histo_data[idx].sumEntries()

        nexp[idx] = ROOT.RooRealVar(
            f"nexp_{cond}_({bin[0]}|{bin[1]})", f"nexp {cond}",
            n_events, n_events-5*(n_events**0.5), n_events+5*(n_events**0.5))

        '''
        Nsig[idx] = ROOT.RooRealVar(
            f"nsig_{cond}_({bin[0]}|{bin[1]})", "#signal events",
            nexp[idx].getVal(), lims_num[idx][0]*nexp[idx].getVal(), nexp[idx].getMax())
        Nbkg[idx] = ROOT.RooRealVar(
            f"nbkg_{cond}_({bin[0]}|{bin[1]})", "#background events",
            lims_num[idx][1]*nexp[idx].getVal()/2., 0.0, lims_num[idx][1]*nexp[idx].getVal())
        '''

        Nsig[idx] = ROOT.RooRealVar(
            f"nsig_{cond}_({bin[0]}|{bin[1]})", "#signal events",
            0.9*nexp[idx].getVal(), 0.5, 1.5*nexp[idx].getVal())
        Nbkg[idx] = ROOT.RooRealVar(
            f"nbkg_{cond}_({bin[0]}|{bin[1]})", "#background events",
            0.1*nexp[idx].getVal(), 0.5, 1.5*nexp[idx].getVal())

        sum_func[idx] = ROOT.RooAddPdf(
                            f"sum_{cond}_({bin[0]}|{bin[1]})", "Signal+Bkg",
                            ROOT.RooArgList(conv_pdf[idx], background[idx]),
                            ROOT.RooArgList(Nsig[idx], Nbkg[idx]))

        # sum_func[idx].setNormRange("fitRange")

        model[idx] = ROOT.RooAddPdf(
            sum_func[idx], f'PDF_{cond}_({bin[0]}|{bin[1]})')
        # model[idx] = ROOT.RooExtendPdf("PDF_{cond}_({bin[0]}|{bin[1]}", "extend", conv_func[idx], nexp[idx])
        model[idx].setNormRange("fitRange")

        results[idx] = model[idx].fitTo(histo_data[idx],
                                        ROOT.RooCmdArg("Extended", 1),
                                        # ROOT.RooCmdArg("Range", 0, 0, 0, 0, 'fitRange'),
                                        ROOT.RooCmdArg("Minimizer", 0, 0, 0, 0,
                                                       "Minuit2", "migrad"),
                                        ROOT.RooCmdArg("Strategy", fit_strategy[idx]),
                                        # ROOT.RooCmdArg("MaxCalls", 100000),
                                        ROOT.RooCmdArg("Save", 1),
                                        ROOT.RooCmdArg("PrintLevel", verb)
                                        )
        results[idx].SetName(f"results_{cond}_({bin[0]}|{bin[1]})")

        chi2val, ndof = pearson_chi2_eval(histo_data[idx], model[idx], 
                                          histo_data[idx].numEntries(), results[idx])

        status = fit_quality_old(results[idx], (chi2val, ndof))

        low_nbkg = (Nbkg[idx].getVal() < 0.005*Nsig[idx].getVal())

        if (status is False) and refit_numbkg and low_nbkg and bkg_pdf=='expo':

            if cond == 'pass':
                tau_p.setVal(1)
                tau_p.setConstant()
            else:
                tau_f.setVal(1)
                tau_f.setConstant()
            Nbkg[idx].setVal(0)
            Nbkg[idx].setConstant()
            print("REFITTING WITHOUT BACKGROUND")
            print("\n\n\n\n")
            results[idx] = model[idx].fitTo(histo_data[idx],
                                            # ROOT.RooCmdArg("Extended", 1),
                                            ROOT.RooCmdArg("Range", 0, 0, 50., 130., 'fitRange'),
                                            ROOT.RooCmdArg("Minimizer", 0, 0, 0, 0,
                                                           "Minuit2", "migrad"),
                                            ROOT.RooCmdArg("Strategy", fit_strategy[idx]),
                                            # ROOT.RooCmdArg("MaxCalls", 100000),
                                            ROOT.RooCmdArg("Save", 1),
                                            ROOT.RooCmdArg("PrintLevel", verb)
                                            )
    
        results[idx].SetName(f"results_{cond}_({bin[0]}|{bin[1]})")
        
        if fit_quality_old(results[idx], (chi2val, ndof)) is False:
            results[idx].setStatus(5)

    res_pass, res_fail = results

    # ----------------------------
    #  Quality checks (and plots)
    # ----------------------------
    '''
    if fit_quality(res_pass) and fit_quality(res_fail):
        # ws.Import(model[0]), ws.Import(model[1])
        # ws.Import(results[0]), ws.Import(results[1])
        if figs:
            print(background[1].GetName())
            plot_distr_with_fit(
                axis[1], histo_data[1], model[1], bkg_pdf=background[1].GetName(),
                name=f"figs/fit_iso/pass_{bin[0]},{bin[1]}.png")
            plot_distr_with_fit(
                axis[0], histo_data[0], model[0], bkg_pdf=background[0].GetName(),
                name=f"figs/fit_iso/fail_{bin[0]},{bin[1]}.png")
    '''
    '''
    if test_bkg is True:
        null_bkg = llr_test_bkg(histo_data, model)
        present_bkg = not null_bkg
        print(f"Background is accepted? {present_bkg}")
    '''

    return res_pass, res_fail 


def fit_on_bin(type_eff, workspace, cond, bin, bkg_pdf, test_bkg=False,
               verb=-1, figs=False):
    """
    """

    if type(workspace.obj(f'PDF_{cond}_({bin[0]}|{bin[1]})')) is ROOT.RooAddPdf:

        print("Not possible to refit an existing PDF! \
              \nReturning the results obtained previously")
        return workspace.pdf(f"results_{cond}_({bin[0]}|{bin[1]})")

    elif type(workspace.obj(f'PDF_{cond}_({bin[0]}|{bin[1]})')) is ROOT.TObject:

        histo_data = workspace.data(f"Minv_data_{cond}_({bin[0]}|{bin[1]})")
        histo_mc = workspace.data(f"Minv_mc_{cond}_({bin[0]}|{bin[1]})")

        axis = ROOT.RooRealVar(workspace.var(f"x_{cond}_({bin[0]}|{bin[1]})"))

        axis.setRange("fitRange", 50, 130)

        NBINS = 6000
        binning = ROOT.RooUniformBinning(
            axis.getRange("fitRange")[0], axis.getRange("fitRange")[1], NBINS)
        axis.setBinning(binning, "cache")

        pdf_mc = ROOT.RooHistPdf(f"pdf_mc_{cond}_({bin[0]}|{bin[1]})",
                                 "pdf_mc", axis, histo_mc)
        # pdf_mc.setNormRange("fitRange")

        mean = ROOT.RooRealVar(
            f"mean_{cond}_({bin[0]}|{bin[1]})", "mean", 0, -5.0, 5.0)
        sigma = ROOT.RooRealVar(
            f"sigma_{cond}_({bin[0]}|{bin[1]})", "sigma", 2, 0.5, 6.0)

        smearing = ROOT.RooGaussian(f"smearing_{cond}_({bin[0]}|{bin[1]})",
                                    "Gaussian smearing", axis, mean, sigma)
        # smearing.setNormRange("fitRange")

        conv_func = ROOT.RooFFTConvPdf(f"conv_{cond}_({bin[0]}|{bin[1]})",
                                       f"Convolution {cond}", axis, pdf_mc, smearing, 3)
        conv_func.setBufferFraction(0.1)
        # conv_func.setBufferStrategy(0)
        # conv_func.setOperMode(3)
        # conv_func.setNormRange("fitRange")

        if (bkg_pdf == 'expo') or (bkg_pdf == 'mixed' and cond == 'pass'):
            tau = ROOT.RooRealVar(
                f"tau_{cond}_({bin[0]}|{bin[1]})", "tau", 0.0, -1.5, 1.5)
            background = ROOT.RooExponential(
                f"expo_bkg_{cond}_({bin[0]},{bin[1]})",
                "Exponential background", axis, tau)
            # background.setNormRange("fitRange")

        elif bkg_pdf == 'cmsshape' or (bkg_pdf == 'mixed' and cond == 'fail'):
            alpha = ROOT.RooRealVar(
                f"alpha_{cond}_({bin[0]}|{bin[1]})", "alpha", 60.0, 40.0, 130.0)
            beta = ROOT.RooRealVar(
                f"beta_{cond}_({bin[0]}|{bin[1]})", "beta", 5.0, 0.1, 40.0)
            gamma = ROOT.RooRealVar(
                f"gamma_{cond}_({bin[0]}|{bin[1]})", "gamma", 0.1, 0.0, 1.0)
            peak = ROOT.RooRealVar(
                f"peak_{cond}_({bin[0]}|{bin[1]})", "peak", 90.0)  # ,88.0, 92.0)

            background = ROOT.RooCMSShape(
                f"cmsshape_bkg_{cond}_({bin[0]}|{bin[1]})",
                "CMSShape background", axis, alpha, beta, gamma, peak)
            # background.setNormRange("fitRange")

        else:
            print("BKG shape given is not implemented! Retry with 'expo' or 'cmsshape'")
            sys.exit()

        events_data = histo_data.sumEntries()

        expected_num = ROOT.RooRealVar(
            "nexp", "nexp", events_data, events_data-10*(events_data**0.5), events_data+10*(events_data**0.5))

        Nsig = ROOT.RooRealVar(
            f"nsig_{cond}_({bin[0]}|{bin[1]})", "#signal events",
            events_data, 0.5*events_data, events_data + 5*ROOT.TMath.Sqrt(events_data))
        Nbkg = ROOT.RooRealVar(
            f"nbkg_{cond}_({bin[0]}|{bin[1]})", "#background events",
            0, 0.0, 0.05*events_data)

        sum_func = ROOT.RooAddPdf(f"sum_{cond}_({bin[0]}|{bin[1]})", "Signal+Bkg",
                                  ROOT.RooArgList(conv_func, background),
                                  ROOT.RooArgList(Nsig, Nbkg))
        sum_func.setNormRange("fitRange")

        model = ROOT.RooAddPdf(sum_func, f'PDF_{cond}_({bin[0]}|{bin[1]}')

        # model = ROOT.RooExtendPdf("extend", "extend", conv_func, expected_num)
        model.setNormRange("fitRange")
        '''
        expected_num = ROOT.RooRealVar(
            "nexp", "nexp", ROOT.RooArgSet(Nsig, Nbkg))
        '''
        # IN QUESTA VERSIONE DI ROOT IL METODO "fitTo" NON HA L'OPZIONE "MaxCalls"... SE SI VUOLE GESTIRE IL TUTTO IN MANIERA PIÙ TRASPARENTE SI PUÒ TRANQUILLAMENTE OPERARE TRAMITE UN OGGETTO "ROOMINIMIZER", PRENDENDO COME ESEMPIO IL SOURCEFILE DI FITTO
        res = model.fitTo(histo_data,
                          ROOT.RooCmdArg("Extended", 1),
                          # ROOT.RooCmdArg("Range", 0, 0, 0, 0, 'fitRange'),
                          # ROOT.RooCmdArg("ExternalConstraints", 0, 0, 0, 0, " ", " ", ROOT.RooArgSet(expected_num)),
                          ROOT.RooCmdArg("Minimizer", 0, 0, 0,
                                         0, "Minuit2", "migrad"),
                          ROOT.RooCmdArg("Strategy", 2),
                          # ROOT.RooCmdArg("MaxCalls", 100000),
                          ROOT.RooCmdArg("Save", 1),
                          ROOT.RooCmdArg("PrintLevel", verb)
                          )

        res.SetName(f"results_{cond}_({bin[0]}|{bin[1]})")

        # pearson_chi2_eval(histo_data, model, histo_data.numEntries(), res)

        if fit_quality(res) is True:
            workspace.Import(model)
            workspace.Import(res)

    else:
        print("******\nERROR in PDF types\n*******")
        sys.exit()

    '''
    if test_bkg is True:
        null_bkg = llr_test_bkg(histo_data, model)
        present_bkg = not null_bkg
        print(f"Background is accepted? {present_bkg}")
    '''

    if figs is True and fit_quality(res) is True:
        print(background.GetName())
        plot_distr_with_fit(axis, histo_data, model,
                            bkg_pdf=background,
                            pull=False,
                            name=f"figs/fit_iso/{cond}_{bin[0]},{bin[1]}.png")

    return res


def independent_efficiency(type_eff, bins,
                           bin_combinations=True,
                           bkg_strategy='expo', test_bkg=False,
                           verbose=0, figs=False):
    """
    """
    path = os.path.dirname(__file__)
    ROOT.gSystem.cd(path)

    file_ws = ROOT.TFile(
        f"root_files/ws/ws_{type_eff}_indep_{bkg_strategy}.root")
    ws = file_ws.Get("w")

    if bin_combinations is True:
        bins_pt = []
        bins_eta = []
        for i in range(len(bins[0])):
            for j in range(len(bins[1])):
                bins_pt.append(bins[0][i])
                bins_eta.append(bins[1][j])
        bins_list = [bins_pt, bins_eta]
    else:
        bins_list = bins

    bkg_pdf = bkg_strategy

    Nproblems = 0
    bins_with_problems = []

    for idx in range(len(bins_list[0])):
        bin_pt, bin_eta = bins_list[0][idx], bins_list[1][idx]

        res_pass = fit_on_bin(type_eff, ws, 'pass', (bin_pt, bin_eta), bkg_pdf,
                              test_bkg=test_bkg, verb=verbose, figs=figs)

        res_fail = fit_on_bin(type_eff, ws, 'fail', (bin_pt, bin_eta), bkg_pdf,
                              test_bkg=test_bkg, verb=verbose, figs=figs)

        # results.add_result(bin_pt, bin_eta, res_pass, res_fail, check=True)

        status = bool(fit_quality(res_pass)*fit_quality(res_fail))
        if status is False:
            print(f"\nBin {bin_pt},{bin_eta} has problems!\n")
            Nproblems += 1
            bins_with_problems.append(f"{bin_pt},{bin_eta}")
            res_pass.Print()
            #res_pass.correlationMatrix().Print()
            print("****")
            res_fail.Print("")
            ##res_fail.correlationMatrix().Print()
            print("****")
            print(res_pass.status(), res_fail.status())
            print(res_pass.covQual(), res_fail.covQual())
            print(res_pass.edm(), res_fail.edm())
            print(' ')

    print(f"NUM of problematic bins = {Nproblems}")
    # ws.Print()
    print(bins_with_problems)
    ws.writeToFile(f"root_files/{type_eff}_workspace.root")


if __name__ == '__main__':

    ROOT.gROOT.SetBatch(True)
    ROOT.PyConfig.IgnoreCommandLineOptions = True

    t0 = time.time()

    custom_pdfs = ['RooCBExGaussShape',
                   'RooDoubleCBFast', 'RooCMSShape', 'my_double_CB']
    import_pdf_library(custom_pdfs[2])

    type_eff = ("sa", "global", "ID", "iso", "trigger", "veto")
    t = type_eff[3]

    # results = results_manager("indep")

    '''
    bins_pt = [num for num in range(1, 2)]
    bins_eta = [num for num in range(1, 49)]
    '''

    '''
    # BIN NEI QUALI NON SI È RIUSCITO AD AVERE UN BUON FIT CON STRATEGIA "MIXED"
    bin_keys = ['1,1', '1,4', '1,11', '1,14', '1,17', '1,19', '1,20',
                '1,24', '1,27', '1,30', '1,31', '1,40', '1,41', '1,42', '1,45']
    '''

    bin_keys = ['1,16', '1,18', '1,21', '1,24',
                '1,29', '1,30', '1,32', '1,43', '1,47']

    bins_pt, bins_eta = [], []

    for key in bin_keys:
        bins_pt.append(key.split(',')[0])
        bins_eta.append(key.split(',')[1])

    bins = (bins_pt, bins_eta)
    independent_efficiency(t, bins, bin_combinations=False,
                           bkg_strategy='expo', verbose=-1, figs=True)

    t1 = time.time()

    print(f"TIME ELAPSED = {t1-t0}")

    # sys.exit()
    # sys.exit()
    # sys.exit()
    # sys.exit()
