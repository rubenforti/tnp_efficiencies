"""
"""

import ROOT
import time
import os
import sys
from utilities.results_utils import results_manager
from utilities.plot_utils import plot_distr_with_fit, plot_pass_and_fail
from utilities.general_utils import import_pdf_library, fit_quality, check_chi2


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
        return (res_pass, res_fail)
    else:
        return (0, 0)


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
        f"sigma_pass_({bin[0]}|{bin[1]})", "sigma", 0.5, 0.1, 5.0)
    mean_f = ROOT.RooRealVar(
        f"mean_fail_({bin[0]}|{bin[1]})", "mean", 0, -5.0, 5.0)
    sigma_f = ROOT.RooRealVar(
        f"sigma_fail_({bin[0]}|{bin[1]})", "sigma", 0.5, 0.1, 5.0)

    if bkg_pdf == "expo":
        tau_p = ROOT.RooRealVar(
            f"tau_pass_({bin[0]}|{bin[1]})", "tau", 0.0, -5, 5)
        tau_f = ROOT.RooRealVar(
            f"tau_fail_({bin[0]}|{bin[1]})", "tau", 0.0, -5, 5)
        tau = [tau_f, tau_p]
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
    mean, sigma = [mean_f, mean_p], [sigma_f, sigma_p]
    axis, histo_data, pdf_mc = [0, 0], [0, 0], [0, 0]
    smearing, conv_pdf, background = [0, 0], [0, 0], [0, 0]

    nexp, Nsig, Nbkg = [0, 0], [0, 0], [0, 0]
    sum_func, model, results, fit_status = [0, 0], [0, 0], [0, 0], [0, 0]

    for cond in ["pass", "fail"]:  # PASS has index 1, FAIL has index 0

        idx = 1 if cond == "pass" else 0

        axis[idx] = ws.var(f"x_{cond}_({bin[0]}|{bin[1]})")
        axis[idx].setRange("fitRange", 60, 120)

        # binning = ROOT.RooUniformBinning(
        #     axis.getRange("fitRange")[0], axis.getRange("fitRange")[1], NBINS[idx])
        axis[idx].setBinning(
            ROOT.RooUniformBinning(axis[idx].getRange("fitRange")[0],
                                   axis[idx].getRange("fitRange")[1], NBINS[idx]), "cache")

        histo_data[idx] = ws.data(f"Minv_data_{cond}_({bin[0]}|{bin[1]})")

        pdf_mc[idx] = ROOT.RooHistPdf(
            f"pdf_mc_{cond}_({bin[0]}|{bin[1]})", "pdf MC", axis[idx],
            ws.data(f"Minv_mc_{cond}_({bin[0]}|{bin[1]})"), 3)

        smearing[idx] = ROOT.RooGaussian(f"smearing_{cond}_({bin[0]}|{bin[1]})",
                                         "Gaussian smearing", axis[idx], mean[idx], sigma[idx])

        conv_pdf[idx] = ROOT.RooFFTConvPdf(f"conv_{cond}_({bin[0]}|{bin[1]})", f"Convolution pdf",
                                           axis[idx], pdf_mc[idx], smearing[idx])
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
 
        Nsig[idx] = ROOT.RooRealVar(f"nsig_{cond}_({bin[0]}|{bin[1]})", "#signal events",
                                    0.9*nexp[idx].getVal(), 0.5, 1.5*nexp[idx].getVal())
        Nbkg[idx] = ROOT.RooRealVar(f"nbkg_{cond}_({bin[0]}|{bin[1]})", "#background events",
                                    0.1*nexp[idx].getVal(), 0.5, 1.5*nexp[idx].getVal())           

        sum_func[idx] = ROOT.RooAddPdf(
                            f"sum_{cond}_({bin[0]}|{bin[1]})", "Signal+Bkg",
                            ROOT.RooArgList(conv_pdf[idx], background[idx]),
                            ROOT.RooArgList(Nsig[idx], Nbkg[idx]))

        # sum_func[idx].setNormRange("fitRange")

        model[idx] = ROOT.RooAddPdf(
            sum_func[idx], f'PDF_{cond}_({bin[0]}|{bin[1]})')

        model[idx].setNormRange("fitRange")

        results[idx] = model[idx].fitTo(histo_data[idx],
                                        # ROOT.RooFit.Extended(1),
                                        ROOT.RooFit.Range("fitRange"),
                                        ROOT.RooFit.Minimizer("Minuit2"),
                                        ROOT.RooFit.Strategy(fit_strategy[idx]),
                                        # ROOT.RooFit.MaxCalls(100000),
                                        ROOT.RooFit.Save(1),
                                        ROOT.RooFit.PrintLevel(verb)
                                        )

        status_chi2 = check_chi2(histo_data[idx], model[idx], results[idx])
        status = bool(status_chi2*fit_quality(results[idx], old_checks=True))

        print(status_chi2)
        print(status)

        low_nbkg = (Nbkg[idx].getVal() < 0.005*Nsig[idx].getVal())
        print(Nbkg[idx].getVal())
        print(Nsig[idx].getVal())
        print(low_nbkg)

        if (status is False) and refit_numbkg and low_nbkg and bkg_pdf=='expo':
            results[idx].Print()
            tau[idx].setVal(1)
            tau[idx].setConstant()
            Nbkg[idx].setVal(0)
            Nbkg[idx].setConstant()
            print("REFITTING WITHOUT BACKGROUND")
            print("\n\n\n\n")
            results[idx] = model[idx].fitTo(histo_data[idx],
                                            # ROOT.RooFit.Extended(1),
                                            ROOT.RooFit.Range("fitRange"),
                                            ROOT.RooFit.Minimizer("Minuit2"),
                                            ROOT.RooFit.Strategy(fit_strategy[idx]),
                                            # ROOT.RooFit.MaxCalls(100000),
                                            ROOT.RooFit.Save(1),
                                            ROOT.RooFit.PrintLevel(verb)
                                            )
    
        if check_chi2(histo_data[idx], model[idx], results[idx]) is False:
            results[idx].SetTitle("Chi2_not_passed")

        fit_status[idx] = fit_quality(results[idx], old_checks=True)

        results[idx].SetName(f"results_{cond}_({bin[0]}|{bin[1]})")

    # ----------------------------------
    #  Quality checks (and plots) - NEW
    # ----------------------------------
    if bool(fit_status[0]*fit_status[1]) is True:
        ws.Import(model[0]), ws.Import(model[1])
        ws.Import(results[0]), ws.Import(results[1])
        if figs:
            bkg_names =[background[0].GetName(), background[1].GetName()]
            plot_pass_and_fail(axis, histo_data, model, bkg_names, 
                               name=f"figs/fit_iso/bin_{bin[0]},{bin[1]}.pdf")
    else:
        bkg_names =[background[0].GetName(), background[1].GetName()]
        plot_pass_and_fail(axis, histo_data, model, bkg_names, 
                            name=f"figs/check_fits/bin_{bin[0]},{bin[1]}.pdf")
    
    res_pass, res_fail = results

    '''
    if test_bkg is True:
        null_bkg = llr_test_bkg(histo_data, model)
        present_bkg = not null_bkg
        print(f"Background is accepted? {present_bkg}")
    '''

    return res_pass, res_fail

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

    bin_keys = ['6,7', '6,17', '6,23', '9,22', '10,21']




    bins_pt, bins_eta = [], []

    for key in bin_keys:
        bins_pt.append(key.split(',')[0])
        bins_eta.append(key.split(',')[1])

    bins = (bins_pt, bins_eta)
    independent_efficiency(t, bins, bin_combinations=False,
                           bkg_strategy='expo', verbose=-1, figs=True)

    t1 = time.time()

    print(f"TIME ELAPSED = {t1-t0}")
