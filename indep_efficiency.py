"""
"""

import ROOT
import time
import os
import sys
from results_utils import res_manager_indep
from plot_functions import makeAndSavePlot
from utilities import import_pdf_library, fit_quality, pearson_chi2_eval, llr_test_bkg


def fit_on_bin(type_eff, workspace, cond, bin, bkg_pdf, test_bkg=False,
               verb=-1, figs=False):
    """
    """

    if type(workspace[f'PDF_{cond}_({bin[0]}|{bin[1]})']) is ROOT.RooAddPdf:

        print("Not possible to refit an existing PDF! \
              \nReturning the results obtained previously")
        return workspace[f"results_{cond}_({bin[0]}|{bin[1]})"]

    elif type(workspace[f'PDF_{cond}_({bin[0]}|{bin[1]})']) is ROOT.TObject:

        histo_data = workspace[f"Minv_data_{cond}_({bin[0]}|{bin[1]})"]
        histo_mc = workspace[f"Minv_mc_{cond}_({bin[0]}|{bin[1]})"]

        axis = ROOT.RooRealVar(workspace[f"x_{cond}_({bin[0]}|{bin[1]})"])

        axis.setRange("fitRange", 50, 130)

        NBINS = 8000
        binning = ROOT.RooUniformBinning(axis.getRange(
            "fitRange")[0], axis.getRange("fitRange")[1], NBINS)
        axis.setBinning(binning, "cache")

        pdf_mc = ROOT.RooHistPdf(f"pdf_mc_{cond}_({bin[0]}|{bin[1]})",
                                 "pdf_mc", axis, histo_mc)

        mean = ROOT.RooRealVar(
            f"mean_{cond}_({bin[0]}|{bin[1]})", "mean", 0, -2, 2)
        sigma = ROOT.RooRealVar(
            f"sigma_{cond}_({bin[0]}|{bin[1]})", "sigma", 2, 0.1, 5)

        smearing = ROOT.RooGaussian(f"smearing_{cond}_({bin[0]}|{bin[1]})",
                                    "Gaussian smearing", axis, mean, sigma)

        if (bkg_pdf == 'expo') or (bkg_pdf == 'mixed' and cond == 'pass'):
            tau = ROOT.RooRealVar(
                f"tau_{cond}_({bin[0]}|{bin[1]})", "tau", -0.1, -3, 0.0001)
            background = ROOT.RooExponential(
                f"expo_bkg_{cond}_({bin[0]},{bin[1]})",
                "Exponential background", axis, tau)

        elif bkg_pdf == 'cmsshape' or (bkg_pdf == 'mixed' and cond == 'fail'):
            alpha = ROOT.RooRealVar(
                f"alpha_{cond}_({bin[0]}|{bin[1]})", "alpha", 60.0, 40.0, 130.0)
            beta = ROOT.RooRealVar(
                f"beta_{cond}_({bin[0]}|{bin[1]})", "beta", 2.5, 0.01, 15.0)
            gamma = ROOT.RooRealVar(
                f"gamma_{cond}_({bin[0]}|{bin[1]})", "gamma", 0.1, 0.0001, 0.2)
            peak = ROOT.RooRealVar(
                f"peak_{cond}_({bin[0]}|{bin[1]})", "peak", 90.0)  # ,88.0, 92.0)

            background = ROOT.RooCMSShape(
                f"cmsshape_bkg_{cond}_({bin[0]}|{bin[1]})",
                "CMSShape background", axis, alpha, beta, gamma, peak)

        else:
            print("BKG shape given is not implemented! Retry with 'expo' or 'cmsshape'")
            sys.exit()

        conv_func = ROOT.RooFFTConvPdf(f"conv_{cond}_({bin[0]}|{bin[1]})",
                                       f"Convolution {cond}", axis, pdf_mc, smearing, 3)
        conv_func.setBufferFraction(0.2)
        conv_func.setBufferStrategy(2)
        # conv_func.prepareFFTBinning(axis)

        events_data = histo_data.sumEntries()

        Nsig = ROOT.RooRealVar(
            f"nsig_{cond}_({bin[0]}|{bin[1]})", "#signal events",
            events_data, 0.5*events_data, events_data + 4*ROOT.TMath.Sqrt(events_data))
        Nbkg = ROOT.RooRealVar(
            f"nbkg_{cond}_({bin[0]}|{bin[1]})", "#background events",
            0.001*events_data, 0.0, 0.2*events_data)

        sum_func = ROOT.RooAddPdf(f"sum_{cond}_({bin[0]}|{bin[1]})", "Signal+Bkg",
                                  [conv_func, background], [Nsig, Nbkg])

        model = ROOT.RooAddPdf(sum_func, f'PDF_{cond}_({bin[0]}|{bin[1]}')

        '''
        res0 = model.fitTo(histo_data,
                           Extended=True,
                           # Range='fitRange',
                           # ExternalConstraints=ROOT.RooArgSet("nexp"),
                           # BatchMode=False,
                           # Offset=True,
                           Strategy=0,
                           MaxCalls=1000,
                           Save=True,
                           PrintLevel=-1)
        '''
        res = model.fitTo(histo_data,
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

    else:
        print("******\nERROR in PDF types\n*******")
        sys.exit()

    if test_bkg is True:
        null_bkg = llr_test_bkg(histo_data, model)
        present_bkg = not null_bkg
        print(f"Background is accepted? {present_bkg}")

    if figs is True:
        makeAndSavePlot(axis, histo_data, model,
                        bkg_name=background.GetName(), pull=False,
                        name=f"{cond}_{bin[0]}_{bin[1]}.png")

    return res


def independent_efficiency(type_eff, bins, results, bin_combinations=True,
                           test_bkg=False, verbose=0):
    """
    """

    path = os.path.dirname(__file__)
    ROOT.gSystem.cd(path)

    ROOT.Math.MinimizerOptions.SetDefaultTolerance(2e-2)
    ROOT.Math.MinimizerOptions.SetDefaultMaxIterations(100000)
    ROOT.Math.MinimizerOptions.SetDefaultErrorDef(0.5)
    ROOT.Math.MinimizerOptions.SetDefaultStrategy(2)

    file_ws = ROOT.TFile(f"root_files/ws/{type_eff}_workspace_indep.root")
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

    bkg_pdf = 'mixed'

    Nproblems = 0

    for idx in range(len(bins_list[0])):
        bin_pt, bin_eta = bins_list[0][idx], bins_list[1][idx]
        res_pass = fit_on_bin(type_eff, ws, 'pass', (bin_pt, bin_eta), bkg_pdf,
                              test_bkg=test_bkg, verb=verbose)
        res_fail = fit_on_bin(type_eff, ws, 'fail', (bin_pt, bin_eta), bkg_pdf,
                              test_bkg=test_bkg, verb=verbose)

        # results.add_result(res_pass, res_fail, bin_pt, bin_eta, check=True)
        status = bool(fit_quality(res_pass)*fit_quality(res_fail))

        if status is False or 1.0 > 0.0:
            print(f"\nBin {bin_pt},{bin_eta} has problems!\n")
            Nproblems += 1
            res_pass.Print()
            res_pass.correlationMatrix().Print()
            print("****")
            res_fail.Print("")
            res_fail.correlationMatrix().Print()
            print("****")
            print(res_pass.status(), res_fail.status())
            print(res_pass.covQual(), res_fail.covQual())
            print(res_pass.edm(), res_fail.edm())
            print(' ')

    print(f"NUM of problematic bins = {Nproblems}")
    # ws.Print()
    ws.writeToFile(f"root_files/{type_eff}_workspace.root")


if __name__ == '__main__':

    t0 = time.time()
    custom_pdfs = ['RooCBExGaussShape',
                   'RooDoubleCBFast', 'RooCMSShape', 'my_double_CB']
    import_pdf_library(custom_pdfs[2])

    type_eff = ("sa", "global", "ID", "iso", "trigger", "veto")
    t = type_eff[3]

    results = res_manager_indep()

    bins_pt = [num for num in range(1, 2)]
    bins_eta = [num for num in range(1, 2)]

    results.open("indep_eff_results.pkl")

    '''
    bin_keys = results.get_problematic_bins()

    bins_pt, bins_eta = [], []

    for key in bin_keys:
        bins_pt.append(key.split(',')[0])
        bins_eta.append(key.split(',')[1])

    print(bins_pt)
    print(bins_eta)
    '''

    bins = (bins_pt, bins_eta)
    independent_efficiency(t, bins, results, bin_combinations=True, verbose=1)

    '''
    file_ws = ROOT.TFile(f"root_files/{t}_workspace.root")
    ws = file_ws.Get("w")
    ws.Print()
    '''
    # results.view_efficiencies()

    #results.write("indep_eff_results.pkl")
    #print("RISULTATI SCRITTI SU PICKLE FILE")

    t1 = time.time()

    print(f"TEMPO = {t1-t0}")
