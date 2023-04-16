"""
"""

import ROOT
import time
import os
import sys
from results_utilities import res_manager_indep
from stat_functions import pearson_chi2_eval, llr_test_bkg
from plot_functions import makeAndSavePlot
from utilities import import_pdf_library
from array import array


def fit_on_bin(type_eff, workspace, cond, bin, bkg, test_bkg=False, verb=-1):
    """
    """

    histo_data = workspace[f"Minv_data_{cond}_({bin[0]},{bin[1]})"]
    histo_mc = workspace[f"Minv_mc_{cond}_({bin[0]},{bin[1]})"]

    axis = workspace[f"x_{cond}_({bin[0]},{bin[1]})"]

    if type(workspace[f'PDF_{cond}_({bin[0]},{bin[1]})']) is ROOT.RooAddPdf:

        '''
        print("Fit with existing PDF not implemented yet!")
        sys.exit()
<<<<<<< Updated upstream
=======
        '''

        init_pdf = workspace[f'PDF_{cond}_({bin[0]},{bin[1]})']

        model = ROOT.RooAddPdf(init_pdf, f"{init_pdf.GetName()}_copy")

        res = model.fitTo(histo_data, Extended=True,
                          Save=True, PrintLevel=verb)

        print("PARAMETERS AFTER FIT")
        pars = model.getParameters(histo_data)
        pars.Print("v")

        pearson_chi2_eval(histo_data, model, histo_data.numEntries(), res)

        if res.status() == 0 and res.covQual() == 3 and res.edm() < 1e-4:
            workspace.Import(model, RenameAllNodes='NEW',
                             RenameAllVariables='NEW')
            workspace.writeToFile(f"root_files/{type_eff}_workspace.root")
>>>>>>> Stashed changes

    elif type(workspace[f'PDF_{cond}_({bin[0]},{bin[1]})']) is ROOT.TObject:

        pdf_mc = ROOT.RooHistPdf(f"pdf_mc_{cond}_({bin[0]},{bin[1]})",
                                 "pdf_mc", axis, histo_mc)

        mean = ROOT.RooRealVar(
            f"mean_{cond}_({bin[0]},{bin[1]})", "mean", 0, -2, 2)
        sigma = ROOT.RooRealVar(
            f"sigma_{cond}_({bin[0]},{bin[1]})", "sigma", 0.5, 0.001, 2)
        smearing = ROOT.RooGaussian(f"gaus_smearing_{cond}_({bin[0]},{bin[1]})",
                                    "gaussian smearing", axis, mean, sigma)

        if bkg == 'expo':
            tau = ROOT.RooRealVar(
                f"tau_{cond}_({bin[0]},{bin[1]})", "tau", -10, 0)
            background = ROOT.RooExponential(
                f"expo_bkg_{cond}_({bin[0]},{bin[1]})",
                "exponential background", axis, tau)

        elif bkg == 'cmsshape':
            alpha = ROOT.RooRealVar(
                f"alpha_{cond}_({bin[0]},{bin[1]})", "alpha", 45, 10, 60)
            beta = ROOT.RooRealVar(
                f"beta_{cond}_({bin[0]},{bin[1]})", "beta", 6, 0, 10)
            gamma = ROOT.RooRealVar(
                f"gamma_{cond}_({bin[0]},{bin[1]})", "gamma", 2, 0, 20)
            peak = ROOT.RooRealVar(
                f"peak_{cond}_({bin[0]},{bin[1]})", "peak", 10, 0, 25)

            background = ROOT.RooCMSShape(
                f"cmsshape_bkg_{cond}_({bin[0]},{bin[1]})",
                "CMSShape background", axis, alpha, beta, gamma, peak)

        else:
            print("BKG shape given is not implemented! Retry with 'expo' or 'cmsshape'")
            sys.exit()

        axis.setBins(10000, "cache")
        conv_func = ROOT.RooFFTConvPdf(f"conv_{cond}_({bin[0]}_{bin[1]})",
                                       "Convolution", axis, pdf_mc, smearing, 3)
        conv_func.setBufferFraction(0.5)

        events_data = histo_data.sumEntries()

        Nsig = ROOT.RooRealVar(
            f"nsig_{cond}_({bin[0]},{bin[1]})", "#signal events",
            events_data, 0, events_data + 4*ROOT.TMath.Sqrt(events_data))
        Nbkg = ROOT.RooRealVar(
            f"nbkg_{cond}_({bin[0]},{bin[1]})", "#background events",
            0.005*events_data, -0., 0.5*events_data)

        sum_func = ROOT.RooAddPdf(f"sum_{cond}_({bin[0]}_{bin[1]})", "Signal+Bkg",
                                  [conv_func, background], [Nsig, Nbkg])

        model = ROOT.RooAddPdf(sum_func, f'PDF_{cond}_({bin[0]},{bin[1]})')

        res = model.fitTo(histo_data, Extended=True,
                          Save=True, PrintLevel=verb)

        pearson_chi2_eval(histo_data, model, histo_data.numEntries(), res)

        if res.status() == 0 and res.covQual() == 3 and res.edm() < 1e-4:
            workspace.Import(model)
            # sum_func.SetName(f"PDF_{cond}_({bin[0]}_{bin[1]})_init")
            # workspace.Import(sum_func)

    else:
        print("******\nERROR in PDF types\n*******")
        sys.exit()

    if test_bkg is True:
        null_bkg = llr_test_bkg(histo_data, model)
        present_bkg = not null_bkg
        print(f"Background is accepted? {present_bkg}")

    if verb != -1:
        makeAndSavePlot(axis, histo_data, model, pull=False,
                        name=f"{cond}_{bin[0]}_{bin[1]}.pdf")

    return res


def independent_efficiency(type_eff, bins_pt, bins_eta, results, background,
                           test_bkg=False, verbose=-1):
    """
    """

    path = os.path.dirname(__file__)
    ROOT.gSystem.cd(path)

    file_ws = ROOT.TFile(f"root_files/{type_eff}_workspace.root")
    ws = file_ws.Get("w")

    for bin_pt in bins_pt:
        for bin_eta in bins_eta:
            res_pass = fit_on_bin(type_eff, ws, 'pass', (bin_pt, bin_eta),
                                  background, test_bkg=test_bkg, verb=verbose)
            res_fail = fit_on_bin(type_eff, ws, 'fail', (bin_pt, bin_eta),
                                  background, test_bkg=test_bkg, verb=verbose)

            results.add_result(res_pass, res_fail, bin_pt, bin_eta)
            status = results.check_fit_status(
                bin_pt, bin_eta, conditions='all')

            print(" ")
            if status != 0:
                print(f"Bin {bin_pt},{bin_eta} has {status} problems!")
            if verbose != 0 or status != 0 or status == 0:
                print(res_pass.status(), res_fail.status())
                print(res_pass.covQual(), res_fail.covQual())
                print(res_pass.edm(), res_fail.edm())
                print(' ')

    #ws.Print()
    ws.writeToFile(f"root_files/{type_eff}_workspace.root")


if __name__ == '__main__':

    t0 = time.time()
    custom_pdfs = ['RooCBExGaussShape',
                   'RooDoubleCBFast', 'RooCMSShape', 'my_double_CB']
    import_pdf_library(custom_pdfs[2])

    type_eff = ("sa", "global", "ID", "iso", "trigger", "veto")
    t = type_eff[3]

    results = res_manager_indep()

    bins_pt = array('I', [num for num in range(1, 2)])
<<<<<<< Updated upstream
    bins_eta = array('I', [num for num in range(1, 2)])

    independent_efficiency(t, bins_pt, bins_eta,
                           'cmsshape', results, verbose=0)
=======
    bins_eta = array('I', [num for num in range(4, 5)])

    independent_efficiency(t, bins_pt, bins_eta, results, verbose=0)
>>>>>>> Stashed changes

    '''
    file_ws = ROOT.TFile(f"root_files/{t}_workspace.root")
    ws = file_ws.Get("w")
    ws.Print()
    '''
    print(" ")
    probs = results.get_problematic_bins()
    print(f'NUM PROBLEMI = {len(probs)}')
    print(" ")

    '''
    results.write("indep_eff_results.pkl")
    print("RISULTATI SCRITTI SU PICKLE FILE")
    '''
    t1 = time.time()

    print(f"TEMPO = {t1-t0}")
