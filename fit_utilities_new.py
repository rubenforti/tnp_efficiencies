"""
"""

import sys
import ROOT
from plot_functions import makeAndSavePlot
from utilities import import_Steve_histos
from stat_functions import pearson_chi2_eval, llr_test_bkg
from workspace_config import ws_init_std_pdf

'''
def fit_without_bkg(axis, t, histo_data, histo_mc, bins_pt_eta, events_data,
                    saveplot=False, verb=-1):
    """
    """

    if 'pass' in histo_data.GetTitle():
        id_flag = 'pass'
    elif 'fail' in histo_data.GetTitle():
        id_flag = 'fail'
    else:
        print("*******\nERROR in condition flag passed to fit function\n*******")

    pdf_mc = ROOT.RooHistPdf("pdf_mc", "pdf_mc", axis, histo_mc)
    mean = ROOT.RooRealVar("mean", "mean", 0, -2, 2)
    sigma = ROOT.RooRealVar("sigma", "sigma", 0.5, 0.001, 2)
    smearing = ROOT.RooGaussian("smearing", "smearing", axis, mean, sigma)

    axis.setBins(1000, "cache")
    conv_func = ROOT.RooFFTConvPdf("conv", "conv", axis, pdf_mc, smearing, 3)
    conv_func.setBufferFraction(0.1)

    res = conv_func.fitTo(histo_data, Extended=True,
                          Save=True, PrintLevel=verb)

    pearson_chi2_eval(histo_data, conv_func, histo_data.numEntries(), res)

    if saveplot is True:
        makeAndSavePlot(axis, histo_data, conv_func, pull=False,
                        name=f"figs/fit_{t}/{id_flag}_{bins_pt_eta[0]}_{bins_pt_eta[1]}.pdf")

    return res
'''


def fit_with_bkg(type_eff, cond, bin, test_bkg=False, verb=-1,):
    """
    """
    file_ws = ROOT.TFile(f"root_files/{type_eff}_workspace.root")
    workspace = file_ws.Get("w")

    ws_init_std_pdf(workspace, cond, bin)

    axis = workspace[f"x_{cond}_({bin[0]},{bin[1]})"]
    # axis.SetTitle(f"TP M_inv  {cond} ({bin[0]}, {bin[1]})")
    print(type(axis))

    '''
    if not ((type(histo_data) is ROOT.RooDataHist) and (type(histo_data) is ROOT.RooDataHist)):
        print("*******\nERROR in loading the RooDataHist objects of MC data\n*******")
        sys.exit()
    '''

    print(f"IMPORTING: {cond}, {bin[0]}, {bin[1]}")
    print(type(workspace[f'PDF_{cond}_({bin[0]},{bin[1]})']))

    if type(workspace[f'PDF_{cond}_({bin[0]},{bin[1]})']) is ROOT.RooAddPdf:

        histo_data = workspace[f"Minv_data_{cond}_({bin[0]},{bin[1]})"]

        model = ROOT.RooAddPdf(workspace[f'PDF_{cond}_({bin[0]},{bin[1]})'],
                               f'PDF_{cond}_({bin[0]},{bin[1]})_copy')

        axis.setBins(10000, "cache")

    elif type(workspace[f'PDF_{cond}_({bin[0]},{bin[1]})']) is ROOT.TObject:

        '''
        histo_data = ROOT.RooDataHist(workspace[f"Minv_data_{cond}_({bin[0]},{bin[1]})"],
                                      f"Minv_data_{cond}_({bin[0]},{bin[1]})_copy")
        histo_mc = ROOT.RooDataHist(workspace[f"Minv_mc_{cond}_({bin[0]},{bin[1]})"],
                                    f"Minv_mc_{cond}_({bin[0]},{bin[1]})_copy")
        '''
        histo_data = workspace[f"Minv_data_{cond}_({bin[0]},{bin[1]})"]
        histo_mc = workspace[f"Minv_mc_{cond}_({bin[0]},{bin[1]})"]
        pdf_mc = ROOT.RooHistPdf(
            f"pdf_mc_{cond}_({bin[0]}_{bin[1]})", "pdf_mc", axis, histo_mc)

        smearing = workspace[f"gaus_smearing_{cond}_({bin[0]},{bin[1]})"]
        background = workspace[f"expo_bkg_{cond}_({bin[0]},{bin[1]})"]

        axis.setBins(10000, "cache")
        conv_func = ROOT.RooFFTConvPdf(f"conv_{cond}_({bin[0]}_{bin[1]})", "Convolution",
                                       axis, pdf_mc, smearing, 3)
        conv_func.setBufferFraction(0.5)

        events_data = histo_data.sumEntries()

        Nsig = ROOT.RooRealVar(
            "nsig", "#signal events",
            events_data, 0, events_data + 4*ROOT.TMath.Sqrt(events_data))
        Nbkg = ROOT.RooRealVar(
            "nbkg", "#background events",
            0.005*events_data, -0., 0.5*events_data)

        sum_func = ROOT.RooAddPdf(f"Sum_({bin[0]}_{bin[1]})", "Signal+Bkg",
                                  [conv_func, background], [Nsig, Nbkg])

        model = ROOT.RooAddPdf(sum_func)
        model.SetName(f'PDF_{cond}_({bin[0]},{bin[1]})')

    else:
        print("******\nERROR in PDF types\n*******")
        sys.exit()

    res = model.fitTo(histo_data, Extended=True, Save=True, PrintLevel=verb)

    pearson_chi2_eval(histo_data, model, histo_data.numEntries(), res)

    if res.status() == 0 and res.covQual() == 3 and res.edm() < 1e-2:
        workspace.Import(model, RenameConflictNodes="NEW")
        # axis.SetName(f"x_{cond}_({bin[0]}, {bin[1]})")
        # workspace.Import(axis)
        workspace.writeToFile(f"root_files/{type_eff}_workspace.root")

    if test_bkg is True:
        null_bkg = llr_test_bkg(histo_data, model)
        present_bkg = not null_bkg
        print(f"Background is accepted? {present_bkg}")

    '''
    if verb != -1:
        makeAndSavePlot(axis, histo_data, model, pull=False,
                        name=f"figs/fit_{t}/{flag}_{bin[0]}_{bin[1]}.pdf")
    '''

    return model, res
