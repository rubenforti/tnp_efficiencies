"""
"""

import sys
import ROOT
from plot_functions import makeAndSavePlot
from utilities import import_Steve_histos
from stat_functions import pearson_chi2_eval, llr_test_bkg
from workspace_config import ws_init_std_pdf


def fit_on_bin(type_eff, workspace, cond, bin, test_bkg=False, verb=-1):
    """
    """

    histo_data = workspace[f"Minv_data_{cond}_({bin[0]},{bin[1]})"]
    histo_mc = workspace[f"Minv_mc_{cond}_({bin[0]},{bin[1]})"]

    axis = workspace[f"x_{cond}_({bin[0]},{bin[1]})"]

    if type(workspace[f'PDF_{cond}_({bin[0]},{bin[1]})']) is ROOT.RooAddPdf:

        print("Fit with existing PDF not implemented yet!")
        sys.exit()

        '''
        init_pdf = workspace[f'PDF_{cond}_({bin[0]},{bin[1]})']

        model = ROOT.RooAddPdf(init_pdf, f"{init_pdf.GetName()}_NEW")

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
        '''

    elif type(workspace[f'PDF_{cond}_({bin[0]},{bin[1]})']) is ROOT.TObject:

        pdf_mc = ROOT.RooHistPdf(f"pdf_mc_{cond}_({bin[0]},{bin[1]})",
                                 "pdf_mc", axis, histo_mc)

        # ws_init_std_pdf(workspace, cond, bin)

        mean = ROOT.RooRealVar(
            f"mean_{cond}_({bin[0]},{bin[1]})", "mean", 0, -2, 2)
        sigma = ROOT.RooRealVar(
            f"sigma_{cond}_({bin[0]},{bin[1]})", "sigma", 0.5, 0.001, 2)
        smearing = ROOT.RooGaussian(f"gaus_smearing_{cond}_({bin[0]},{bin[1]})",
                                    "gaussian smearing", axis, mean, sigma)

        tau = ROOT.RooRealVar(f"tau_{cond}_({bin[0]},{bin[1]})", "tau", -10, 0)
        background = ROOT.RooExponential(f"expo_bkg_{cond}_({bin[0]},{bin[1]})",
                                         "exponential background", axis, tau)

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

    else:
        print("******\nERROR in PDF types\n*******")
        sys.exit()

    if test_bkg is True:
        null_bkg = llr_test_bkg(histo_data, model)
        present_bkg = not null_bkg
        print(f"Background is accepted? {present_bkg}")

    '''
    if verb != -1:
        makeAndSavePlot(axis, histo_data, model, pull=False,
                        name=f"figs/fit_{t}/{flag}_{bin[0]}_{bin[1]}.pdf")
    '''

    return res


if __name__ == '__main__':

    file_ws = ROOT.TFile(f"root_files/iso_workspace.root")
    workspace = file_ws.Get("w")

    histo_data = workspace['Minv_data_pass_(1,1)']

    model_new = workspace['PDF_pass_(1,1)_NEW_NEW']
    model = workspace['PDF_pass_(1,1)']

    print("PARAMETERS OLD")
    print(workspace['nsig_pass_(1,1)'])
    print(" ")
    print("PARAMETERS NEW")
    pars_new = model_new.getParameters(histo_data)
    pars_new.Print("v")
