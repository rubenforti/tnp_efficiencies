"""
"""

import ROOT
from utilities import makeAndSavePlot


def import_pdf_library(*functions):
    """
    """
    for function in functions:
        header_incl = ' #include "libCpp/'+function+'.h"'
        sourcefile = 'libCpp/'+function+'.cc'

        ctrl_head = ROOT.gInterpreter.Declare(header_incl)
        ctrl_source = ROOT.gSystem.CompileMacro(sourcefile, opt="ks")

        if not ctrl_head:
            print("ERROR in header loading")
            quit()
        if not ctrl_source == 1:
            print("ERROR in sourcefile compiling and loading")
            quit()


def import_histo(histo_th3, bin_pt, bin_eta):
    """
    Returns a RooDataHist of the variable TP_mass, in a (pt, eta) bin, selected
    from the TH3 given as input.
    """
    if len(bin_pt) == 1:
        bin_pt.append(bin_pt[0])
        bin_pt[0] -= 1
    if len(bin_eta) == 1:
        bin_eta.append(bin_eta[0])
        bin_eta[0] -= 1

    histo_th1 = histo_th3.ProjectionX(
        "histo_mass", bin_pt[0], bin_pt[1], bin_eta[0], bin_eta[1])
    xAxis = histo_th1.GetXaxis()

    x = ROOT.RooRealVar("x", "x", xAxis.GetXmin(),
                        xAxis.GetXmax(), unit="GeV/c^2")

    roohist = ROOT.RooDataHist("roohist", "roohist", [x], Import=histo_th1)
    
    n_events = histo_th1.GetEntries()
    print(n_events)

    return roohist, x, n_events


def fit_composite(histo, axis, signal, background, NEvents=10000):
    """
    """
    Nsig = ROOT.RooRealVar("nsig", "#signal events", 0, NEvents)
    Nbkg = ROOT.RooRealVar("nbkg", "#background events", 0, NEvents)
    sig_funcname = signal.GetName()
    bkg_funcname = background.GetName()

    model = ROOT.RooAddPdf("model", "model", [signal, background], [Nsig, Nbkg])

    model.fitTo(histo, Extended=True)
    makeAndSavePlot(axis, histo, model, name="provafit_composite.png")


if __name__ == '__main__':

    custom_pdfs = ['RooCBExGaussShape', 'RooDoubleCBFast', 'RooCMSShape', 'my_double_CB']

    import_pdf_library(custom_pdfs[2], custom_pdfs[3])

    f = ROOT.TFile("root_files/tnp_iso_data.root")
    histo_pass = f.pass_mu_RunGtoH
    histo_fail = f.fail_mu_RunGtoH

    dh, x, n_events = import_histo(histo_pass, [1], [1])
    # c = ROOT.TCanvas()
    # c.cd()

    # a = ROOT.RooRealVar("a", "a", 90, 80, 100)
    mean = ROOT.RooRealVar("b", "b", 91, 80, 100)
    sigma = ROOT.RooRealVar("c", "c", 2, 0.1, 10)
    a1 = ROOT.RooRealVar("d", "d", 1.5, 0, 3)
    n1 = ROOT.RooRealVar("e", "e", 1, -5, 5)
    a2 = ROOT.RooRealVar("f", "f", 1.5, 0, 3)
    n2 = ROOT.RooRealVar("g", "g", 1, 5, 5)


    '''
    tau = ROOT.RooRealVar("tau", "tau", -0.5, -1, 1)
    
    func1 = ROOT.RooCBExGaussShape("strange_CB", "strange_cb", x, b, c, d, e, f, g)
    func2 = ROOT.RooExponential("expo", "expo", x, tau)
    '''
    func3 = ROOT.RooVoigtian("CB", "CB", x, mean, a1, sigma) 
    func3.fitTo(dh)
    makeAndSavePlot(x, dh, func3, name="provafit_voigtian.png")

    # fit_composite(dh, x, func1, func2, n_events) 
    # makeAndSavePlot(x, dh, func1, func2, name="prova_custompdf.png")
