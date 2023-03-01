"""
"""

import ROOT
from utilities import makeAndSavePlot, init_CrystalBall, roodouble


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


def fit_composite(axis, histo, signal, background, NEvents=10000):
    """
    """
    Nsig = ROOT.RooRealVar("nsig", "#signal events", 0, NEvents)
    Nbkg = ROOT.RooRealVar("nbkg", "#background events", 0, NEvents)
    sig_funcname = signal.GetName()
    bkg_funcname = background.GetName()

    model = ROOT.RooAddPdf("model", "model", [signal, background], [Nsig, Nbkg])
    print(type(model))
    fitres = signal.fitTo(histo, Extended=True, Save=True)
    makeAndSavePlot(axis, histo, model, name="provafit_composite.png")


if __name__ == '__main__':

    custom_pdfs = ['RooCBExGaussShape', 'RooDoubleCBFast', 'RooCMSShape', 'my_double_CB']

    # import_pdf_library(custom_pdfs[2], custom_pdfs[3])

    f = ROOT.TFile("root_files/tnp_iso_data.root")
    histo_pass = f.pass_mu_RunGtoH
    histo_fail = f.fail_mu_RunGtoH

    dh, x, n_events = import_histo(histo_pass, [1], [1])
    
    c = ROOT.TCanvas()
    c.cd()
   
    # a = ROOT.RooRealVar("a", "a", 90, 80, 100)
    mean = ROOT.RooRealVar("b", "b", 91, 80, 100)
    sigma = ROOT.RooRealVar("c", "c", 2, 0.1, 10)
    a1 = ROOT.RooRealVar("d", "d", 1.5, 0, 3)
    n1 = ROOT.RooRealVar("e", "e", 1, -5, 5)
    a2 = ROOT.RooRealVar("f", "f", 1.5, 0, 3)
    n2 = ROOT.RooRealVar("g", "g", 1, 5, 5)
   
    cb_shape = ROOT.RooCBShape("cb", "cb", x, mean, sigma, a1, n1)
    
    tau = ROOT.RooRealVar("tau", "tau", -0.5, -1, 1)
    expo = ROOT.RooExponential("expo", "expo", x, tau)

    gamma = roodouble(1, 0.5, 3)
    voigt = ROOT.RooVoigtian("voigt", "voigt", x, mean, gamma, sigma)

    mean = ROOT.RooRealVar("mean", "mean", 91, 85, 97)
    sigma = ROOT.RooRealVar("sigma", "sigma", 0.1, 5)
    gauss = ROOT.RooGaussian("gauss", "gauss", x, mean, sigma)
    
    a2 = ROOT.RooRealVar("a2", "a2", -10, 10)
    a3 = ROOT.RooRealVar("a3", "a3", -10, 10)
    pol = ROOT.RooPolynomial("pol", "pol", x, [a2, a3])


    Nsig = ROOT.RooRealVar("nsig", "#signal events", 0, n_events)
    Nbkg = ROOT.RooRealVar("nbkg", "#background events", 0, n_events)
    sum_func = ROOT.RooAddPdf("sum", "sum", [gauss, pol], [Nsig, Nbkg])

    res = sum_func.fitTo(dh, Save=True)
    
    frame = x.frame("Titolo")
    sum_func.plotOn(frame, Components={expo}, LineStyle=':')
    dh.plotOn(frame)
    sum_func.plotOn(frame)
    frame.Draw()
    c.SaveAs("provafit_composite.png")


    # makeAndSavePlot(x, dh, sum_func, name="provafit_voigtian.png")
    
    # total_func = fit_composite(x, dh, cb_func, expo, n_events) 
    # makeAndSavePlot(x, dh, total_func, name="provafit_composite.png")
