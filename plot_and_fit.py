"""
"""

import ROOT
from utilities import import_pdf_library, import_histo, makeAndSavePlot


# Needs to be improved !!!
#def fit_composite(axis, histo, signal, background, NEvents=10000):
#
#    Nsig = ROOT.RooRealVar("nsig", "#signal events", 0, NEvents)
#    Nbkg = ROOT.RooRealVar("nbkg", "#background events", 0, NEvents)
#    sig_funcname = signal.GetName()
#    bkg_funcname = background.GetName()
#
#    model = ROOT.RooAddPdf(
#        "model", "model", [signal, background], [Nsig, Nbkg])
#    print(type(model))
#    fitres = signal.fitTo(histo, Extended=True, Save=True)
#    makeAndSavePlot(axis, histo, model, name="provafit_composite.png")


if __name__ == '__main__':

    custom_pdfs = ['RooCBExGaussShape',
                   'RooDoubleCBFast', 'RooCMSShape', 'my_double_CB']

    import_pdf_library(custom_pdfs[2])

    # Import of the 3D histograms
    f = ROOT.TFile("root_files/tnp_iso_data.root")
    histo_pass = f.pass_mu_RunGtoH
    histo_fail = f.fail_mu_RunGtoH
    
    # Histogram in the first bin (eta, pt)
    dh, x, n_events = import_histo(histo_pass, [1], [1])


    #  -----------------------------------------------------------------------
    # | ~~~~~~~~~~ PDF definition ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ |
    #  -----------------------------------------------------------------------

    y = ROOT.RooRealVar("y", "y", -5, 5)    

    # General parameters - used everywhere
    mean = ROOT.RooRealVar("b", "b", 0, -5, 5)
    sigma = ROOT.RooRealVar("c", "c", 1, 0.1, 5)
    gauss = ROOT.RooGaussian("gauss", "gauss", y, mean, sigma)

    # Crystal Ball PDF
    a1 = ROOT.RooRealVar("d", "d", 1.5, 0, 3)
    n1 = ROOT.RooRealVar("e", "e", 1, -5, 5)
    a2 = ROOT.RooRealVar("f", "f", 1.5, 0, 3)
    n2 = ROOT.RooRealVar("g", "g", 1, 5, 5)
    cb_shape = ROOT.RooCBShape("cb", "cb", x, mean, sigma, a1, n1)

    # Exponential PDF
    tau = ROOT.RooRealVar("tau", "tau", -0.5, -2, 2)
    expo = ROOT.RooExponential("expo", "expo", y, tau)

    # Voigtian PDF
    gamma = ROOT.RooRealVar("gamma", "gamma", 1, 0.5, 3)
    voigt = ROOT.RooVoigtian("voigt", "voigt", x, mean, gamma, sigma)

    # Polynomial PDF (order 3)
    a2 = ROOT.RooRealVar("a2", "a2", -10, 10)
    a3 = ROOT.RooRealVar("a3", "a3", -10, 10)
    pol = ROOT.RooPolynomial("pol", "pol", x, [a2, a3])

    # CMS shape
    alpha = ROOT.RooRealVar("alpha", "alpha", 0.5, -10, 10)
    beta = ROOT.RooRealVar("beta", "beta", 1, 0, 2)
    gamma = ROOT.RooRealVar("gamma", "gamma", 6, 0, 10)
    peak = ROOT.RooRealVar("peak", "peak", 2, -10, 10)
    cmsshape = ROOT.RooCMSShape("func", "func", y, alpha, beta, gamma, peak)

    #  -----------------------------------------------------------------------
    # | ~~~~~~~~~~ Fit and Plot section ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ |
    #  -----------------------------------------------------------------------

    c = ROOT.TCanvas()
    c.cd()

    Nsig = ROOT.RooRealVar("nsig", "#signal events", 0, 10000)
    Nbkg = ROOT.RooRealVar("nbkg", "#background events", 0, 10000)
    
    
    sum_func = ROOT.RooAddPdf("sum", "sum", [gauss, expo], [Nsig, Nbkg])

    # res = sum_func.fitTo(dh, Save=True)

    data = sum_func.generate({y}, 10000)
   
    frame = y.frame("Gauss+expo")
    # sum_func.plotOn(frame)
    # data.plotOn(frame)
    # model = ROOT.RooAddPdf(sum_func)
    r = sum_func.fitTo(data, Save=True, Extended=True, Verbose=False)
    # sum_func.plotOn(frame)
    data.plotOn(frame)
    sum_func.plotOn(frame, VisualizeError=(r, 2))
    sum_func.plotOn(frame)
    # model.plotOn(frame, Components="expo", LineStyle=':')
    sum_func.paramOn(frame)
    
    frame.Draw()
    #alpha.Print()
    #beta.Print()
    #gamma.Print()
    #peak.Print()
    c.SaveAs("gauss-expo_plot.png")

    # makeAndSavePlot(x, dh, sum_func, name="provafit_voigtian.png")
