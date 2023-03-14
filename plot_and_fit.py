"""
"""

import ROOT
from utilities import import_pdf_library, import_Steve_histos, makeAndSavePlot


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

    type_eff = ("sa", "global", "ID", "iso", "trigger", "veto")
    t = type_eff[3]

    h_data, h_mc, x = import_Steve_histos(t, [1], [1])

    #  -----------------------------------------------------------------------
    # | ~~~~~~~~~~ PDF definition ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ |
    #  -----------------------------------------------------------------------

    y = ROOT.RooRealVar("y", "y", -5, 5)
    z = ROOT.RooRealVar("z", "z", 50, 70)

    # General parameters - used everywhere
    mean = ROOT.RooRealVar("mean", "mean", 91, 80, 101)
    sigma = ROOT.RooRealVar("sigma", "sigma", 2, 0.2, 5)
    gauss = ROOT.RooGaussian("gauss", "gauss", x, mean, sigma)

    # Crystal Ball PDF
    a1 = ROOT.RooRealVar("d", "d", 1.5, 0, 3)
    n1 = ROOT.RooRealVar("e", "e", 1, -5, 5)
    a2 = ROOT.RooRealVar("f", "f", 1.5, 0, 3)
    n2 = ROOT.RooRealVar("g", "g", 1, 5, 5)
    cb_shape = ROOT.RooCBShape("cb", "cb", x, mean, sigma, a1, n1)

    # Breit-Wigner PDF
    gamma = ROOT.RooRealVar("gamma", "gamma", 2.5, 0.5, 8)
    breitwigner = ROOT.RooBreitWigner(
            "breitwigner", "breitwigner", x, mean, gamma)

    # Exponential PDF
    tau = ROOT.RooRealVar("tau", "tau", -0.5, -2, 2)
    expo = ROOT.RooExponential("expo", "expo", y, tau)

    # Voigtian PDF
    gamma = ROOT.RooRealVar("gamma", "gamma", 1, 0.5, 3)
    voigt = ROOT.RooVoigtian("voigt", "voigt", x, mean, gamma, sigma)

    # Polynomial PDF (order 3)
    a1 = ROOT.RooRealVar("a1", "a1", -1e-2, 1e-2)
    a2 = ROOT.RooRealVar("a2", "a2", -2, 2)
    pol = ROOT.RooPolynomial("pol", "pol", x, [a1])

    # CMS shape
    alpha = ROOT.RooRealVar("alpha", "alpha", 0.5, -10, 10)
    beta = ROOT.RooRealVar("beta", "beta", 1, 0, 2)
    gamma = ROOT.RooRealVar("gamma", "gamma", 6, 0, 10)
    peak = ROOT.RooRealVar("peak", "peak", 2, -10, 10)
    cmsshape = ROOT.RooCMSShape("func", "func", y, alpha, beta, gamma, peak)

    #  -----------------------------------------------------------------------
    # | ~~~~~~~~~~ Fit and Plot section ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ |
    #  -----------------------------------------------------------------------

    Nsig = ROOT.RooRealVar("nsig", "#signal events", 0, h_data[0].GetEntries())
    Nbkg = ROOT.RooRealVar("nbkg", "#background events",
                           0, h_data[0].GetEntries())

    sum_func = ROOT.RooAddPdf("sum", "sum", [gauss, pol], [Nsig, Nbkg])

    # res = sum_func.fitTo(dh, Save=True)

    # data = sum_func.generate({y}, 10000)
    model = ROOT.RooAddPdf(sum_func)
    model.fitTo(h_data[0], Save=True, Extended=True,
                Verbose=False, Hesse=False)

    makeAndSavePlot(x, h_data[0], model, name="provafit_voigtian.png")
    # frame = x.frame("Expo_bkg")
    # sum_func.plotOn(frame)
    # data.plotOn(frame)

    # pol.setStringAttribute("fitrange", "")
    '''
    h_data[0].plotOn(frame)
    model.plotOn(frame)
    model.plotOn(frame, Components="pol", LineStyle=':')
    '''

    # sum_func.plotOn(frame)
    # model.plotOn(frame, Components="expo", LineStyle=':')
    # sum_func.paramOn(frame)

    # frame.Draw()
    #alpha.Print()
    #beta.Print()
    #gamma.Print()
    #peak.Print()
    c.SaveAs("prova.png")

    #
