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
    # | ~~~~~~~~~~ Fit and Plot section ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ |
    #  -----------------------------------------------------------------------

    Nsig = ROOT.RooRealVar("nsig", "#signal events", 0, 10000)
    Nbkg = ROOT.RooRealVar("nbkg", "#background events", 0, 10000)

    sum_func = ROOT.RooAddPdf("sum", "sum", [gauss, expo], [Nsig, Nbkg])

    y.setBins(10000, name='cache')
    conv_func = ROOT.RooFFTConvPdf("conv", "conv", y, breitwigner, gauss, 2)

    name = conv_func.Class_Name()
    print(name)
    # res = sum_func.fitTo(dh, Save=True)

    data = conv_func.generate({y}, 10000)
    print(data.Class_Name())
    frame = y.frame("Gauss+expo")

    # sum_func.plotOn(frame)
    # data.plotOn(frame)

    # model = ROOT.RooAbsPdf(conv_func)

    r = conv_func.fitTo(data, Save=True, Verbose=False)
    # sum_func.plotOn(frame)
    data.plotOn(frame)
    conv_func.plotOn(frame, VisualizeError=(r, 2))
    # model.plotOn(frame)
    # model.plotOn(frame, Components="expo", LineStyle=':')
    conv_func.paramOn(frame)

    frame.Draw()
    #alpha.Print()
    #beta.Print()
    #gamma.Print()
    #peak.Print()
    c.SaveAs("fit/convolution_pdf_plot.png")

    # makeAndSavePlot(x, dh, sum_func, name="provafit_voigtian.png")
