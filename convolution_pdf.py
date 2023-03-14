"""
"""

import ROOT
from utilities import import_pdf_library, import_Steve_histos, makeAndSavePlot


def fit_convolution(axis, histo, template_pdf, smearing, nbins=1000, buffer_frac=0.1, int_order=3):
    """
    """
    axis.setBins(1000, "cache")
    conv_func = ROOT.RooFFTConvPdf("conv", "conv", axis, template_pdf, smearing, int_order)
    conv_func.setBufferFraction(buffer_frac)
    name = conv_func.Class_Name()
    print(name)
    model = ROOT.RooFFTConvPdf(conv_func)
    
    res = model.fitTo(histo, Save=True, Verbose=False)
    
    return model





if __name__ == '__main__':

    custom_pdfs = ['RooCBExGaussShape',
                   'RooDoubleCBFast', 'RooCMSShape', 'my_double_CB']
    # import_pdf_library(custom_pdfs[2])

    type_eff = ("sa", "global", "ID", "iso", "trigger", "veto")
    t = type_eff[3]

    h_data, h_mc, x = import_Steve_histos(t, [1], [1])
    print(type(h_data[0]), type(h_mc[0]))
    print(f"Num RooDataHist entries = {h_data[0].numEntries()}")

    pdf_mc = ROOT.RooHistPdf("pdf_mc", "pdf_mc", x, h_mc[0])
    
    mean = ROOT.RooRealVar("mean", "mean", 0, -2, 2)
    sigma = ROOT.RooRealVar("sigma", "sigma", 0.5, 0.001, 2)
    smearing = ROOT.RooGaussian("smearing", "smearing", x, mean, sigma)

        
    x.setBins(1000, "cache")
    '''
    conv_func = ROOT.RooFFTConvPdf("conv", "conv", x, pdf_mc, smearing, 3)
    conv_func.setBufferFraction(0.1)
    name = conv_func.Class_Name()
    print(name)

    # fit_func = ROOT.RooAddPdf("fitfunc", "fitfunc", x

    # res = sum_func.fitTo(dh, Save=True)
    
    #model = pdf_mc
    model = ROOT.RooFFTConvPdf(conv_func) 
    
    # data = pdf_mc.generate({x}, 500000)
    res = model.fitTo(h_data[0], Save=True, Verbose=False)
    '''

    model = fit_convolution(x, h_data[0], pdf_mc, smearing)


    makeAndSavePlot(x, h_data[0], model, name="fit_on_data_smearing_nobkg.png", pull=False)

    '''
    c = ROOT.TCanvas()
    c.cd()
    frame = x.frame()
    # sum_func.plotOn(frame)
    # h_mc[0].plotOn(frame, LineColor="kBlue")
    h_data[0].plotOn(frame, LineColor="kRed")
   # data.plotOn(frame)
    #conv_func.plotOn(frame, VisualizeError=(r, 2))
    model.plotOn(frame)
    # model.plotOn(frame, Components="expo", LineStyle=':')
    
    #alpha.Print()
    #beta.Print()
    #gamma.Print()
    #peak.Print()
    
    frame.Draw()
    c.SaveAs("prova.png")
    '''
