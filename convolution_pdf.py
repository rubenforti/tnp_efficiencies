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
    
    return model, res





if __name__ == '__main__':

    custom_pdfs = ['RooCBExGaussShape',
                   'RooDoubleCBFast', 'RooCMSShape', 'my_double_CB']
    # import_pdf_library(custom_pdfs[2])

    type_eff = ("sa", "global", "ID", "iso", "trigger", "veto")
    t = type_eff[3]
    
    idx_cond = 1  # Zero for pass, one for fail (HAS TO BE FIXED, IS AWFUL)
    
    NBINS_MASS = 80

    h_data, h_mc, x = import_Steve_histos(t, [1], [1])
    print(type(h_data[idx_cond]), type(h_mc[idx_cond]))
    print(f"Num RooDataHist entries = {h_data[idx_cond].numEntries()}")

    pdf_mc = ROOT.RooHistPdf("pdf_mc", "pdf_mc", x, h_mc[idx_cond])
    
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

    model, res = fit_convolution(x, h_data[idx_cond], pdf_mc, smearing)
    
    npars = NBINS_MASS - res.floatParsFinal().getSize()
    chi2_sqrtvar = (2*npars)**(1/2)
    print(f"Expected chi2 pars: mu={npars}, sqrt(var)={chi2_sqrtvar}")

    chi2_obj = ROOT.RooChi2Var("chi2", "chi2", model, h_data[idx_cond], Verbose=True)
    print(f"Measured chi2 = {chi2_obj.getVal()}")

    print(f"Distance in sigma = {(chi2_obj.getVal()-npars)/chi2_sqrtvar}")
    
    id_flag = "pass" if idx_cond==1 else "fail"

    makeAndSavePlot(x, h_data[idx_cond], model, name=f"figs/iso/fit_{id_flag}_smearing_nobkg.pdf", pull=False)

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
