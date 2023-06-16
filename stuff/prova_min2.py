

import ROOT
import os
import sys
from utilities.general_utils import import_pdf_library


def prova_minuit2():
    """
    """

    import_pdf_library('RooCMSShape')

    axis = ROOT.RooRealVar("x", "x", 50, 130)
    axis.setBins(100)
    axis.setRange(50, 130)

    # file_ws = ROOT.TFile("ws_prova.root")
    # ws = file_ws.Get("w")

    # axis = ROOT.RooRealVar(ws["x"])
    # data1 = ROOT.RooDataHist(ws["Data_cmsshape"])
    # data2 = ROOT.RooDataHist(ws["Data_gauss"])

    alpha = ROOT.RooRealVar("alpha", "alpha", 60.0, 40.0, 130.0)
    beta = ROOT.RooRealVar("beta", "beta", 5.0, 0.1, 20)
    gamma = ROOT.RooRealVar("gamma", "gamma", 0.07, 0.0001, 0.2)
    peak = ROOT.RooRealVar("peak", "peak", 91.0)  # 88.0, 92.0)

    mu = ROOT.RooRealVar("mu", "mu", 91.0, 85.0, 97.0)
    sigma = ROOT.RooRealVar("sigma", "sigma", 4.0, 0.5, 10.0)

    background1 = ROOT.RooCMSShape("cmsshape_bkg", "CMSShape background",
                                   axis, alpha, beta, gamma, peak)
    '''
    background2 = ROOT.RooGaussian("gauss_bkg", "Gaussian background",
                                   axis, mu, sigma)
    '''

    NDATA = 100000

    data = background1.generateBinned(ROOT.RooArgSet(axis), NDATA)
    data.SetName("Data_gauss")

    mu2 = ROOT.RooRealVar(
        "mu2", "mu2", 97.0, 50.0, 130.0)  # 91
    sigma2 = ROOT.RooRealVar("sigma2", "sigma2", 6.0, 0.5, 10.0)  # 4

    alpha2 = ROOT.RooRealVar("alpha2", "alpha", 62.0, 40.0, 130.0)
    beta2 = ROOT.RooRealVar("beta2", "beta", 5.4, 0.1, 20)
    gamma2 = ROOT.RooRealVar("gamma2", "gamma", 0.077, 0.0001, 0.2)
    bkg_fit = ROOT.RooCMSShape("cmsshape_bkg", "CMSShape background",
                               axis, alpha2, beta2, gamma2, peak)
    '''
    bkg_fit = ROOT.RooGaussian("bkg_fit", "Gaussian background for fit",
                               axis, mu2, sigma2)
    '''
    expected_num = ROOT.RooRealVar(
        "nexp", "nexp", NDATA, 90000, 110000)

    model = ROOT.RooExtendPdf(
        "Extended", "extend", bkg_fit, expected_num)

    res2 = model.fitTo(
        data, ROOT.RooCmdArg("Minimizer", 0, 0, 0, 0, "Minuit2", "migrad"),
        ROOT.RooCmdArg("Save", 1), 
        ROOT.RooCmdArg("PrintLevel", -1))

    print("***********************************")
    res2.Print()
    res2.correlationMatrix().Print()
    print(res2.status())
    print(res2.covQual())

    '''
    w = ROOT.RooWorkspace("w")
    w.Import(data1)
    w.Import(data2)
    w.writeToFile("ws_prova.root")
    '''

    c = ROOT.TCanvas()
    c.cd()
    frame = axis.frame()
    data.plotOn(frame)
    model.plotOn(frame)
    frame.Draw()
    c.SaveAs("prova.png")
    frame.Delete()
    #print(f"TIME ELAPSED = {t1-t0}")


def prova_roominimizer():

    axis = ROOT.RooRealVar("x", "x", 50, 130)
    axis.setBins(100)
    axis.setRange(50, 130)

    mu = ROOT.RooRealVar("mu", "mu", 91.0, 85.0, 97.0)  # 91
    sigma = ROOT.RooRealVar("sigma", "sigma", 4.0, 0.5, 10.0)  # 4

    background2 = ROOT.RooGaussian("gauss_bkg", "Gaussian background",
                                   axis, mu, sigma)

    NDATA = 100000

    data = background2.generateBinned({axis}, NDATA)
    data.SetName("Data_gauss")

    mu2 = ROOT.RooRealVar("mu2", "mu2", 100.0, 70.0, 110.0)  # 91
    sigma2 = ROOT.RooRealVar("sigma2", "sigma2", 7.0, 0.5, 10.0)  # 4

    bkg_fit = ROOT.RooGaussian("bkg_fit", "Gaussian background for fit",
                               axis, mu2, sigma2)

    expected_num = ROOT.RooRealVar(
        "nexp", "nexp", NDATA, NDATA-10*(NDATA**0.5), NDATA+10*(NDATA**0.5))

    model = ROOT.RooExtendPdf(
        "Extended", "extend", bkg_fit, expected_num)

    nll = ROOT.RooNLLVar("nll", "nll", bkg_fit, data)

    ROOT.Math.MinimizerOptions.SetDefaultMinimizer("Minuit2")

    minimizer = ROOT.RooMinimizer(nll)
    minimizer.setStrategy(2)
    minimizer.setPrintLevel(2)
    minimizer.setMaxFunctionCalls(10000)

    minimizer.migrad()
    minimizer.hesse()

    '''
    minfcn = ROOT.RooMinimizerFcn(nll, minimizer)

    min2 = ROOT.Minuit2.Minuit2Minimizer(0)
    min2.SetFunction(minfcn)
    min2.Minimize()
    min2.PrintResults()

    '''
    res2 = minimizer.lastMinuitFit()
    print("***********************************")
    res2.Print()
    res2.correlationMatrix().Print()
    print(res2.status())
    print(res2.covQual())

    minimizer.cleanup()


if __name__ == '__main__':

    '''
    path = os.path.dirname(__file__)
    ROOT.gSystem.cd(path)
    '''

    ROOT.RooRandom.randomGenerator().SetSeed(42)
    ROOT.gRandom.SetSeed(42)
    prova_minuit2()

    '''
    print(ROOT.RooRandom.randomGenerator().GetSeed())
    print(ROOT.gRandom.GetSeed())

    # prova_roominimizer()

    #sys.exit()
    '''
