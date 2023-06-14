
import ROOT
import time


def prova_min2():
    """
    """

    axis = ROOT.RooRealVar("x", "x", 50, 130)
    axis.setBins(100)
    axis.setRange(50, 130)

    mu = ROOT.RooRealVar("mu", "mu", 91.0, 85.0, 97.0)  # 91
    sigma = ROOT.RooRealVar("sigma", "sigma", 4.0, 0.5, 10.0)  # 4

    background = ROOT.RooGaussian("gauss_bkg", "Gaussian background",
                                  axis, mu, sigma)

    NDATA = 100000

    data = background.generateBinned({axis}, NDATA)
    data.SetName("Data_gauss")

    mu2 = ROOT.RooRealVar("mu2", "mu2", 97.0, 50.0, 130.0)  # 91
    sigma2 = ROOT.RooRealVar("sigma2", "sigma2", 6.0, 0.5, 10.0)  # 4

    bkg_fit = ROOT.RooGaussian("bkg_fit", "Gaussian background for fit",
                               axis, mu2, sigma2)

    expected_num = ROOT.RooRealVar(
        "nexp", "nexp", NDATA, 90000, 110000)

    model = ROOT.RooExtendPdf(
        "Extended", "extend", bkg_fit, expected_num)

    res2 = model.fitTo(data, Extended=True, Minimizer=("Minuit2"), Save=True, PrintLevel=-1)

    print("***********************************")
    res2.Print()
    res2.correlationMatrix().Print()
    print(res2.status())
    print(res2.covQual())

    '''
    c = ROOT.TCanvas()
    c.cd()
    frame = axis.frame(Title="Title")
    data.plotOn(frame)
    model.plotOn(frame)
    frame.Draw()
    c.SaveAs("gaussian.png")
    '''

if __name__ == '__main__':

    ROOT.RooRandom.randomGenerator().SetSeed(42)
    prova_min2()

    print(ROOT.RooRandom.randomGenerator().GetSeed())

