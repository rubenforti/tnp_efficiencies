
import ROOT
import time
import unittest


NDATA = 100000
mu_gen = 91.0
sigma_gen = 4.0


def generate_data(seed=0):
    """
    """
    ROOT.RooRandom.randomGenerator().SetSeed(seed)

    axis = ROOT.RooRealVar("x", "x", 50, 130)
    axis.setBins(100)
    axis.setRange(50, 130)

    mu = ROOT.RooRealVar("mu", "mu", mu_gen, 85.0, 97.0)  # 91
    sigma = ROOT.RooRealVar("sigma", "sigma", sigma_gen, 0.5, 10.0)  # 4

    func = ROOT.RooGaussian("gauss_bkg", "Gaussian background",
                                  axis, mu, sigma)

    data = func.generateBinned(ROOT.RooArgSet(axis), NDATA)
    data.SetName("Data_gauss")

    return data


class TestSimpleFit(unittest.TestCase):

    def test_single_gausfit(self):

        data = generate_data(seed=42)

        axis = ROOT.RooRealVar("x", "x", 50, 130)
        axis.setBins(100)
        axis.setRange(50, 130)

        mu = ROOT.RooRealVar("mu", "mu", 100.0, 50.0, 130.0)  # 91
        sigma = ROOT.RooRealVar("sigma", "sigma", 7.5, 0.5, 10.0)  # 4

        gaus_fitfunc = ROOT.RooGaussian("gaus_fitfunc", "Gaussian background for fit",
                                         axis, mu, sigma)
        
        expected_num = ROOT.RooRealVar("nexp", "nexp", NDATA, 0.8*NDATA, 1.2*NDATA)
        
        model = ROOT.RooExtendPdf("Extended", "extend", gaus_fitfunc, expected_num)

        res = model.fitTo(data, 
                          ROOT.RooFit.Extended(True), 
                          ROOT.RooFit.Minimizer("Minuit2"), 
                          ROOT.RooFit.Save(True), 
                          ROOT.RooFit.PrintLevel(-1))

        self.assertEqual(res.status(), 0)
        self.assertEqual(res.covQual(), 3)
        self.assertLess(res.edm(), 1e-4)

        pars = res.floatParsFinal()
        self.assertAlmostEqual(pars.find("mu").getVal(), mu_gen, delta=pars.find("mu").getError())
        self.assertAlmostEqual(pars.find("sigma").getVal(), sigma_gen, delta=pars.find("sigma").getError())
        self.assertAlmostEqual(pars.find("nexp").getVal(), NDATA, delta=pars.find("nexp").getError())
    

    def test_doublefit_same_data(self):

        data = generate_data(seed=42)

        expected_num = ROOT.RooRealVar("nexp", "nexp", NDATA, 0.8*NDATA, 1.2*NDATA)

        axis1 = ROOT.RooRealVar("x2", "x2", 50, 130)
        axis1.setBins(100)
        axis1.setRange(50, 130)

        mu1 = ROOT.RooRealVar("mu1", "mu", 100.0, 50.0, 130.0)  # 91
        sigma1 = ROOT.RooRealVar("sigma1", "sigma", 7.5, 0.5, 10.0)  # 4

        gaus_fitfunc1 = ROOT.RooGaussian("gaus_fitfunc1", "Gaussian background for fit",
                                         axis1, mu1, sigma1)
        
        model1 = ROOT.RooExtendPdf("Extended1", "extend", gaus_fitfunc1, expected_num)

        res1 = model1.fitTo(data, 
                            ROOT.RooFit.Extended(True), 
                            ROOT.RooFit.Minimizer("Minuit2"), 
                            ROOT.RooFit.Save(True), 
                            ROOT.RooFit.PrintLevel(-1))
    
        axis2 = ROOT.RooRealVar("x2", "x2", 50, 130)
        axis2.setBins(100)
        axis2.setRange(50, 130)

        mu2 = ROOT.RooRealVar("mu2", "mu2", 100.0, 50.0, 130.0)  # 91
        sigma2 = ROOT.RooRealVar("sigma2", "sigma2", 7.5, 0.5, 10.0)  # 4

        gaus_fitfunc2 = ROOT.RooGaussian("gaus_fitfunc2", "Gaussian background for fit",
                                         axis2, mu2, sigma2)
        
        model2 = ROOT.RooExtendPdf("Extended2", "extend", gaus_fitfunc2, expected_num)

        res2 = model2.fitTo(data, 
                            ROOT.RooFit.Extended(True), 
                            ROOT.RooFit.Minimizer("Minuit2"), 
                            ROOT.RooFit.Save(True), 
                            ROOT.RooFit.PrintLevel(-1))

        pars1 = res1.floatParsFinal()
        pars2 = res2.floatParsFinal()

        self.assertAlmostEqual(pars1.find("mu1").getVal(), pars2.find("mu2").getVal())
        self.assertAlmostEqual(pars1.find("sigma1").getVal(), pars2.find("sigma2").getVal())
        self.assertAlmostEqual(pars1.find("nexp").getVal()/NDATA, pars2.find("nexp").getVal()/NDATA)


if __name__ == '__main__':
    unittest.main()

