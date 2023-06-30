"""
"""

import ROOT
import unittest


NDATA = 100000
NBINS = 100


def generate_th1_from_roodatahist(sumw2=False):

    axis = ROOT.RooRealVar("x", "x", 50, 130)
    axis.setBins(NBINS)
    axis.setRange(50, 130)

    '''
    mu = ROOT.RooRealVar("mu", "mu", 90, 85.0, 97.0)  # 91
    sigma = ROOT.RooRealVar("sigma", "sigma", 4, 0.5, 10.0)  # 4
    gaus_pdf = ROOT.RooGaussian("gaus", "gaus", axis, mu, sigma)
    '''

    unif_pdf = ROOT.RooUniform("unif", "unif", axis)

    data = unif_pdf.generateBinned(ROOT.RooArgSet(axis), NDATA)

    # histo = ROOT.TH1D("histo", "histo", NBINS, 50, 130)
    

    histo = data.createHistogram("histo", axis)
    

    #histo = ROOT.TH1D("histo", "histo", 100, 0, 1)

    return histo

###############################################################################

def generate_th1_from_randgen(name, weight=1):

    axis = ROOT.RooRealVar("x", "x", 50, 130)
    axis.setBins(NBINS)
    axis.setRange(50, 130)

    '''
    mu = ROOT.RooRealVar("mu", "mu", 90, 85.0, 97.0)  # 91
    sigma = ROOT.RooRealVar("sigma", "sigma", 4, 0.5, 10.0)  # 4
    gaus_pdf = ROOT.RooGaussian("gaus", "gaus", axis, mu, sigma)
    '''
    unif_pdf = ROOT.RooUniform("unif", "unif", axis)

    histo = ROOT.TH1D(name, "histo", NBINS, 50, 130)

    ROOT.RooRandom.randomGenerator().SetSeed(42)

    dataset = unif_pdf.generate(ROOT.RooArgSet(axis), NDATA)

    for i in range(NDATA):
        val_set = dataset.get(i)
        histo.Fill(val_set.find("x").getVal(), weight)

    return histo

###############################################################################

class TestHistoErrors(unittest.TestCase):

    def test_poisson_roohisto(self):
        """
        Check that the errors in a RooDataHist are the sqrt(bincontent)
        """
        histo = generate_th1_from_roodatahist()

        for i in range (1, NBINS+1):
            bin_val = histo.GetBinContent(i)
            bin_error = histo.GetBinError(i)
            self.assertAlmostEqual(bin_error, bin_val**0.5)


    def test_errors_th1(self):
        """
        Checks on how the errors in a TH1 are treated when using some methods like "Add", "Scale",
        "Sumw2" or the histogram is filled with weighted events.
        """

        histo_pois = generate_th1_from_randgen("histo_pois")
        histo_weighted = generate_th1_from_randgen("histo_weighted", weight=2.0)

        histo_sumw2 = histo_pois.Clone("histo_sumw2")
        # histo_sumw2.SetName("histo_sumw2")
        histo_sumw2.Sumw2()

        ## Check that the "Add" method (per se) leaves the bin contents and errors as expected
        histo_pois2 = histo_pois.Clone("histo_pois_2")
        histo_pois2.Add(histo_pois)
        [self.assertAlmostEqual(histo_pois2.GetBinError(i), histo_pois2.GetBinContent(i)**0.5) 
         for i in range(1, NBINS+1)]

        ## Check that with the "Scale" method, the errors returned are also scaled by the factor 
        ## indicated (=2), so that the bin content are NOT the sqrt(bincontent), but are the same of 
        ## the errors present in an histogram filled with weighted events
        histo_sumw2.Scale(2)

        for i in range(1, NBINS+1):
            self.assertAlmostEqual(histo_pois2.GetBinContent(i), histo_sumw2.GetBinContent(i))
            self.assertNotAlmostEqual(histo_sumw2.GetBinError(i), histo_sumw2.GetBinContent(i)**0.5)
            self.assertAlmostEqual(histo_sumw2.GetBinError(i), histo_pois.GetBinError(i)*2)
            self.assertAlmostEqual(histo_sumw2.GetBinContent(i), histo_weighted.GetBinContent(i))
            self.assertAlmostEqual(histo_sumw2.GetBinError(i), histo_weighted.GetBinError(i))

        print(histo_pois.GetBinError(43), histo_pois.GetBinContent(43)**0.5)
        print(histo_pois.GetBinContent(43), histo_sumw2.GetBinContent(43))


###############################################################################

class TestHistoIntegral(unittest.TestCase):

    def test_integral_scaled_th1(self):
        """
        Check that the integral of a TH1 is scaled by the factor indicated in the "Scale" method
        """
        histo_pois = generate_th1_from_randgen("histo_pois")
        histo_pois.Scale(2.0)
        self.assertAlmostEqual(histo_pois.Integral(), 2*NDATA)
        self.assertAlmostEqual((2*histo_pois.Integral())**0.5, 2*(NDATA**0.5))  # dummy, just for clarity


    def test_sumEntries_scaled(self):
        """
        Check that the sumEntries() method of a RooDataHist returns the same value as the integral of the
        corresponding TH1, in the case the TH1 is scaled        """
        h_th1 = generate_th1_from_randgen("h_th1", weight=1.5)
        axis = ROOT.RooRealVar("x", "x", 50, 130)
        h_roohist = ROOT.RooDataHist("h_roohist", "h_roohist", ROOT.RooArgList(axis), h_th1)
        self.assertAlmostEqual(h_roohist.sumEntries(), h_th1.Integral())
    
    
    def test_roohisto_add(self):
        """
        Check that the add() method for a RooDataHist works as expected
        """
        h_th1_1 = generate_th1_from_randgen("h_th1_1", weight=1.5)
        h_th1_2 = generate_th1_from_randgen("h_th1_2", weight=2.7)
        axis = ROOT.RooRealVar("x", "x", 50, 130)
        h_roohist_1 = ROOT.RooDataHist("h_roohist_1", "h_roohist", ROOT.RooArgList(axis), h_th1_1)
        h_roohist_2 = ROOT.RooDataHist("h_roohist_22", "h_roohist", ROOT.RooArgList(axis), h_th1_2)

        h_roohist = ROOT.RooDataHist("h_roohist", "h_roohist", ROOT.RooArgList(axis))
        h_roohist.add(h_roohist_1)
        h_roohist.add(h_roohist_2)

        self.assertAlmostEqual(h_roohist.sumEntries(), h_th1_1.Integral()+h_th1_2.Integral())

        for i in range(0, h_roohist.numEntries()):
            h_roohist.get(i)
            self.assertAlmostEqual(h_roohist.weight(i), h_th1_1.GetBinContent(i+1)+h_th1_2.GetBinContent(i+1))




###############################################################################  
###############################################################################


if __name__ == '__main__':
    unittest.main()