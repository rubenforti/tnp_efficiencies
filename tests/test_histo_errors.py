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

    unif_pdf = ROOT.RooUniform("unif", "unif", axis)

    histo = ROOT.TH1D(name, "histo", NBINS, 50, 130)

    ROOT.RooRandom.randomGenerator().SetSeed(42)

    dataset = unif_pdf.generate(ROOT.RooArgSet(axis), NDATA)

    for i in range(NDATA):
        val_set = dataset.get(i)
        histo.Fill(val_set.find("x").getVal(), weight)

    return histo


def generate_roodatahist():
    
    axis = ROOT.RooRealVar("x", "x", 50, 130)
    axis.setBins(NBINS)
    axis.setRange(50, 130)

    mu = ROOT.RooRealVar("mu", "mu", 90, 85.0, 97.0)  # 91
    sigma = ROOT.RooRealVar("sigma", "sigma", 4, 0.5, 10.0)  # 4
    gaus_pdf = ROOT.RooGaussian("gaus", "gaus", axis, mu, sigma)

    data_gaus = gaus_pdf.generateBinned(ROOT.RooArgSet(axis), NDATA)

    unif_pdf = ROOT.RooUniform("unif", "unif", axis)
    data_unif = unif_pdf.generateBinned(ROOT.RooArgSet(axis), 2*NDATA)

    th1_gaus = data_gaus.createHistogram("th1_gaus", axis)
    th1_gaus.Scale(1.573)

    th1_unif = data_unif.createHistogram("th1_unif", axis)
    th1_unif.Scale(1.234)

    tot_roohisto = ROOT.RooDataHist("tot_roohisto", "tot_roohisto", ROOT.RooArgList(axis))
    tot_roohisto.add(ROOT.RooDataHist("roohisto_gaus", "roohisto_gaus", ROOT.RooArgList(axis), th1_gaus))
    tot_roohisto.add(ROOT.RooDataHist("roohisto_unif", "roohisto_unif", ROOT.RooArgList(axis), th1_unif))


    return tot_roohisto, axis


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

    
    def test_errors_summed_weighted_roohisto(self):
        """
        Check that the bin errors in a Roodatahist produced by subtracting two other
        RooDataHists are the sqrt(sum of the errors squared)
        """

        roohisto_in, axis = generate_roodatahist()

        mu_s = ROOT.RooRealVar("mu_s", "mu_s", 90, 85.0, 97.0)  # 91
        sigma_s = ROOT.RooRealVar("sigma_s", "sigma_s", 10, 0.5, 10.0)  # 4
        gaus_pdf_s = ROOT.RooGaussian("gaus_s", "gaus_s", axis, mu_s, sigma_s)

        subt_histo = gaus_pdf_s.generateBinned(ROOT.RooArgSet(axis), int(NDATA/2.))
        print(type(subt_histo))
        roohisto_out = roohisto_in.Clone("roohisto_out")

        for i in range(0, NBINS):
            subt_histo.get(i)
            wsquared_in = subt_histo.weightSquared(i)
            w_in = subt_histo.weight(i)
            weight, weight_error = subt_histo.weight(i), subt_histo.weightError(ROOT.RooAbsData.SumW2)
            subt_histo.set(i, -1*weight, weight_error)
            # print(w_in, wsquared_in, subt_histo.weight(i), subt_histo.weightSquared(i))


        roohisto_out.add(subt_histo)

        for i in range(0, NBINS):
            roohisto_out.get(i)
            roohisto_in.get(i)
            subt_histo.get(i)

            err_tot = roohisto_out.weightError(ROOT.RooAbsData.SumW2)
            err_exp = (roohisto_in.weightError(ROOT.RooAbsData.SumW2)**2 + subt_histo.weightError(ROOT.RooAbsData.SumW2)**2)**0.5
            self.assertAlmostEqual(err_tot, err_exp)
            self.assertAlmostEqual(roohisto_out.weight(i), roohisto_in.weight(i)+subt_histo.weight(i))
            self.assertAlmostEqual(roohisto_out.weightSquared(i), roohisto_in.weightSquared(i)+subt_histo.weightSquared(i))
            self.assertAlmostEqual(roohisto_out.weightSquared(i), roohisto_out.weightError(ROOT.RooAbsData.SumW2)**2)



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
        corresponding TH1, in the case the TH1 is scaled 
        """
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

        errors_etype0 = 0
        errors_etype1 = 0
        for i in range(0, h_roohist.numEntries()):
            h_roohist.get(i)
            self.assertAlmostEqual(h_roohist.weight(i), h_th1_1.GetBinContent(i+1)+h_th1_2.GetBinContent(i+1))
            errors_etype0 += h_roohist.weightError(0)**2
            errors_etype1 += h_roohist.weightError(1)**2
        
        # errors_etype0 = errors_etype0**0.5
        # errors_etype1 = errors_etype1**0.5

        # print(h_roohist.sumEntries()**0.5, errors_etype0, h_roohist.get_sumw2()**0.5, errors_etype1)
        # print(h_roohist.sumEntries(), h_roohist.get_sumw2())
            




###############################################################################  
###############################################################################


if __name__ == '__main__':
    unittest.main()