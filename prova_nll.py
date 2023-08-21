
import ROOT
from utilities.fit_utils import llr_eval, pearson_chi2_eval
import matplotlib.pyplot as plt
import numpy as np

axis = ROOT.RooRealVar("x", "x", 60, 120)

def compare_nll(NBINS):

    axis.setBins(NBINS)

    mu = ROOT.RooRealVar("mu", "mu", 91, 85, 95)
    sigma = ROOT.RooRealVar("sigma", "sigma", 2.5, 0.5, 5)

    gauss = ROOT.RooGaussian("gauss", "gauss", axis, mu, sigma)

    EVTS = int(NBINS*100)

    nevents = ROOT.RooRealVar("numev", "numev", EVTS, 0, 5*EVTS)

    gaus_extended = ROOT.RooExtendPdf("gaus_extended", "gaus_extended", gauss, nevents)

    data = gaus_extended.generateBinned(ROOT.RooArgSet(axis), EVTS)
    data_unbinned = gaus_extended.generate(ROOT.RooArgSet(axis), EVTS)

    binning = axis.getBinning()
    bin_centers = [binning.binCenter(i) for i in range(NBINS)]
    print(bin_centers)

    weights = [data.weight(i) for i in range(NBINS)]
    print(weights)

    root_nll = gaus_extended.createNLL(data, ROOT.RooFit.Extended(1))
    nll_val = root_nll.getVal()
    minimizer = ROOT.RooMinimizer(root_nll)
    minimizer_func = ROOT.RooMinimizerFcn(root_nll, minimizer)


    res = gaus_extended.fitTo(data, 
                    ROOT.RooFit.Save(1), 
                    ROOT.RooFit.Minimizer("Minuit2", "Migrad"),
                    ROOT.RooFit.PrintLevel(3))
                    #ROOT.RooFit.Extended(1))



    def getGaus(x):

        expo = ROOT.TMath.Exp(-0.5*(x-mu.getVal())*(x-mu.getVal())/(sigma.getVal()*sigma.getVal()))
        denom = 1./ROOT.TMath.Sqrt(2*ROOT.TMath.Pi()*sigma.getVal()*sigma.getVal())

        return expo*denom


    def getLogLikelihood(weight, x):
        mu = EVTS*(axis.getMax()-axis.getMin())*getGaus(x)/NBINS
        if mu > 0:
            return weight*ROOT.TMath.Log(mu) - mu
        else:
            return 0.0

    sum_ll = 0
    max_ll = 0
    for i in range(NBINS):
        sum_ll += 2*getLogLikelihood(weights[i], bin_centers[i])
        if weights[i] > 0:
            max_ll += 2*weights[i]*ROOT.TMath.Log(weights[i]) - 2*weights[i]


    sum_ll += 2*EVTS*ROOT.TMath.Log(int(nevents.getVal())) - 2*nevents.getVal()
    max_ll += 2*EVTS*ROOT.TMath.Log(EVTS) - 2*EVTS

    # chi2val, ndof = pearson_chi2_eval(data, gaus_extended, NBINS, res)



    hist_pdf = ROOT.RooHistPdf("roohistpdf", "roohistpdf", ROOT.RooArgSet(axis), data, 3)

    roostat_var = ROOT.RooStats.MinNLLTestStat(gaus_extended)


    pars = ROOT.RooArgSet()

    for par in res.floatParsFinal():
        #par.setConstant(True)
        pars.addClone(par)

    chi2_roostat = -2*roostat_var.Evaluate(data, pars)

    offset = minimizer_func.getOffset()
    print(offset)
    
    '''
    print(max_ll, sum_ll)
    print((max_ll-sum_ll)/(2*(NBINS-3)))
    print(chi2_roostat)
    '''
    return (sum_ll-chi2_roostat)


if __name__=="__main__":

    bins = [30 + 5*i for i in range(1)]
    nlls = [compare_nll(b) for b in bins]

    bins = np.array(bins)
    nlls = np.array(nlls)

    ang_coeff = (np.max(nlls)-np.min(nlls))/(np.max(bins)-np.min(bins))
    # print(ang_coeff)


    dev = []
    for i in range(len(nlls)):
        dev.append(nlls[i] - ang_coeff*bins[i])

    exp_nlls = [ang_coeff*b for b in bins]

    plt.figure(1)
    plt.plot(bins, nlls, marker="o", linestyle="")
    plt.plot(bins, exp_nlls, marker="", linestyle="--")
    plt.grid()
    plt.show()
    print(nlls)
    
    plt.figure(2)
    plt.hist(dev, bins=10)
    plt.show()