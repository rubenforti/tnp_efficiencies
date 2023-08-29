
import ROOT
from utilities.fit_utils import llr_eval, pearson_chi2_eval
import matplotlib.pyplot as plt
import numpy as np

axis = ROOT.RooRealVar("x", "x", 50, 130)

def compare_nll(NBINS):

    EVTS = int(NBINS*100)
    nevents = ROOT.RooRealVar("numev", "numev", EVTS, 0, 5*EVTS)

    axis.setBins(NBINS)
    binning = axis.getBinning()

    
    mu = ROOT.RooRealVar("mu", "mu", 82, 80, 100)
    sigma = ROOT.RooRealVar("sigma", "sigma", 2.5, 0.5, 5)
    gaus = ROOT.RooGaussian("gaus", "gaus", axis, mu, sigma)
    '''
    tau = ROOT.RooRealVar("tau", "tau", -1, -2, 0)  
    expo = ROOT.RooExponential("expo", "expo", axis, tau)
    '''

    fit_func = ROOT.RooExtendPdf("fit_func", "fit_func", gaus, nevents)
    data = fit_func.generateBinned(ROOT.RooArgSet(axis), EVTS)

    '''
    root_nll = fit_func.createNLL(data, ROOT.RooFit.Extended(1))
    nll_val = root_nll.getVal()
    minimizer = ROOT.RooMinimizer(root_nll)
    minimizer_func = ROOT.RooMinimizerFcn(root_nll, minimizer)
    '''

    res = fit_func.fitTo(data, 
                    ROOT.RooFit.Save(1), 
                    ROOT.RooFit.Minimizer("Minuit2", "Migrad"),
                    ROOT.RooFit.PrintLevel(-1))
                    #ROOT.RooFit.Extended(1))

    '''
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
    '''

    sum_ll = 0
    max_ll = 0
    for i in range(NBINS):

        axis.setVal(binning.binCenter(i))
        weight = data.weight(i)

        pdf_val = fit_func.getVal(ROOT.RooArgSet(axis))
        mu = EVTS*(axis.getMax()-axis.getMin())*pdf_val/NBINS
        mu = round(mu, 5)

        new_sumll = weight*ROOT.TMath.Log(mu) - mu if mu>0 else 0.0
        sum_ll += 2*new_sumll

        if weight > 0:
            max_ll += 2*weight*ROOT.TMath.Log(weight) - 2*weight

    sum_ll += 2*EVTS*ROOT.TMath.Log(int(nevents.getVal())) - 2*nevents.getVal()
    max_ll += 2*EVTS*ROOT.TMath.Log(EVTS) - 2*EVTS


    '''
    roostat_var = ROOT.RooStats.MinNLLTestStat(fit_func)
    pars = ROOT.RooArgSet()
    for par in res.floatParsFinal():
        #par.setConstant(True)
        pars.addClone(par)
    chi2_roostat = -2*roostat_var.Evaluate(data, pars)

    '''    

    chi2_roostat = -2*res.minNll()

    print(max_ll, sum_ll)
    return sum_ll-chi2_roostat


if __name__=="__main__":

    bins = [20 + 5*i for i in range(10)]
    nlls = [compare_nll(b) for b in bins]

    bins = np.array(bins)
    nlls = np.array(nlls)

    ang_coeff = (np.max(nlls)-np.min(nlls))/(np.max(bins)-np.min(bins))
    offset = np.min(nlls) - ang_coeff*np.min(bins)
    print(ang_coeff, offset)


    exp_nlls = [ang_coeff*b for b in bins]

    plt.figure(1)
    plt.plot(bins, nlls, marker="o", linestyle="")
    plt.plot(bins, exp_nlls, marker="", linestyle="--")
    plt.grid()
    plt.show()
    print(nlls)
    
    '''
    plt.figure(2)
    plt.hist((nlls-exp_nlls)/nlls, bins=10)
    plt.show()
    '''