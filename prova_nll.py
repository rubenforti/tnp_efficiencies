
import ROOT
from utilities.fit_utils import llr_eval, pearson_chi2_eval

axis = ROOT.RooRealVar("x", "x", 60, 120)

NBINS = 10
axis.setBins(NBINS)

mu = ROOT.RooRealVar("mu", "mu", 91, 85, 95)
sigma = ROOT.RooRealVar("sigma", "sigma", 2.5, 0.5, 5)

gauss = ROOT.RooGaussian("gauss", "gauss", axis, mu, sigma)

EVTS = 100

nevents = ROOT.RooRealVar("numev", "numev", EVTS, 0, 5*EVTS)

gaus_extended = ROOT.RooExtendPdf("gaus_extended", "gaus_extended", gauss, nevents)

data = gaus_extended.generateBinned(ROOT.RooArgSet(axis), EVTS)

binning = axis.getBinning()
bin_centers = [binning.binCenter(i) for i in range(NBINS)]
print(bin_centers)

weights = [data.weight(i) for i in range(NBINS)]
print(weights)

res = gaus_extended.fitTo(data, 
                  ROOT.RooFit.Save(1), 
                  ROOT.RooFit.Minimizer("Minuit2", "Migrad"))
                  #ROOT.RooFit.Extended(1))



def getGaus(x):

    expo = ROOT.TMath.Exp(-0.5*(x-mu.getVal())*(x-mu.getVal())/(sigma.getVal()*sigma.getVal()))
    denom = 1./ROOT.TMath.Sqrt(2*ROOT.TMath.Pi()*sigma.getVal()*sigma.getVal())

    return expo*denom


def getLogLikelihood(weight, x):
    mu = EVTS*(axis.getMax()-axis.getMin())*getGaus(x)/NBINS
    #mu = EVTS*getGaus(x)/NBINS
    if mu > 0:
        return weight*ROOT.TMath.Log(mu) - mu

sum_ll = 0
max_ll = 0
for i in range(NBINS):
    sum_ll += 2*getLogLikelihood(weights[i], bin_centers[i])
    if weights[i] > 0:
        max_ll += 2*weights[i]*ROOT.TMath.Log(weights[i]) - 2*weights[i] 

#sum_ll = 2*EVTS*ROOT.TMath.Log(nevents.getVal()) - 2*nevents.getVal()
#max_ll += 2*EVTS*ROOT.TMath.Log(EVTS) - 2*EVTS - (2*EVTS*ROOT.TMath.Log(nevents.getVal()) - 2*nevents.getVal())

chi2val, ndof = pearson_chi2_eval(data, gaus_extended, NBINS, res)



print(max_ll, sum_ll)
print(max_ll-sum_ll, chi2val, ndof)
print(-res.minNll())


