

import ROOT
from utilities import import_pdf_library


def enableBinIntegrator(func, num_bins):
    """
    Force numeric integration and do this numeric integration with the
    RooBinIntegrator, which sums the function values at the bin centers.
    """
    custom_config = ROOT.RooNumIntConfig(func.getIntegratorConfig())
    custom_config.method1D().setLabel("RooBinIntegrator")
    custom_config.getConfigSection(
        "RooBinIntegrator").setRealValue("numBins", num_bins)
    func.setIntegratorConfig(custom_config)
    func.forceNumInt(True)


import_pdf_library('RooCMSShape')

axis = ROOT.RooRealVar("x", "x", 50, 130)
axis.setBins(1000)

alpha = ROOT.RooRealVar("alpha", "alpha", 102, 40.0, 130.0)
beta = ROOT.RooRealVar("beta", "beta", 5.0, 0.1, 40.0)
gamma = ROOT.RooRealVar("gamma", "gamma", 0.67, 0, 1)
peak = ROOT.RooRealVar("peak", "peak", 90.0)  # 88.0, 92.0)

#background = ROOT.RooCMSShape("cmsshape_bkg", "CMSShape background",
#                              axis, alpha, beta, gamma, peak)
#enableBinIntegrator(background, 500)

axis.setRange("fitRange", 75, 105)


mu = ROOT.RooRealVar("mu", "mu", 91, 85, 97)
sigma = ROOT.RooRealVar("sigma", "sigma", 4, 0.5, 10)
background = ROOT.RooGaussian("gaus", "gaus", axis, mu, sigma)


expected_num = ROOT.RooRealVar("nexp", "nexp", 10000, 0, 20000)
extend = ROOT.RooExtendPdf(
    "Extended", "extend", background, expected_num, rangeName="fitRange")

data = background.generate({axis}, 10000)


model = ROOT.RooExtendPdf(extend)
# enableBinIntegrator(model, 1000)

res = model.fitTo(data,
                  # NumCPU=(10, 1),
                  Extended=True,
                  Range="fitRange",
                  #Constrain=ROOT.RooArgSet("nexp"),
                  Save=True,
                  IntegrateBins=1e-6,
                  #Minimizer=("Minuit2", "migrad"),
                  PrintLevel=1)

res.Print()
res.correlationMatrix().Print()

print(res.status())
print(res.covQual())
# covm.Print()
