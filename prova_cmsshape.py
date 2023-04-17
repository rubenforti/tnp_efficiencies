

import ROOT
from utilities import import_pdf_library

import_pdf_library('RooCMSShape')

axis = ROOT.RooRealVar("x", "x", 50, 130)

alpha = ROOT.RooRealVar("alpha", "alpha", 60.0, 40.0, 130.0)
beta = ROOT.RooRealVar("beta", "beta", 5.0, 0.1, 40.0)
gamma = ROOT.RooRealVar("gamma", "gamma", 0.1, 0, 1)
peak = ROOT.RooRealVar("peak", "peak", 90.0)  # 88.0, 92.0)

background = ROOT.RooCMSShape("cmsshape_bkg", "CMSShape background",
                              axis, alpha, beta, gamma, peak)

expected_num = ROOT.RooRealVar("nexp", "nexp", 10000, 0, 20000)
extend = ROOT.RooExtendedTerm("extend", "extend", expected_num)

data = background.generate({axis}, 10000)

model = ROOT.RooProduct("model", "model", [extend, background])

res = model.fitTo(data, Extended=True, Save=True)

pars = res.floatParsFinal()
pars.Print("v")
