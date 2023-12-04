import ROOT

from utilities.base_library import bin_dictionary

from utilities.dataset_utils import import_pdf_library

import_pdf_library("RooCMSShape")

bin_dict = bin_dictionary("pt_tracking", "eta")

filen = ROOT.TFile("/scratchnvme/rforti/tnp_efficiencies_results/tracking/numEst_fail/ws_tracking_indep_numEstF.root")
ws = filen.Get("w")

statuses = []

for bin_key in bin_dict.keys():
    
    axis = ws.var(f"x_fail_{bin_key}")
    
    alpha = ROOT.RooRealVar("a", "a", 60, 40, 130)
    beta = ROOT.RooRealVar("b", "b", 5, 0.1, 40)
    gamma = ROOT.RooRealVar("g", "g", 0.1, 0.0, 1.0)
    peak = ROOT.RooRealVar("p", "p", 91)

    cmsshape = ROOT.RooCMSShape("cmsshape", "cmsshape", axis, alpha, beta, gamma, peak)

    mu_alpha = ROOT.RooRealVar("mu_alpha", "mu_alpha", 90)
    sigma_alpha = ROOT.RooRealVar("sigma_alpha", "sigma_alpha", 50)
    constr_alpha = ROOT.RooGaussian("constr_alpha", "constr_alpha", alpha, mu_alpha, sigma_alpha)

    mu_beta = ROOT.RooRealVar("mu_beta", "mu_beta", 5)
    sigma_beta = ROOT.RooRealVar("sigma_beta", "sigma_beta", 2)
    constr_beta = ROOT.RooGaussian("constr_beta", "constr_beta", beta, mu_beta, sigma_beta)

    mu_gamma = ROOT.RooRealVar("mu_gamma", "mu_gamma", 0.1)
    sigma_gamma = ROOT.RooRealVar("sigma_gamma", "sigma_gamma", 0.05)
    constr_gamma = ROOT.RooGaussian("constr_gamma", "constr_gamma", gamma, mu_gamma, sigma_gamma)

    constr = ROOT.RooArgSet(constr_alpha, constr_beta, constr_gamma)
    
    data = ws.data(f"Minv_bkg_fail_{bin_key}_total")
    
    res = cmsshape.fitTo(data, Minimizer="Minuit2", Save=True, Strategy=2, Verbose=0, PrintLevel=-1, 
                         SumW2Error=False, ExternalConstraints=constr)
    

    c = ROOT.TCanvas()
    c.cd()

    frame = axis.frame()
    data.plotOn(frame)
    cmsshape.plotOn(frame)
    cmsshape.paramOn(frame)
    frame.Draw()
    c.SaveAs(f"/scratchnvme/rforti/tnp_efficiencies_results/tracking/fit_cmsshape/cmsshape_fail_{bin_key}.pdf")
    
    statuses.append(res.status())


    
    
print(statuses)
    

