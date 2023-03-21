"""
"""

import sys
import ROOT


def pearson_chi2_eval(histo, pdf, nbins, res):
    """
    """
    npars = nbins - res.floatParsFinal().getSize()
    print(npars)
    chi2_sqrtvar = (2*npars)**(1/2.)
    print(f"Expected chi2 pars: mu={npars}, sqrt(var)={chi2_sqrtvar}")

    chi2_obj = ROOT.RooChi2Var("chi2", "chi2", pdf, histo, Verbose=False)
    print(f"Measured chi2 = {chi2_obj.getVal()}")
    print(f"Distance in sigma = {(chi2_obj.getVal()-npars)/chi2_sqrtvar}")


def llr_test_bkg(histo, pdf):
    """
    """

    pars_set = pdf.getParameters(histo)

    profiled_var = ROOT.RooRealVar()

    for idx in range(pars_set.getSize()):
        if pars_set[idx].GetName() == 'nbkg':
            profiled_var = pars_set[idx]

    null_profiled_var = ROOT.RooRealVar("nbkg", "nbkg", 0, 0, 0)

    llr_obj = ROOT.RooStats.ProfileLikelihoodCalculator(
            histo, pdf, ROOT.RooArgSet(profiled_var), 0.05, ROOT.RooArgSet(null_profiled_var))

    test_res = llr_obj.GetHypoTest()

    pval = test_res.NullPValue()

    print(f"p-value for null hypo is: {pval}")

    null = True if pval > 0.05 else False

    return null
