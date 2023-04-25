"""
"""

import sys
import os
import pickle
import ROOT


def import_pdf_library(*functions):
    """
    """
    current_path = os.path.dirname(__file__)

    for function in functions:
        header_incl = ' #include "libCpp/'+function+'.h"'
        sourcefile = os.path.join(current_path, 'libCpp', f'{function}.cc')

        print(sourcefile)
        ctrl_head = ROOT.gInterpreter.Declare(header_incl)
        ctrl_source = ROOT.gSystem.CompileMacro(sourcefile, opt="ks")

        if ctrl_head is not True:
            print("ERROR in header loading")
            sys.exit()
        if ctrl_source != 1:
            print("ERROR in sourcefile compiling and loading")
            sys.exit()


def fit_quality(res):
    """
    """
    check_migrad = (res.status() == 0)
    check_covm = (res.covQual() == 3)
    check_edm = (res.edm() < 1e-4)

    return bool(check_migrad*check_covm*check_edm)


def eval_efficiency(npass, nfail, sigma_npass, sigma_nfail):
    """
    """
    eff = npass/(npass+nfail)
    var1 = ((1-npass)**2)*(sigma_npass**2)
    var2 = (npass**2)*(sigma_nfail**2)
    sigma_eff = ROOT.TMath.Sqrt(var1+var2)/((npass+nfail)**2)

    return eff, sigma_eff


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


def llr_test_bkg(histo, pdf, alpha=0.05):
    """
    """

    pars_set = pdf.getParameters(histo)

    profiled_var = ROOT.RooRealVar()

    for idx in range(pars_set.getSize()):
        if pars_set[idx].GetName() == 'nbkg':
            profiled_var = pars_set[idx]

    null_profiled_var = ROOT.RooRealVar("nbkg", "nbkg", 0, 0, 0)

    llr_obj = ROOT.RooStats.ProfileLikelihoodCalculator(
        histo, pdf, ROOT.RooArgSet(profiled_var), alpha,
        ROOT.RooArgSet(null_profiled_var))

    test_res = llr_obj.GetHypoTest()

    pval = test_res.NullPValue()

    print(f"p-value for null hypo is: {pval}")

    null = True if pval > 0.05 else False

    return null


if __name__ == '__main__':
    print("CIAO")
