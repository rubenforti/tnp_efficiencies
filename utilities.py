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

def pearson_chi2_eval(histo, pdf, nbins, res):
    """
    """
    ndof = nbins - res.floatParsFinal().getSize()
    print(ndof)

    chi2_obj = ROOT.RooChi2Var("chi2", "chi2", pdf, histo)
    chi2_val = chi2_obj.getVal() 

    return chi2_val, ndof


def fit_quality(res):
    """
    """
    check_migrad = (res.status() == 0)
    check_covm = (res.covQual() == 3)
    check_edm = (res.edm() < 1e-4)

    return bool(check_migrad*check_covm*check_edm)

def fit_quality_old(res, chi2stat):
    """
    """
    check_migrad = (res.status() == 0 or res.status() == 1)
    check_covm = (res.covQual() == 3)
    check_chi2 = (chi2stat[0] - chi2stat[1] < 10*((2*chi2stat[1])**0.5))

    return bool(check_migrad*check_covm*check_chi2)


def eval_efficiency(npass, nfail, sigma_npass, sigma_nfail):
    """
    """
    eff = npass/(npass+nfail)
    var1 = ((1-npass)**2)*(sigma_npass**2)
    var2 = (npass**2)*(sigma_nfail**2)
    sigma_eff = ROOT.TMath.Sqrt(var1+var2)/((npass+nfail)**2)

    return eff, sigma_eff


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
