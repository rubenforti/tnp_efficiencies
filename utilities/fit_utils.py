"""
"""

import sys
import os
import ROOT
from array import array


def check_existing_fit(type_analysis, ws, bin_key):

    if type_analysis == 'indep':

        isFittedPass = type(ws.obj(f'PDF_pass_{bin_key}')) is ROOT.RooAddPdf
        isFittedFail = type(ws.obj(f'PDF_fail_{bin_key}')) is ROOT.RooAddPdf
        if isFittedPass and isFittedFail:
            print("Not possible to refit an existing PDF! \nReturning the results obtained previously")
            res_pass = ROOT.RooFitResult(ws.obj(f'results_pass_{bin_key}'))
            res_fail = ROOT.RooFitResult(ws.obj(f'results_fail_{bin_key}'))
            return (res_pass, res_fail)
        else:
            return (0, 0)
    
    if type_analysis == 'sim':
        if type(ws.obj(f'simPDF_({bin[0]}|{bin[1]})')) is ROOT.RooSimultaneous:
            print("Not possible to refit an existing PDF! \nReturning the results obtained previously")
            res = ROOT.RooFitResult(ws.obj(f'results_({bin[0]}|{bin[1]})'))
            return res
        else:
            return 0

###############################################################################

def pearson_chi2_eval(histo, pdf, nbins, res):
    """
    """
    ndof = nbins - res.floatParsFinal().getSize()
    # print(ndof)
    chi2_obj = ROOT.RooChi2Var("chi2", "chi2", pdf, histo)
    chi2_val = chi2_obj.getVal() 

    return chi2_val, ndof

###############################################################################

def check_chi2(histo, pdf, res):
    """
    """
    chi2val, ndof = pearson_chi2_eval(histo, pdf, histo.numEntries(), res)
    # print(chi2val, ndof)
    status_chi2 = bool(abs(chi2val - ndof) < 15*((2*ndof)**0.5))

    return status_chi2

###############################################################################

def fit_quality(res, type_checks="benchmark"):
    """
    """
    check_covm = (res.covQual() == 3)

    if type_checks == "benchmark":
        check_migrad = (res.status() == 0 or res.status() == 1)
        check_chi2 = (res.GetTitle() != "Chi2_not_passed")
        return bool(check_migrad*check_covm*check_chi2)

    elif type_checks == "bkg_fit":
        check_migrad = (res.status()==0) or (res.status()==3 and res.edm()<0.01)
        check_chi2 = (res.GetTitle() != "Chi2_not_passed")
        return bool(check_migrad*check_covm)
    
    elif type_checks == "new_checks":
        check_migrad = res.status() == 0
        check_edm = (res.edm() < 1e-4)
        return bool(check_migrad*check_covm*check_edm)

###############################################################################

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