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
            res_pass = ws.obj(f'results_pass_{bin_key}')
            res_fail = ws.obj(f'results_fail_{bin_key}')
            return (res_pass, res_fail)
        else:
            return 0
    
    if type_analysis == 'sim':
        if type(ws.obj(f'PDF_{bin_key}')) is ROOT.RooSimultaneous:
            print("Not possible to refit an existing PDF! \nReturning the results obtained previously")
            res = ws.obj(f'results_{bin_key}')
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

def llr_eval(histo, res):
    """
    """
    max_ll_data = 0

    for idx in range(histo.numEntries()):
        k = histo.weight(idx)
        if k == 0:
            max_ll_data += 0
        elif k>0:
            max_ll_data += 2*(k*ROOT.TMath.Log(k))

    llr_val = max_ll_data #  + 2*res.minNll()

    pars = res.floatParsFinal()
    ndof = histo.numEntries() - pars.getSize()

    return llr_val, ndof



###############################################################################

def status_chi2(histo, pdf, res, type="pearson", nsigma=15):
    """
    """

    if "pseudodata" in histo.GetName():
        return True

    else:
        if type == "pearson":
            chi2val, ndof = pearson_chi2_eval(histo, pdf, histo.numEntries(), res)
            print(chi2val, ndof)
        elif type == "llr":
            chi2val, ndof = llr_eval(histo, res)
        
        chi2_status = bool(abs(chi2val - ndof) < nsigma*((2*ndof)**0.5))

        return chi2_status

###############################################################################

def fit_quality(fit_obj, type_checks="benchmark"):
    """
    """
    check_covm = (fit_obj["res"].covQual() == 3)

    if type_checks == "benchmark":
        check_migrad = (fit_obj["res"].status() == 0 or fit_obj["res"].status() == 1)
        check_chi2 = status_chi2(fit_obj["histo"], fit_obj["pdf"], fit_obj["res"], type="pearson")
        return bool(check_covm*check_migrad*check_chi2)
    
    elif type_checks == "new_checks":
        check_migrad = fit_obj["res"].status() == 0
        check_edm = (fit_obj["res"].edm() < 1e-3)
        check_chi2 = status_chi2(fit_obj["histo"], fit_obj["pdf"], fit_obj["res"], 
                                 type="llr", nsigma=5)
        return bool(check_covm*check_migrad*check_edm*check_chi2)
    '''
    # Not used anymore, could be useful if Sumw2Error turns out to be needed
    elif type_checks == "pseudodata":
        check_migrad = fit_obj["res"].status()==0
        check_covm = (fit_obj["res"].covQual()==2 or fit_obj["res"].covQual()==3)
        return bool(check_migrad*check_covm)
    '''

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
