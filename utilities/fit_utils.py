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

def pearson_chi2_eval(histo, pdf):
    """
    """
    chi2_obj = ROOT.RooChi2Var("chi2", "chi2", pdf, histo)
    chi2_val = chi2_obj.getVal() 

    return chi2_val

###############################################################################

def llr_eval(histo, pdf, axis):
    """
    """

    binning = axis.getBinning()
    NBINS = binning.numBins()   
    EVTS = histo.sumEntries()

    sum_ll = 0
    max_ll = 0

    for i in range(NBINS):
    
        axis.setVal(binning.binCenter(i))
        weight = histo.weight(i)

        pdf_val = pdf.getVal(ROOT.RooArgSet(axis))
        mu = EVTS*(axis.getMax()-axis.getMin())*pdf_val/NBINS
        mu = round(mu, 5)

        new_sumll = weight*ROOT.TMath.Log(mu) - mu if mu>0 else 0.0
        sum_ll += 2*new_sumll

        if weight > 0:
            max_ll += 2*weight*ROOT.TMath.Log(weight) - 2*weight

    n_fitted_events=0
    for server in pdf.servers():
        if "nsig" in server.GetName() or "nbkg" in server.GetName():
            n_fitted_events += server.getVal()

    sum_ll += 2*EVTS*ROOT.TMath.Log(n_fitted_events) - 2*n_fitted_events
    max_ll += 2*EVTS*ROOT.TMath.Log(EVTS) - 2*EVTS

    llr = max_ll - sum_ll

    print("***")
    print("Nbins, Nevents, NfittedEvents:", NBINS, EVTS, n_fitted_events)
    print("SumLL, MaxLL", sum_ll, max_ll)
    print("***")

    return llr


###############################################################################

def status_chi2(axis, histo, pdf, res, type_chi2="pearson", nsigma=15):
    """
    """

    type_an = "indep" if "pass" in res.GetName() or "fail" in res.GetName() else None
    type_an = "sim" if "sim" in res.GetName() else None

    if type_an == "indep":
        flag = "pass" if "pass" in res.GetName() else "fail"
        chi2val = llr_eval(histo[flag], pdf[flag], axis[flag]) \
            if type_chi2=="llr" else pearson_chi2_eval(histo[flag], pdf[flag])
        ndof = histo[flag].numEntries() - res.floatParsFinal().getSize()
    elif type_an == "sim":
        chi2val = 0
        for flag in ["pass", "fail"]:
            if type_chi2=="llr": chi2val = chi2val + llr_eval(histo[flag], pdf[flag], axis[flag]) 
            else: chi2val = chi2val + pearson_chi2_eval(histo[flag], pdf[flag])
        ndof = histo["pass"].numEntries() + histo["fail"].numEntries() - res.floatParsFinal().getSize()

    else:
        print("ERROR: status_chi2() function is not implemented for this type of fit")
        sys.exit()
     
    chi2_status = bool(abs(chi2val - ndof) < nsigma*((2*ndof)**0.5))
    print(chi2val, ndof, chi2_status) 
    return chi2_status
    

###############################################################################

def fit_quality(fit_obj, type_checks="benchmark"):
    """
    """
    check_edm = True
    check_migrad = True
    check_chi2 = True
    check_covm = True

    if type_checks == "benchmark":
        check_migrad = (fit_obj["res"].status() == 0 or fit_obj["res"].status() == 1)
        check_covm = (fit_obj["res"].covQual() == 3)
        check_chi2 = status_chi2(fit_obj["axis"], fit_obj["histo"], fit_obj["pdf"],
                                 fit_obj["res"], type_chi2="pearson", nsigma=15)
    elif type_checks == "new_checks":
        check_migrad = (fit_obj["res"].status() == 0)
        check_covm = (fit_obj["res"].covQual() == 3)
        check_edm = (fit_obj["res"].edm() < 1e-3)
        check_chi2 = status_chi2(fit_obj["axis"], fit_obj["histo"], fit_obj["pdf"], 
                                 fit_obj["res"], type_chi2="llr", nsigma=5)
    elif type_checks == "pseudodata":
        check_migrad = (fit_obj["res"].status() == 0 or fit_obj["res"].status() == 1)
        check_covm = (fit_obj["res"].covQual() == 3)
    '''
    # Not used, could be useful if Sumw2Error turns out to be needed
    elif type_checks == "pseudodata_sumw2":
        check_migrad = fit_obj["res"].status()==0
        check_covm = (fit_obj["res"].covQual()==2 or fit_obj["res"].covQual()==3)
        return bool(check_migrad*check_covm)
    '''

    
    return bool(check_migrad*check_covm*check_edm*check_chi2)


    
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



if __name__ == "__main__":


    NBINS = 60

    axis = ROOT.RooRealVar("x_pass", "x", 60, 120)
    axis.setBins(NBINS)

    mu = ROOT.RooRealVar("mu", "mu", 91, 85, 95)
    sigma = ROOT.RooRealVar("sigma", "sigma", 2.5, 0.5, 5)

    gaus = ROOT.RooGaussian("gaus", "gaus", axis, mu, sigma)

    EVTS = int(NBINS*100)

    nevents = ROOT.RooRealVar("nsig", "numev", EVTS, 0, 5*EVTS)

    gaus_extended = ROOT.RooExtendPdf("gaus_extended", "gaus_extended", gaus, nevents)

    data = gaus_extended.generateBinned(ROOT.RooArgSet(axis), EVTS)

    res = gaus_extended.fitTo(data, 
                    ROOT.RooFit.Save(1), 
                    ROOT.RooFit.Minimizer("Minuit2", "Migrad"),
                    ROOT.RooFit.PrintLevel(-1))
    

    print(-2*res.minNll())
    
    llr, ndof = llr_eval(data, gaus_extended)

    print(llr, ndof)