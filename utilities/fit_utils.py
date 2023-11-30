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

def pearson_chi2_eval(histo, pdf, axis):
    """
    """
    binning = axis.getBinning("x_binning")
    NBINS = binning.numBins()   
    # EVTS = histo.sumEntries()
    BIN_VOLUME = (axis.getMax("fitRange") - axis.getMin("fitRange"))/NBINS

    PRECISION = 1

    n_fitted_events=0
    for server in pdf.servers():
        if "nsig" in server.GetName() or "nbkg" in server.GetName():
            n_fitted_events += server.getVal()

    chi2_val =0
    used_bins = 0


    for i in range(NBINS):
    
        axis.setVal(binning.binCenter(i))
        weight = histo.weight(i)

        pdf_val = pdf.getVal(ROOT.RooArgSet(axis))
        mu = n_fitted_events*BIN_VOLUME*pdf_val
        mu = round(mu, PRECISION)

        if mu > 0:
            chi2_val += (weight - mu)**2/mu
            used_bins += 1


    '''
    # Old version: can be problematic if the pdf assumes values near 0 but not exactly 0
    chi2_obj = ROOT.RooChi2Var("chi2", "chi2", pdf, histo,
                                ROOT.RooFit.Range("fitRange"),
                                ROOT.RooFit.DataError(ROOT.RooAbsData.Expected))
    chi2_val = chi2_obj.getVal() 
    '''
    return chi2_val, used_bins

###############################################################################

def llr_eval(histo, pdf, axis):
    """
    """

    binning = axis.getBinning()
    NBINS = binning.numBins()   
    EVTS = histo.sumEntries()
    BIN_VOLUME = (axis.getMax("fitRange") - axis.getMin("fitRange"))/NBINS

    sum_ll = 0
    max_ll = 0

    n_fitted_events=0
    for server in pdf.servers():
        if "nsig" in server.GetName() or "nbkg" in server.GetName():
            n_fitted_events += server.getVal()

    used_bins = 0
    for i in range(NBINS):
    
        axis.setVal(binning.binCenter(i))
        weight = histo.weight(i)

        pdf_val = pdf.getVal(ROOT.RooArgSet(axis))
        mu = n_fitted_events*BIN_VOLUME*pdf_val
        mu = round(mu, 5)

        new_sumll = weight*ROOT.TMath.Log(mu) - mu if mu>0 else 0.0
        sum_ll += 2*new_sumll

        if weight > 0:
            max_ll += 2*weight*ROOT.TMath.Log(weight) - 2*weight
            used_bins += 1

    
    # sum_ll += 2*EVTS*ROOT.TMath.Log(n_fitted_events) - 2*n_fitted_events
    # max_ll += 2*EVTS*ROOT.TMath.Log(EVTS) - 2*EVTS

    llr = max_ll - sum_ll

    print("***")
    print("Nbins, Nevents, NfittedEvents:", NBINS, EVTS, n_fitted_events)
    print("SumLL, MaxLL", sum_ll, max_ll)
    print("***")

    return llr, used_bins


###############################################################################

def status_chi2(axis, histo, pdf, res, type_chi2="pearson", nsigma=15):
    """
    """

    if ("pass" in res.GetName()) or ("fail" in res.GetName()):
        flag = "pass" if "pass" in res.GetName() else "fail"
        chi2val, used_bins = llr_eval(histo, pdf, axis) \
            if type_chi2 == "llr" else pearson_chi2_eval(histo, pdf, axis)
        ndof = used_bins - res.floatParsFinal().getSize()
        res.SetTitle(str(chi2val))
    
    elif "sim" in res.GetName():
        chi2val = 0
        used_bins = 0
        for flag in ["pass", "fail"]:
            if type_chi2=="llr": 
                chi2val += llr_eval(histo[flag], pdf[flag], axis[flag])[0]
                used_bins += llr_eval(histo[flag], pdf[flag], axis[flag])[1]
            else: 
                chi2val += pearson_chi2_eval(histo[flag], pdf[flag], axis[flag])[0]
                used_bins += pearson_chi2_eval(histo[flag], pdf[flag], axis[flag])[1]
        ndof = used_bins - res.floatParsFinal().getSize()
        res.SetTitle(str(chi2val))

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
    elif type_checks == "new":
        check_migrad = (fit_obj["res"].status() == 0)
        check_covm = (fit_obj["res"].covQual() == 3)
        check_edm = (fit_obj["res"].edm() < 1e-3)
        check_chi2 = status_chi2(fit_obj["axis"], fit_obj["histo"], fit_obj["pdf"], 
                                 fit_obj["res"], type_chi2="llr", nsigma=5)
    elif type_checks == "pseudodata":
        check_migrad = (fit_obj["res"].status() == 0)
        check_covm = (fit_obj["res"].covQual() == 3)
        check_edm = (fit_obj["res"].edm() < 1e-3)
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

    norm = ROOT.RooRealVar()
    tau = ROOT.RooRealVar()

    for par in pdf.getParameters(histo):
        print(par)
        if 'nbkg' in par.GetName(): norm = par
        if 'tau' in par.GetName(): tau = par


    poi_set = ROOT.RooArgSet("POI_set")
    poi_set.add(norm)
    poi_set.add(tau)

    null_par_set = ROOT.RooArgSet("null_params")

    null_par_set.addClone(norm)
    null_par_set.addClone(tau)
    null_par_set.setRealValue(norm.GetName(), 0)
    null_par_set.setRealValue(tau.GetName(), 2)

    llr_obj = ROOT.RooStats.ProfileLikelihoodCalculator(histo, pdf, ROOT.RooArgSet(norm, tau), 1-alpha, null_par_set)
    #llr_obj.SetConfidenceLevel(1-alpha)

    # Null hypotesis is that the background is 0.
    test_res = llr_obj.GetHypoTest()
    test_res.Print()
    pval = test_res.NullPValue()

    print(f"p-value for null hypo is: {pval}")

    # If the p-value is smaller than alpha, we reject the null hypothesis
    null = True if pval > alpha else False

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