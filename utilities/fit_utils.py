"""
"""

import sys
import os
import math
import ROOT
import json
import utilities.base_lib as base_lib
from array import array

# list of bkg strategies that use 
bkg_template_types = ["mc_raw", "BB_light", "num_estimation", "cmsshape_prefitBkg", "cmsshape_prefitBkg_SS"]

# dictionary that maps the bkg model to the corresponding custom C++ class to be loaded
dict_classes = {
        "cmsshape" : "RooCMSShape",
        "cmsshape_prefitBkg" : "RooCMSShape",
        "cmsshape_prefitBkg_SS" : "RooCMSShape",
        "cmsshape_new" : "RooCMSShape_mod",
        "CB" : "my_double_CB"
    }

fit_settings_args = ["type_analysis", "bkg_categories", "fitOnlyBkg", "fitPseudodata", "fit_verb", "refit_nobkg", "useMinos"]

###############################################################################

def base_parser_fit(parser):
    """
    """
    parser.add_argument("-n", "--process_name", default="", 
                        help="Name that characterizes the fit strategy")
    '''
    # The argument 'analysis' is already present in the base parser, better to use that one
    parser.add_argument("-f", "--type_fit", default="indep", choices=["indep", "sim", "sim_sf"],
                        help="Type of fit to be performed")
    '''
    parser.add_argument("-p", "--par_fit_settings", default="custom_fit_settings.json",
                        help="Name of the json file containing the fit settings (to be stored in the 'configs' folder)")

    parser.add_argument("--fitOnlyBkg", action="store_true",
                        help="Perform only the fit on background template")

    parser.add_argument("--fitPseudodata", action="store_true",
                        help="Perform the fit on pseudodata")

    parser.add_argument("--extended_sig_template_for_fail", action="store_true",
                        help="For failing probes, the signal template is built by using both pass and fail MC template, evaluated with SA variables")

    parser.add_argument("--auto_run", action="store_true",
                        help="Run the fits automatically, updating the fit settings according to the json file")

    parser.add_argument("--import_mc_SS", action="store_true",
                        help="Import the same-sign MC")

    parser.add_argument("--fit_verb", type=int, default=-1,
                        help="Verbosity of the fit output")

    parser.add_argument("--parallel", action="store_true",
                        help="Perform the fits in parallel")

    parser.add_argument("--refit_nobkg", action="store_true",
                        help="Refit the signal templates without the background")

    parser.add_argument("--useMinos", action="store_true",
                        help="Use Minos for the fit")

    parser.add_argument("--no_importPdfs", action="store_true",
                        help="Don't import the pdfs into the workspace")

    parser.add_argument("--no_saveFigs", action="store_true",
                        help="Don't save the fit plots")
    
    return parser


###############################################################################

def finalize_fit_parsing(args):
    """
    Apply control conditions on the arguments returned by the parser, regarding
    the fit settings
    """
    # Import of the fit settings
    with open(f"configs/{args.par_fit_settings.replace('.json', '')}.json") as file: fit_settings_json = json.load(file)
    if "legacy" in args.par_fit_settings:
        if args.eff in ["idip", "trigger", "iso"]:
            par_fit_settings = fit_settings_json["idip_trig_iso"]
        else:
            par_fit_settings = fit_settings_json[args.eff]
    else:
        if args.auto_run:
            par_fit_settings = fit_settings_json["run1"]
        else:
            nRUN = input("Insert the run number: ")
            par_fit_settings = fit_settings_json["run1"]
            [par_fit_settings.update(fit_settings_json["run"+str(run_idx)]) for run_idx in range(2, int(nRUN)+1)]
    args.par_fit_settings = par_fit_settings

    # Control on background categories
    if "all" in args.bkg_categories: args.bkg_categories = base_lib.bkg_categories

    # Apposite argument for import of background template
    args.importBkg = True if (args.fitPseudodata is True or args.fitOnlyBkg is True or \
                              args.par_fit_settings["bkg_model"]["pass"] in bkg_template_types or \
                              args.par_fit_settings["bkg_model"]["fail"] in bkg_template_types) else False 

    # Control on the mergedbins option
    if (args.mergedbins_bkg[0] != "" or args.mergedbins_bkg[1] != ""):
        if args.eff not in [ "idip", "trigger", "iso"]: 
            sys.exit("ERROR: mergedbins_bkg can be used only for idip, trigger, iso efficiencies")
        if args.binning_pt != "pt" or args.binning_eta != "eta":
            sys.exit("ERROR: Evaluation of background in merged bins for its comparison on data is allowed only w.r.t. standard bins of pt and eta for data")

    # Control on the output folder
    outpath = args.output if args.process_name in args.output else f"{args.output}/{args.process_name}"
    if not os.path.exists(outpath): os.makedirs(outpath)
    args.output = outpath

    # Setting the output folders for figures
    args.figpath = { "good"  : f"{args.output}/fit_plots",
                     "check" : f"{args.output}/fit_plots/check" }
    if "prefit" in args.par_fit_settings["bkg_model"]["pass"] or "prefit" in args.par_fit_settings["bkg_model"]["fail"]:
        args.figpath.update({ "prefit": f"{args.output}/fit_plots/prefit_bkg" })

    # Fit settings dictionary (to be passed to the 'fitter' object)
    print(vars(args))
    args.fit_settings = { key : vars(args)[key] for key in fit_settings_args }
    args.fit_settings.update(args.par_fit_settings)

###############################################################################


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


def checkImport_custom_pdf(bkg_models):
    """
    """
    dict_classes = {
        "cmsshape" : "RooCMSShape",
        "cmsshape_prefitBkg" : "RooCMSShape",
        "cmsshape_prefitBkg_SS" : "RooCMSShape",
        "cmsshape_new" : "RooCMSShape_mod",
        "CB" : "my_double_CB"
    }
    for flag in ["pass", "fail"]:
        if bkg_models[flag] in dict_classes.keys():
            base_lib.import_pdf_library(dict_classes[bkg_models[flag]])
            print(f"Imported {dict_classes[bkg_models[flag]]} from pdf_library")

###############################################################################


def printFitStatus(type_analysis, bin_key, res, status):

    if status is True: print(f"Bin {bin_key} is OK!\n\n")
    else:
        print(f"Bin {bin_key} has problems!\n")
        list_flags = ["pass", "fail"] if type_analysis == "indep" else ["sim"]
        print("****")
        for fl in list_flags:
            res[fl].Print()
            res[fl].correlationMatrix().Print()
            print("****")
        print('\n')
        print("Minimizer status")
        print("----------------")
        print("       ", *list_flags)
        print("----------------")
        print("Status ", *[res[fl].status() for fl in list_flags])
        print("CovQual", *[res[fl].covQual() for fl in list_flags])
        print("EDM    ", *[res[fl].edm() for fl in list_flags])
        print('\n')



###############################################################################


def getSidebands(histo, axis, cut=0.68):
    """
    """
    binning = axis.getBinning()
    NBINS = binning.numBins()

    print(binning)
    
    if NBINS%2 != 0:
        print("ERROR: number of bins must be even")
        sys.exit()

    ntot = histo.sumEntries()

    for i in range(int(NBINS/2)):

        idx_left, idx_right = int(NBINS/2-1-i), int(NBINS/2+i)

        lowLim, upLim = binning.binLow(idx_left), binning.binHigh(idx_right)
        ntot_sel = histo.sumEntries(f"x>{lowLim} && x<{upLim}")

        if ntot_sel/ntot > cut:
            break

    return lowLim, upLim

###############################################################################


def pearson_chi2_eval(histo, pdf, axis):
    """
    """
    print("Evaluating pearson chi2")
    binning = axis.getBinning("x_binning")
    NBINS = binning.numBins()   
    # EVTS = histo.sumEntries()
    BIN_VOLUME = (axis.getMax("fitRange") - axis.getMin("fitRange"))/NBINS
    PRECISION = 0

    n_fitted_events=0
    for server in pdf.servers():
        if "nsig" in server.GetName() or "nbkg" in server.GetName():
            n_fitted_events += server.getVal()
    if n_fitted_events == 0:
        print("Warning: not found normalization parameter in the pdf, using the number of events in the dataset as normalization")
        n_fitted_events = histo.sumEntries()

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

    flag = "pass" if "pass" in histo.GetName() else "fail"

    bin_key = histo.GetName().split("_")[-1]
    
    binning = axis.getBinning("x_binning")
    NBINS = binning.numBins()   
    EVTS = histo.sumEntries()
    BIN_VOLUME = (axis.getMax("fitRange") - axis.getMin("fitRange"))/NBINS

    isBBmodel = False

    sum_ll = 0
    max_ll = 0

    sum_nomi_bins = 0 #needed for the BB model

    n_fitted_events=0
    for serv_pdf in pdf.servers():
        if "nsig" in serv_pdf.GetName() or "nbkg" in serv_pdf.GetName():
            n_fitted_events += serv_pdf.getVal()

        if "BB" in serv_pdf.GetName():
            isBBmodel = True
            histConstr = serv_pdf.servers().findByName(f"bkg_histConstr_{flag}_{bin_key}")

            if histConstr.servers().size() != 2*NBINS: 
                return -1, 999
           
            for i in range(NBINS):
                sum_nomi_bins += histConstr.servers().findByName(f"bkg_histConstr_{flag}_{bin_key}_nominal_bin_{i}").getVal() * \
                                 histConstr.servers().findByName(f"bkg_paramHist_{flag}_{bin_key}_gamma_bin_{i}").getVal()

    if n_fitted_events == 0:
        print("Warning: not found normalization parameter in the pdf, using the number of events in the dataset as normalization")
        n_fitted_events = histo.sumEntries()

    used_bins = 0

    for i in range(NBINS):
    
        axis.setVal(binning.binCenter(i))
        weight = histo.weight(i)

        pdf_val = pdf.getVal(ROOT.RooArgSet(axis))

        if isBBmodel:

            nomi_bin = histConstr.servers().findByName(f"bkg_histConstr_{flag}_{bin_key}_nominal_bin_{i}").getVal()
            gamma_bin = histConstr.servers().findByName(f"bkg_paramHist_{flag}_{bin_key}_gamma_bin_{i}").getVal()
            bkg_val = gamma_bin*nomi_bin/sum_nomi_bins
            sig_val = pdf.servers().findByName(f"conv_{flag}_{bin_key}").getVal(ROOT.RooArgSet(axis))

            norm_bkg = pdf.servers().findByName(f"nbkg_{flag}_{bin_key}").getVal()
            norm_sig = pdf.servers().findByName(f"nsig_{flag}_{bin_key}").getVal()
            pdf_val = ((sig_val*norm_sig) + (bkg_val*norm_bkg))/n_fitted_events
            pdf_val = pdf_val

        mu = n_fitted_events*BIN_VOLUME*pdf_val
        mu = round(mu, 2)

        # if isBBmodel: print(sig_val*norm_sig, bkg_val*norm_bkg, mu)

        new_sumll = 2*round(weight*ROOT.TMath.Log(mu) - mu, 2) if mu>0 else 0.0
        sum_ll += new_sumll

        if weight > 0:
            new_maxll = 2*round(weight*ROOT.TMath.Log(weight) - weight, 2)
            max_ll += new_maxll
            used_bins += 1


        # print(round(weight,0), "", round(((sig_val*norm_sig) + bkg_val), 0), "", round(new_maxll-new_sumll, 1), "", round(max_ll - sum_ll, 1))

    
    # sum_ll += 2*EVTS*ROOT.TMath.Log(n_fitted_events) - 2*n_fitted_events
    # max_ll += 2*EVTS*ROOT.TMath.Log(EVTS) - 2*EVTS

    llr = max_ll - sum_ll

    print("***")
    print("Nbins, Nevents, NfittedEvents:", NBINS, EVTS, n_fitted_events)
    print("SumLL, MaxLL", sum_ll, max_ll)
    print("***")

    return llr, used_bins

###############################################################################


def status_parsAtLim(res, absTol=1e-4, relTol=1e-5):
    """
    """
    pars = res.floatParsFinal()
    status = True
    cnt_pars_at_lim = 0
    for par in pars:
        if par.isConstant(): continue
        final_val, par_min, par_max = par.getVal(), par.getMin(), par.getMax()

        if math.isclose(final_val, par_min, abs_tol=absTol, rel_tol=relTol) or \
           math.isclose(final_val, par_max, abs_tol=absTol, rel_tol=relTol):
            cnt_pars_at_lim += 1
            status = False
            print(par.GetName(), par.getVal(), par.getMin(), par.getMax())

            
    
    if status is False: print(f"WARNING: {res.GetName()} has {cnt_pars_at_lim} parameters at limit")

    return status

###############################################################################


def status_chi2(axis, histo, pdf, res, type_chi2="pearson", nsigma=15):
    """
    """

    if ("pass" in res.GetName()) or ("fail" in res.GetName()):
        flag = "pass" if "pass" in res.GetName() else "fail"
        chi2val, used_bins = llr_eval(histo, pdf, axis) \
            if type_chi2 == "llr" else pearson_chi2_eval(histo, pdf, axis)
    
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

    else:
        sys.exit("ERROR: status_chi2() function is not implemented for this type of fit")

    ndof = used_bins
    for par in res.floatParsFinal():
        if not ("_gamma_bin_" in par.GetName()): ndof -= 1  #To not count the gamma parameters in the BB method
    res.SetTitle(f"{chi2val}_{ndof}")

    print(chi2val, used_bins)
    
    chi2_status = bool(abs(chi2val - ndof) < nsigma*((2*ndof)**0.5))
    print(chi2val, ndof, chi2_status) 
    return chi2_status
    
###############################################################################


def fit_quality(fit_obj, type_checks="egm_legacy", isFitPseudodata=False, isBBmodel=False):
    """
    """
    check_edm = True
    check_migrad = True
    check_chi2 = True
    check_covm = True
    check_parsAtLim = True

    if type_checks == "egm_legacy":
        check_migrad = (fit_obj["res"].status()==0 or fit_obj["res"].status()==1)
        check_covm = (fit_obj["res"].covQual() == 3)
        if isBBmodel:
            check_chi2 = True
        else:
            check_chi2 = status_chi2(fit_obj["axis"], fit_obj["histo"], fit_obj["pdf"],
                                     fit_obj["res"], type_chi2="pearson", nsigma=10)
    
    elif type_checks == "new_tight":
        check_migrad = (fit_obj["res"].status() == 0)
        check_covm = (fit_obj["res"].covQual() == 3)
        check_edm = (fit_obj["res"].edm() < 2e-3)
        check_chi2 = status_chi2(fit_obj["axis"], fit_obj["histo"], fit_obj["pdf"], 
                                 fit_obj["res"], type_chi2="llr", nsigma=7)
        check_parsAtLim = status_parsAtLim(fit_obj["res"], absTol=1e-4, relTol=1e-5)
    
    elif type_checks == "new_loose":
        check_migrad = (fit_obj["res"].status()==0 or fit_obj["res"].status()==3)
        check_covm = (fit_obj["res"].covQual() == 3)
        check_edm = (fit_obj["res"].edm() < 5e-2)
        check_chi2 = status_chi2(fit_obj["axis"], fit_obj["histo"], fit_obj["pdf"], 
                                 fit_obj["res"], type_chi2="llr", nsigma=10)
        check_parsAtLim = status_parsAtLim(fit_obj["res"], absTol=1e-4, relTol=1e-5)
    
    elif type_checks == "neglect_atLimPars":
        check_migrad = (fit_obj["res"].status()==0)
        check_covm = (fit_obj["res"].covQual() == 3)
        check_edm = (fit_obj["res"].edm() < 1e-2)
        check_chi2 = status_chi2(fit_obj["axis"], fit_obj["histo"], fit_obj["pdf"],
                                 fit_obj["res"], type_chi2="llr", nsigma=16)

    elif type_checks == "pseudodata":
        check_migrad = (fit_obj["res"].status() == 0)
        check_covm = (fit_obj["res"].covQual() == 3)
        check_edm = (fit_obj["res"].edm() < 5e-2)
        check_parsAtLim = status_parsAtLim(fit_obj["res"], absTol=1e-4, relTol=1e-5)
    
    else:
        sys.exit("ERROR: wrong type of fit quality check indicated")


    print("check chi2", check_chi2)

    if isFitPseudodata: check_chi2 = True

    return bool(check_migrad*check_covm*check_edm*check_chi2*check_parsAtLim)
 
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

###############################################################################


# Useful for parallel fits

def doSingleFit(fitter, ws, flags):
    """
    """
    fitter.manageFit(ws)

    res = {flag : fitter.res_obj[flag] for flag in flags}

    return fitter, res


###############################################################################
###############################################################################


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
