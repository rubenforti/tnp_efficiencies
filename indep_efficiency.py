"""
"""

import ROOT
import time
import os
import sys
import copy
from utilities.results_utils import results_manager, efficiency_from_res
from utilities.plot_utils import plot_fitted_pass_fail
from utilities.base_library import eval_efficiency, sumw2_error
from utilities.fit_utils import fit_quality, check_chi2


def independent_efficiency(ws, bin_key, settings_dict, refit_numbkg=True, verb=-1, import_pdfs=False,
                           figs=False, figpath={"good":"figs/stuff", "check":"figs/check/stuff"}):
    """
    """
    # path = os.path.dirname(__file__)
    # ROOT.gSystem.cd(path)
    
    settings = copy.deepcopy(settings_dict)

    axis, histo_data = {}, {}
    smearing, pdf_mc, conv_pdf, bkg_pdf, sum_pdf, model_pdf = {}, {}, {}, {}, {}, {}
    results, fit_status = {}, {}

    bkg_tothist = {}
    bkg_pdf_ztautau, bkg_pdf_ttsemi = {}, {}
    h_pseudodata = {}
    
    if bin_key!="[24.0to26.0][-2.4to-2.3]":
        sys.exit()

    fit_pseudodata = True if ("fit_on_pseudodata" in settings.keys() and 
                              settings["fit_on_pseudodata"] is True) else False


    for flag in ["pass", "fail"]:

        axis.update({flag : ws.var(f"x_{flag}_{bin_key}")})
        axis[flag].setRange("fitRange", settings["fit_range"][0], settings["fit_range"][1])
        axis[flag].setBinning(ROOT.RooUniformBinning(axis[flag].getRange("fitRange")[0], 
                                                     axis[flag].getRange("fitRange")[1], 
                                                     settings["Nbins"][flag]), "cache")

        histo_data.update({flag : ws.data(f"Minv_data_{flag}_{bin_key}")})            

        # ---------------------------------------------------------------------------------------------------
        # --------------------- Parameters definition -------------------------------------------------------

        for pars_key in settings["pars"].keys():
            pars_obj = settings["pars"][pars_key][flag]
            settings["pars"][pars_key].update({
                flag : ROOT.RooRealVar(f"{pars_key}_{flag}_{bin_key}", f"{pars_key} {flag}", 
                                          pars_obj[0], pars_obj[1], pars_obj[2])})
        
        if fit_pseudodata is False:
            for norm_key in ["nsig", "nbkg"]:
                norm_obj = settings["norm"][norm_key][flag]
                for idx in range(3):
                    if (type(norm_obj[idx]) is str) and ("n" in norm_obj[idx]):
                        norm_obj[idx] = float(norm_obj[idx].replace("n",""))*histo_data[flag].sumEntries()
                    else:
                        norm_obj[idx] = float(norm_obj[idx])
                settings["norm"][norm_key].update({
                    flag : ROOT.RooRealVar(f"{norm_key}_{flag}_{bin_key}", f"{norm_key} {flag}",
                                            norm_obj[0], norm_obj[1], norm_obj[2])})


        # ---------------------------------------------------------------------------------------------------
        # -------------------- Signal PDF -------------------------------------------------------------------

        smearing.update({
            flag : ROOT.RooGaussian(f"smearing_{flag}_{bin_key}", "Gaussian smearing", 
                                    axis[flag], settings["pars"]["mu"][flag], settings["pars"]["sigma"][flag])})

        pdf_mc.update({
            flag : ROOT.RooHistPdf(f"pdf_mc_{flag}_{bin_key}", "pdf MC", 
                                   axis[flag], ws.data(f"Minv_mc_{flag}_{bin_key}"), 3)})
        
        conv_pdf.update({
            flag : ROOT.RooFFTConvPdf(f"conv_{flag}_{bin_key}", f"Convolution pdf", 
                                      axis[flag], pdf_mc[flag], smearing[flag])})
        conv_pdf[flag].setBufferFraction(settings["bufFraction"][flag])


        # ---------------------------------------------------------------------------------------------------
        # -------------------- Background PDF ---------------------------------------------------------------

        if settings["bkg_shape"][flag] == "expo":
            bkg_pdf.update({
                flag : ROOT.RooExponential(f"expo_bkg_{flag}_{bin_key}", "Exponential background",
                                           axis[flag], settings["pars"]["tau"][flag])})
            
        elif settings["bkg_shape"][flag] == "mixed":
            pass
        elif settings["bkg_shape"][flag] == "cmsshape":
            pass

        elif settings["bkg_shape"][flag] == "mc_raw":
            # histo_binning = axis.getBinning()
            bkg_tothist[flag] = ROOT.RooDataHist(f"Minv_bkg_{flag}_{bin_key}_total", "bkg_total_histo",
                                                 ROOT.RooArgSet(axis[flag]), "")
            for cat in settings["bkg_categories"]:
                bkg_tothist[flag].add(ws.data(f"Minv_bkg_{flag}_{bin_key}_{cat}")) 
            bkg_pdf.update({
                flag : ROOT.RooHistPdf(f"mcbkg_{flag}_{bin_key}", "MC-driven background", 
                                    ROOT.RooArgSet(axis[flag]), bkg_tothist[flag])})
            
        elif settings["bkg_shape"][flag] == "mc_double_pdf":
            n_ztautau = ws.data(f"Minv_bkg_{flag}_{bin_key}_Ztautau").sumEntries()
            n_ttsemileptonic = ws.data(f"Minv_bkg_{flag}_{bin_key}_TTSemileptonic").sumEntries()
            settings["norm"][flag].update({
                "ztautau": ROOT.RooRealVar(f"nbkg_ztautau_{flag}_{bin_key}", f"Nbkg {flag} Ztautau", 
                                           n_ztautau, 0.5, 1.5*n_ztautau)})
            settings["norm"][flag].update({
                "ttsemi": ROOT.RooRealVar(f"nbkg_ttsemileptonic_{flag}_{bin_key}", f"Nbkg {flag} TTSemi", 
                                          n_ttsemileptonic, 0.5, 1.5*n_ttsemileptonic)})
            settings["norm"][flag].update({
                "nbkg" : ROOT.RooAddition(f"nbkg_{flag}_{bin_key}", f"Nbkg {flag}", 
                                          ROOT.RooArgList(settings["norm"][flag]["ztautau"], 
                                                          settings["norm"][flag]["ttsemi"]))})
            bkg_pdf_ztautau.update({
                flag: ROOT.RooHistPdf(f"mcbkg_ztautau_{flag}_{bin_key}",
                                      "MC-driven bkg Ztautau", ROOT.RooArgSet(axis[flag]), 
                                       ws.data(f"Minv_bkg_{flag}_{bin_key}_Ztautau"), 3)})
            bkg_pdf_ttsemi.update({
                flag: ROOT.RooHistPdf(f"mcbkg_ttsemileptonic_{flag}_{bin_key}", 
                                      "MC-driven bkg TTSemi", ROOT.RooArgSet(axis[flag]), 
                                      ws.data(f"Minv_bkg_{flag}_{bin_key}_TTSemileptonic"), 3)})
            bkg_pdf.update({
                flag : ROOT.RooAddPdf(
                    f"mcbkg_{flag}_{bin_key}", "MC-driven background", 
                    ROOT.RooArgList(bkg_pdf_ztautau[flag], bkg_pdf_ttsemi[flag]), 
                    ROOT.RooArgList(settings["norm"]["ztautau"][flag], settings["norm"]["ttsemi"][flag]))})
        
        else:
            print("REQUESTED BACKGROUND SHAPE IS NOT SUPPORTED")
            sys.exit()


        # ---------------------------------------------------------------------------------------------------
        # -------------------- Pseudodata generation if requested -------------------------------------------

        if fit_pseudodata is True:
            del histo_data[flag]


            h_pseudodata[flag] = ROOT.RooDataHist(f"Minv_pseudodata_{flag}_{bin_key}", "pseudodata_histo",
                                                  ROOT.RooArgSet(axis[flag]), "")
            h_pseudodata[flag].add(ws.data(f"Minv_mc_{flag}_{bin_key}"))
            for cat in settings["bkg_categories"]:
                h_pseudodata[flag].add(ws.data(f"Minv_bkg_{flag}_{bin_key}_{cat}")) 
            histo_data.update({flag : h_pseudodata[flag]})
            
            
            if fit_pseudodata is False:
                for norm_key in ["nsig", "nbkg"]:
                    norm_obj = settings["norm"][norm_key][flag]
                    for idx in range(3):
                        if (type(norm_obj[idx]) is str) and ("n" in norm_obj[idx]):
                            norm_obj[idx] = float(norm_obj[idx].replace("n",""))*histo_data[flag].sumEntries()
                        else:
                            norm_obj[idx] = float(norm_obj[idx])
                    settings["norm"][norm_key].update({
                        flag : ROOT.RooRealVar(f"{norm_key}_{flag}_{bin_key}", f"{norm_key} {flag}",
                                                norm_obj[0], norm_obj[1], norm_obj[2])})


        # ---------------------------------------------------------------------------------------------------
        # -------------------- Final models and fits --------------------------------------------------------

        sum_pdf.update({
            flag : ROOT.RooAddPdf(f"sum_{flag}_{bin_key}", "Signal+Bkg", 
                                  ROOT.RooArgList(conv_pdf[flag], bkg_pdf[flag]), 
                                  ROOT.RooArgList(settings["norm"]["nsig"][flag], 
                                                  settings["norm"]["nbkg"][flag]))})
        sum_pdf[flag].setNormRange("fitRange")

        model_pdf.update({flag : ROOT.RooAddPdf(sum_pdf[flag], f'PDF_{flag}_{bin_key}')})
        model_pdf[flag].setNormRange("fitRange")

        res = model_pdf[flag].fitTo(histo_data[flag],
                                    # ROOT.RooFit.Extended(1),
                                    ROOT.RooFit.Range("fitRange"),
                                    ROOT.RooFit.Minimizer("Minuit2"),
                                    ROOT.RooFit.Strategy(settings["fit_strategy"][flag]),
                                    ROOT.RooFit.SumW2Error(False),
                                    # ROOT.RooFit.MaxCalls(100000),
                                    ROOT.RooFit.Save(1),
                                    ROOT.RooFit.PrintLevel(verb)
                                    )
        results.update({flag : res})

        if fit_pseudodata is False:
            status_chi2 = check_chi2(histo_data[flag], model_pdf[flag], results[flag], type="pearson")
        else:
            status_chi2 = True

        status = bool(status_chi2*fit_quality(results[flag], type_checks="benchmark"))

        print(status_chi2)
        print(status)
        
        nsig_fitted = settings["norm"]["nsig"][flag]
        nbkg_fitted = settings["norm"]["nbkg"][flag]
       
        low_nbkg = (nbkg_fitted.getVal() < 0.005*nsig_fitted.getVal())
        print(settings["norm"]["nbkg"][flag].getVal())
        print(settings["norm"]["nsig"][flag].getVal())
        print(low_nbkg)

        if (status is False) and refit_numbkg and low_nbkg:
            results[flag].Print()
            settings["pars"]["tau"][flag].setVal(1)
            settings["pars"]["tau"][flag].setConstant()


            if type(settings["norm"]["nbkg"][flag]) is ROOT.RooRealVar:
                    settings["norm"]["nbkg"][flag].setVal(0)
                    settings["norm"]["nbkg"][flag].setConstant()

            print("\n\nREFITTING WITHOUT BACKGROUND\n\n\n\n")

            res = model_pdf[flag].fitTo(histo_data[flag],
                                       # ROOT.RooFit.Extended(1),
                                       ROOT.RooFit.Range("fitRange"),
                                       ROOT.RooFit.Minimizer("Minuit2", "Migrad"),
                                       ROOT.RooFit.Strategy(settings["fit_strategy"][flag]),
                                       ROOT.RooFit.SumW2Error(False),
                                       # ROOT.RooFit.MaxCalls(100000),
                                       ROOT.RooFit.Save(1),
                                       ROOT.RooFit.PrintLevel(verb)
                                       )
            results.update({flag : res})
    
        fit_status[flag] = fit_quality(results[flag], type_checks="benchmark")

        if fit_pseudodata is False:
            fit_status[flag] = fit_status[flag]*check_chi2(histo_data[flag], model_pdf[flag], results[flag])

        results[flag].SetName(f"results_{flag}_{bin_key}")

    status = bool(fit_status["pass"]*fit_status["fail"])

    if status and import_pdfs:
        ws.Import(model_pdf["pass"]), ws.Import(model_pdf["fail"])
        ws.Import(results["pass"]), ws.Import(results["fail"])


    if figs:
        
        eff, d_eff = efficiency_from_res(results["pass"], results["fail"])

        eff_mc, deff_mc = eval_efficiency(ws.data(f"Minv_mc_pass_{bin_key}").sumEntries(), 
                                          ws.data(f"Minv_mc_fail_{bin_key}").sumEntries(),
                                          sumw2_error(ws.data(f"Minv_mc_pass_{bin_key}")),
                                          sumw2_error(ws.data(f"Minv_mc_fail_{bin_key}")))
        
        scale_factor = eff/eff_mc
        d_scale_factor = scale_factor*((d_eff/eff)**2 + (deff_mc/eff_mc)**2)**0.5

        plot_objects={}

        fig_folder = figpath["good"] if status is True else figpath["check"]

        for flag in ["pass", "fail"]:
            plot_objects.update({
                flag : {"axis" : axis[flag],
                        "data" : histo_data[flag],
                        "model" : model_pdf[flag],
                        "res" : results[flag]}
                        })
        plot_objects.update({"efficiency" : [eff, d_eff],
                             "efficiency_mc" : [eff_mc, deff_mc],
                             "scale_factor" : [scale_factor, d_scale_factor]})
        plot_fitted_pass_fail("indep", plot_objects, bin_key, pull=False, figpath=fig_folder)


    return results["pass"], results["fail"], status

      
