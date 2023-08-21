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
from utilities.fit_utils import fit_quality


class EfficiencyFitter():
    """
    """

    def __init__(self, type_analysis, bin_key, settings):
        """
        """
        self.__type_analysis = type_analysis
        self.__settings = copy.deepcopy(settings)
        self.__bin_key = bin_key
        
        self.__axis = {}
        self.__smearing, self.__pdf_mc, self.__conv_pdf = {}, {}, {}
        self.__bkg_pdf, self.__sum_pdf, self.__model_pdf = {}, {}, {}
        self.__histo_data = {}

        self.__fit_pseudodata = True if ("fit_on_pseudodata" in settings.keys() 
                                         and self.__settings["fit_on_pseudodata"] is True) else False

        if self.__settings["bkg_model"] == "mc_raw":
            self.__bkg_categories = self.__settings["bkg_categories"]
            self.__bkg_tothist = {}

    def importAxis(self, flag, ws):
        self.__axis.update({flag : ws.var(f"x_{flag}_{self.__bin_key}")})
        self.__axis[flag].setRange("fitRange", 
                                   self.__settings["fit_range"][0], self.__settings["fit_range"][1])
        self.__axis[flag].setBinning(ROOT.RooUniformBinning(self.__axis[flag].getRange("fitRange")[0], 
                                                            self.__axis[flag].getRange("fitRange")[1], 
                                                            self.__settings["Nbins"][flag]), "cache")

    def initParams(self, flag):
        for pars_key in self.__settings["pars"].keys():
            pars_obj = self.__settings["pars"][pars_key][flag]
            self.__settings["pars"][pars_key].update({
                flag : ROOT.RooRealVar(f"{pars_key}_{flag}_{self.__bin_key}", f"{pars_key} {flag}", 
                                       pars_obj[0], pars_obj[1], pars_obj[2])})

    def initNorm_indep(self, flag):
        for norm_key in ["nsig", "nbkg"]:
                norm_obj = self.__settings["norm"][norm_key][flag]
                for idx in range(3):
                    if (type(norm_obj[idx]) is str) and ("n" in norm_obj[idx]):
                        n_events = self.__histo_data[flag].sumEntries()
                        norm_obj[idx] = float(norm_obj[idx].replace("n",""))*n_events
                    else:
                        norm_obj[idx] = float(norm_obj[idx])
                self.__settings["norm"][norm_key].update({
                    flag : ROOT.RooRealVar(f"{norm_key}_{flag}_{self.__bin_key}", f"{norm_key} {flag}",
                                           norm_obj[0], norm_obj[1], norm_obj[2])})


    def createSigPdf(self, flag, data_mc):
        self.__smearing.update({ flag : ROOT.RooGaussian(
            f"smearing_{flag}_{self.__bin_key}", "Gaussian smearing", self.__axis[flag], 
            self.__settings["pars"]["mu"][flag], self.__settings["pars"]["sigma"][flag]) })

        self.__pdf_mc.update({ flag : ROOT.RooHistPdf(f"pdf_mc_{flag}_{self.__bin_key}", "pdf MC",
                                                      self.__axis[flag], data_mc, 3) })
        
        self.__conv_pdf.update({ flag : ROOT.RooFFTConvPdf(
            f"conv_{flag}_{self.__bin_key}", f"Convolution pdf", self.__axis[flag], 
            self.__pdf_mc[flag], self.__smearing[flag]) })
        self.__conv_pdf[flag].setBufferFraction(self.__settings["bufFraction"][flag])
        self.__conv_pdf[flag].setNormRange("fitRange")

    def createBkgPdf(self, flag, ws=""):
        bkg_model = self.__settings["bkg_model"][flag]
        if bkg_model == "expo":
            self.__bkg_pdf.update({ flag : ROOT.RooExponential(
                f"expo_bkg_{flag}_{self.__bin_key}","Exponential background",
                self.__axis[flag], self.__settings["pars"]["tau"][flag]) })
            self.__bkg_pdf[flag].setNormRange("fitRange")       
            
        elif bkg_model == "mixed":
            pass
        elif bkg_model == "cmsshape":
            pass

        elif bkg_model == "mc_raw":
            # histo_binning = axis.getBinning()
            self.__bkg_tothist[flag] = ROOT.RooDataHist(
                f"Minv_bkg_{flag}_{self.__bin_key}_total", "bkg_total_histo",
                ROOT.RooArgSet(self.__axis[flag]), "")

            for cat in self.__bkg_categories:
                self.__bkg_tothist[flag].add(ws.data(f"Minv_bkg_{flag}_{self.__bin_key}_{cat}")) 
    
            self.__bkg_pdf.update({ flag : ROOT.RooHistPdf(
                f"mcbkg_{flag}_{self.__bin_key}", "MC-driven background", 
                ROOT.RooArgSet(self.__axis[flag]), self.__bkg_tothist[flag]) })
            

    def create_pseudodata(self, flag, ws):
        roohist_pseudodata = ROOT.RooDataHist(
            f"Minv_pseudodata_{flag}_{self.__bin_key}", "pseudodata_histo",
            ROOT.RooArgSet(self.__axis[flag]), "")

        roohist_pseudodata.add(ws.data(f"Minv_mc_{flag}_{self.__bin_key}"))

        for cat in self.__settings["bkg_categories"]:
            roohist_pseudodata.add(ws.data(f"Minv_bkg_{flag}_{self.__bin_key}_{cat}")) 

        self.__histo_data.update({flag : roohist_pseudodata})


    def doFit(self, flag, ws, init_pdfs=True):
        """
        """
        if init_pdfs:
            self.importAxis(flag, ws)
            self.initParams(flag)
            self.createSigPdf(flag, ws.data(f"Minv_mc_{flag}_{self.__bin_key}"))
            self.createBkgPdf(flag, ws)

            if len(self.__histo_data.keys())<2:
                if self.__fit_pseudodata is False:
                    self.__histo_data.update({flag : ws.data(f"Minv_data_{flag}_{self.__bin_key}")})
                else:
                    self.create_pseudodata(flag, ws)
        
            self.initNorm_indep(flag)

            self.__sum_pdf.update({
                flag : ROOT.RooAddPdf(f"sum_{flag}_{self.__bin_key}", "Signal+Bkg", 
                                    ROOT.RooArgList(self.__conv_pdf[flag], self.__bkg_pdf[flag]), 
                                    ROOT.RooArgList(self.__settings["norm"]["nsig"][flag], 
                                                    self.__settings["norm"]["nbkg"][flag]))})
            self.__sum_pdf[flag].setNormRange("fitRange")

            self.__model_pdf.update({flag : ROOT.RooAddPdf(self.__sum_pdf[flag], 
                                                        f'PDF_{flag}_{self.__bin_key}')})
            self.__model_pdf[flag].setNormRange("fitRange")

        res = self.__model_pdf[flag].fitTo(self.__histo_data[flag],
                                           # ROOT.RooFit.Extended(1),
                                           ROOT.RooFit.Range("fitRange"),
                                           ROOT.RooFit.Minimizer("Minuit2"),
                                           ROOT.RooFit.Strategy(self.__settings["fit_strategy"][flag]),
                                           ROOT.RooFit.SumW2Error(False),
                                           # ROOT.RooFit.MaxCalls(100000),
                                           ROOT.RooFit.Save(1),
                                           ROOT.RooFit.PrintLevel(-1)
                                           )
        res.SetName(f"results_{flag}_{self.__bin_key}")
        
        fit_obj = {"histo" : self.__histo_data[flag], "pdf" : self.__model_pdf[flag], "res" : res}
        status = fit_quality(fit_obj, type_checks=self.__settings["fit_checks"])

        return res, status
    

    def refit_noBkg(self, flag, ws):
        """
        """
        self.__settings["pars"]["tau"][flag].setVal(1)
        self.__settings["pars"]["tau"][flag].setConstant()
        
        if type(self.__settings["norm"]["nbkg"][flag]) is ROOT.RooRealVar:
                self.__settings["norm"]["nbkg"][flag].setVal(0)
                self.__settings["norm"]["nbkg"][flag].setConstant()
        
        return self.doFit(flag, ws, init_pdfs=False)
    

    def getFinalPdf(self, flag):
        """
        """
        return self.__model_pdf[flag]


    def saveFig(self, ws, results, status, figpath):
        """
        """
        eff, d_eff = efficiency_from_res(results["pass"], results["fail"])

        eff_mc, deff_mc = eval_efficiency(ws.data(f"Minv_mc_pass_{self.__bin_key}").sumEntries(), 
                                          ws.data(f"Minv_mc_fail_{self.__bin_key}").sumEntries(),
                                          sumw2_error(ws.data(f"Minv_mc_pass_{self.__bin_key}")),
                                          sumw2_error(ws.data(f"Minv_mc_fail_{self.__bin_key}")))
        
        scale_factor = eff/eff_mc
        d_scale_factor = scale_factor*((d_eff/eff)**2 + (deff_mc/eff_mc)**2)**0.5

        plot_objects={}
        fig_folder = figpath["good"] if status is True else figpath["check"]
        for flag in ["pass", "fail"]:
            plot_objects.update({
                flag : {"axis" : self.__axis[flag],
                        "data" : self.__histo_data[flag],
                        "model" : self.__model_pdf[flag],
                        "res" : results[flag]}
                        })
        plot_objects.update({"efficiency" : [eff, d_eff],
                             "efficiency_mc" : [eff_mc, deff_mc],
                             "scale_factor" : [scale_factor, d_scale_factor]})
        plot_fitted_pass_fail("indep", plot_objects, self.__bin_key, pull=False, figpath=fig_folder)
        




'''
def independent_efficiency(ws, bin_key, settings_dict, refit_numbkg=True, verb=-1, import_pdfs=False,
                           figs=False, figpath={"good":"figs/stuff", "check":"figs/check/stuff"}):
    """
    """
    
    if bin_key!="[24.0to26.0][-2.4to-2.3]":
        sys.exit()

    results, status_dict = {}, {}            
   
    fitter = EfficiencyFitter("indep", bin_key, settings_dict)


    for flag in ["pass", "fail"]:
        res, status = fitter.doFit(flag, ws, type_checks = "benchmark")
        
        pars_fitted = res.floatParsFinal()
        nsig_fitted = pars_fitted.find(f"nsig_{flag}_{bin_key}")
        nbkg_fitted = pars_fitted.find(f"nbkg_{flag}_{bin_key}")
       
        low_nbkg = (nbkg_fitted.getVal() < 0.005*nsig_fitted.getVal())

        if (status is False) and refit_numbkg and low_nbkg:
            res, status = fitter.refit_noBkg(flag, ws)

        results.update({flag : res})
        status_dict.update({flag : status})

    status = bool(status_dict["pass"]*status_dict["fail"])

    if status and import_pdfs:
        ws.Import(fitter.getFinalPdf("pass")), ws.Import(fitter.getFinalPdf("fail"))
        ws.Import(results["pass"]), ws.Import(results["fail"])
    

    if figs:
        fitter.saveFig(ws, results, status, figpath)
        
        
    return results["pass"], results["fail"], False
'''