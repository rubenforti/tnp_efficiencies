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


class AbsFitter():

    """
    Parent class that initializes the objects for the efficiency fits and
    contains basic methods. Needs to be inherited by the child classes that 
    actually perform the fits.
    """

    def __init__(self, bin_key, fit_settings):
        """
        Constructor. Initializes the basic attributes of the fitter.
        """
        self.bin_key = bin_key
        settings_obj = copy.deepcopy(fit_settings)
        self.pars = settings_obj.pop("pars")
        self.norm = settings_obj.pop("norm")
        self.settings = settings_obj
        
        self.histo_data = {}
        self.axis = {}
        self.pdfs = {"smearing" : {}, "mc_pdf" : {}, "conv_pdf" : {}, 
                     "bkg_pdf" : {}, "sum_pdf" : {}, "fit_model" : {} }
        self.results = {}


        if self.settings["bkg_model"] == "mc_raw":
            self.bkg_categories = self.settings["bkg_categories"]
            self.bkg_tothist = {}


    def _initParams(self, flag):
        """
        Transforms the list of numbers, representing each parameter (initial
        value, lower bound, upper bound), into a RooRealVar object. 
        """
        for pars_key in self.pars.keys():
            pars_obj = self.pars[pars_key][flag]
            self.pars[pars_key].update({ flag : ROOT.RooRealVar(
                f"{pars_key}_{flag}_{self.bin_key}", f"{pars_key} {flag}", 
                pars_obj[0], pars_obj[1], pars_obj[2]) })


    def _createSigPdf(self, flag, ws):
        """
        Creates the signal pdf and loads it into the pdfs dictionary. The lower
        level pdfs (i.e. the smearing and MC pds) are loaded as well.
        """
        self.pdfs["smearing"].update({ flag : ROOT.RooGaussian(
            f"smearing_{flag}_{self.bin_key}", "Gaussian smearing", 
            self.axis[flag], self.pars["mu"][flag], self.pars["sigma"][flag]) })

        self.pdfs["mc_pdf"].update({ flag : ROOT.RooHistPdf(
            f"mc_pdf_{flag}_{self.bin_key}", "MC pdf",
            self.axis[flag], ws.data(f"Minv_mc_{flag}_{self.bin_key}"), 3) })
        
        if self.settings["type_analysis"] == "sim_sf":
            self.pdfs["mc_pdf"].update({ flag : ROOT.RooExtendPdf(
                f"conv_{flag}_{self.bin_key}", "pdf MC", self.pdfs["mc_pdf"][flag], 
                ws.data(f"Minv_mc_{flag}_{self.bin_key}").sumEntries()) })
        
        self.pdfs["conv_pdf"].update({ flag : ROOT.RooFFTConvPdf(
            f"conv_{flag}_{self.bin_key}", f"Convolution pdf", self.axis[flag], 
            self.pdfs["mc_pdf"][flag], self.pdfs["smearing"][flag]) })
        self.pdfs["conv_pdf"][flag].setBufferFraction(self.settings["bufFraction"][flag])
        self.pdfs["conv_pdf"][flag].setNormRange("fitRange")


    def _createBkgPdf(self, flag, ws):
        """
        Creates the background pdf and loads it into the pdfs dictionary. In 
        case the background is modeled with MC, the histogram of the sum of the
        backgrounds is created and loaded in the apposite attribute.
        """
        bkg_model = self.settings["bkg_model"][flag]
        if bkg_model == "expo":
            self.pdfs["bkg_pdf"].update({ flag : ROOT.RooExponential(
                f"expo_bkg_{flag}_{self.bin_key}","Exponential background",
                self.axis[flag], self.pars["tau"][flag]) })
            self.pdfs["bkg_pdf"][flag].setNormRange("fitRange")       
            
        elif bkg_model == "mixed":
            pass
        elif bkg_model == "cmsshape":
            pass

        elif bkg_model == "mc_raw":
            # histo_binning = axis.getBinning()
            self.bkg_tothist[flag] = ROOT.RooDataHist(
                f"Minv_bkg_{flag}_{self.bin_key}_total", "bkg_total_histo",
                ROOT.RooArgSet(self.axis[flag]), "")

            for cat in self.bkg_categories:
                self.bkg_tothist[flag].add(ws.data(f"Minv_bkg_{flag}_{self.bin_key}_{cat}")) 
    
            self.pdfs["bkg_pdf"].update({ flag : ROOT.RooHistPdf(
                f"mcbkg_{flag}_{self.bin_key}", "MC-driven background", 
                ROOT.RooArgSet(self.axis[flag]), self.bkg_tothist[flag]) })
            

    def _createPseudodata(self, flag, ws):
        """
        Creates the pseudodata histogram, which is the sum of the signal and
        background MC histograms. The histogram is loaded into the histo_data
        attribute.
        """
        roohist_pseudodata = ROOT.RooDataHist(
            f"Minv_pseudodata_{flag}_{self.bin_key}", "pseudodata_histo",
            ROOT.RooArgSet(self.axis[flag]), "")

        roohist_pseudodata.add(ws.data(f"Minv_mc_{flag}_{self.bin_key}"))

        [roohist_pseudodata.add(ws.data(f"Minv_bkg_{flag}_{self.bin_key}_{cat}"))
         for cat in self.settings["bkg_categories"]]

        self.histo_data.update({flag : roohist_pseudodata})

    def checkExistingFit(self, ws):
        """
        Checks if the fit has already been performed, by looking for the 
        results object in the workspace. If so, it loads the results object
        into the results attribute.
        """

        if self.settings["type_analysis"] == 'indep':
            isFittedPass = type(ws.obj(f'fitPDF_pass_{self.bin_key}')) is ROOT.RooAddPdf
            isFittedFail = type(ws.obj(f'fitPDF_fail_{self.bin_key}')) is ROOT.RooAddPdf
            if isFittedPass and isFittedFail:
                print("Not possible to refit an existing PDF! \nReturning the results obtained previously")
                self.results = {
                    "pass": {"res_obj": ws.obj(f'results_pass_{self.bin_key}'), "status" : True},
                    "fail": {"res_obj": ws.obj(f'results_fail_{self.bin_key}'), "status" : True}
                }
                self.existingFit = True
            else:
                self.existingFit = False
    
        elif self.settings["type_analysis"] == 'sim':
            if type(ws.obj(f'fitPDF_sim_{self.bin_key}')) is ROOT.RooSimultaneous:
                print("Not possible to refit an existing PDF! \nReturning the results obtained previously")
                self.results = { "sim" : { "res_obj" : ws.obj(f'results_sim_{self.bin_key}'), "status" : True } }
                self.existingFit = True
            else:
                self.existingFit = False
   
    def importAxis(self, flag, ws):
        """
        Imports the axis from the workspace. If the simultaneous analysis is 
        being performed, the same object is referenced for both the pass and
        fail categories in the dictionary for the sake of simplicity; there 
        are no problems in doing so, since the object in memory is the same. 
        """
        if self.settings["type_analysis"] == "sim" or self.settings["type_analysis"] == "sim_sf":
            flag_ws = "sim" 
        else:
            flag_ws = flag
        self.axis.update({flag : ws.var(f"x_{flag_ws}_{self.bin_key}")})
        self.axis[flag].setRange("fitRange", 
                                   self.settings["fit_range"][0], self.settings["fit_range"][1])
        self.axis[flag].setBinning(ROOT.RooUniformBinning(self.axis[flag].getRange("fitRange")[0], 
                                                          self.axis[flag].getRange("fitRange")[1], 
                                                          self.settings["Nbins"][flag]), "cache")

    def initDatasets(self, ws, import_mc=False):
        """
        """
        for flag in ["pass", "fail"]:
            if self.settings["fit_on_pseudodata"] is False:
                self.histo_data.update({flag : ws.data(f"Minv_data_{flag}_{self.bin_key}")})
                if import_mc is True:
                    self.histo_data.update({f"mc_{flag}" : ws.data(f"Minv_mc_{flag}_{self.bin_key}")})
            else:
                self._createPseudodata(flag, ws)
    

    def initNorm(self):
        """
        Transforms the normalization parameters, for signal and background,
        into RooRealVar objects, with the same strategy of initParams() method.
        If a value is given as a string, with an "n" at the end, it is 
        interpreted as the number of events in the dataset times that number.
        """
        for flag in ["pass", "fail"]:
            for norm_key in ["nsig", "nbkg"]:
                norm_obj = self.norm[norm_key][flag]
                for idx in range(3):
                    if (type(norm_obj[idx]) is str) and ("n" in norm_obj[idx]):
                        n_events = self.histo_data[flag].sumEntries()
                        norm_obj[idx] = float(norm_obj[idx].replace("n",""))*n_events
                    else:
                        norm_obj[idx] = float(norm_obj[idx])
                self.norm[norm_key].update({ flag : ROOT.RooRealVar(
                    f"{norm_key}_{flag}_{self.bin_key}", f"{norm_key} {flag}",
                    norm_obj[0], norm_obj[1], norm_obj[2]) })


    def initFitPdfs(self, ws):
        """
        """
        for flag in ["pass", "fail"]:
            self._initParams(flag)
            self._createSigPdf(flag, ws)
            self._createBkgPdf(flag, ws)
            self.pdfs["sum_pdf"].update({ flag : ROOT.RooAddPdf(
                        f"sum_{flag}_{self.bin_key}", "Signal+Bkg model", 
                        ROOT.RooArgList(self.pdfs["conv_pdf"][flag], self.pdfs["bkg_pdf"][flag]), 
                        ROOT.RooArgList(self.norm["nsig"][flag], self.norm["nbkg"][flag])) })
            self.pdfs["sum_pdf"][flag].setNormRange("fitRange")

            self.pdfs["fit_model"].update({ flag : ROOT.RooAddPdf(
                self.pdfs["sum_pdf"][flag], f'fitPDF_{flag}_{self.bin_key}')})
            self.pdfs["fit_model"][flag].setNormRange("fitRange")
    

    def doFit(self, flag):
        """
        """
        if self.settings["type_analysis"] == "sim":
            if self.settings["useMinos"] is False:
                par_set = "none"
            elif self.settings["useMinos"] == "eff":
                par_set = ROOT.RooArgSet()
                par_set.add(self.efficiency)
            elif self.settings["useMinos"] == "all":
                par_set = self.pdfs["fit_model"][flag].getParameters(self.histo_data[flag])

        res = self.pdfs["fit_model"][flag].fitTo(self.histo_data[flag],
                                                  ROOT.RooFit.Range("fitRange"),
                                                  ROOT.RooFit.Minimizer("Minuit2", "Migrad"),
                                                  ROOT.RooFit.Minos(par_set),
                                                  ROOT.RooFit.SumW2Error(False),
                                                  ROOT.RooFit.Save(1), 
                                                  ROOT.RooFit.PrintLevel(self.settings["fit_verb"]))
        res.SetName(f"results_{flag}_{self.bin_key}")
        
        res.Print("v")

        fit_obj = {"axis" : self.axis, "histo" : self.histo_data, 
                   "pdf" : self.pdfs["fit_model"], "res" : res}
    
        self.results.update({ flag : {
            "res_obj" : res, 
            "status" : fit_quality(fit_obj, type_checks=self.settings["fit_checks"])} })


    def saveFig(self, ws, figpath):
        """
        """
        type_an = self.settings["type_analysis"]

        if self.existingFit is False:

            if type_an == "indep":
                eff, d_eff = efficiency_from_res(self.results["pass"]["res_obj"], 
                                                self.results["fail"]["res_obj"])
                fig_status = bool(self.results["pass"]["status"]*self.results["fail"]["status"])
            if type_an == "sim":
                eff, d_eff = self.efficiency.getVal(), self.efficiency.getError()
                fig_status = bool(self.results["sim"]["status"])
                self.results.update({"pass" : self.results["sim"]["res_obj"], 
                                     "fail" : self.results["sim"]["res_obj"],
                                     "sim" : self.results["sim"]["res_obj"]})
            else:
                pass
            # ATTENZIONE: CONTROLLARE CHE SUCCEDE SE SI FA IL FIT CON L'OPZIONE "SIM_SF"

            eff_mc, deff_mc = eval_efficiency(ws.data(f"Minv_mc_pass_{self.bin_key}").sumEntries(), 
                                            ws.data(f"Minv_mc_fail_{self.bin_key}").sumEntries(),
                                            sumw2_error(ws.data(f"Minv_mc_pass_{self.bin_key}")),
                                            sumw2_error(ws.data(f"Minv_mc_fail_{self.bin_key}")))
            
            scale_factor = eff/eff_mc
            d_scale_factor = scale_factor*((d_eff/eff)**2 + (deff_mc/eff_mc)**2)**0.5

            print(self.results)
            plot_objects={}
            fig_folder = figpath["good"] if fig_status is True else figpath["check"]
            for flag in ["pass", "fail"]:
                plot_objects.update({ flag : {
                            "axis" : self.axis[flag],
                            "data" : self.histo_data[flag],
                            "model" : self.pdfs["fit_model"][flag],
                            "res" : self.results[flag],
                            } })
            print(plot_objects["pass"]["res"])
            plot_objects.update({"efficiency" : [eff, d_eff],
                                "efficiency_mc" : [eff_mc, deff_mc],
                                "scale_factor" : [scale_factor, d_scale_factor]})
            plot_fitted_pass_fail(type_an, plot_objects, self.bin_key, pull=False, figpath=fig_folder)
        else:
            pass


###############################################################################
###############################################################################


class IndepFitter(AbsFitter):
    """
    """
    def __init__(self, bin_key, settings):
        """
        """
        AbsFitter.__init__(self, bin_key, settings)


    def attempt_noBkgFit(self, flag):
        """
        """
        if self.settings["bkg_model"][flag] != "expo":
            # Need to be implemented for, e.g., cmsshape
            pass
        else: 
            nsig_fitted, nbkg_fitted = self.norm["nsig"][flag], self.norm["nbkg"][flag]
            if nbkg_fitted.getVal() < 0.005*nsig_fitted.getVal():
                self.pars["tau"][flag].setVal(1)
                self.pars["tau"][flag].setConstant()
                self.norm["nbkg"][flag].setVal(0)
                self.norm["nbkg"][flag].setConstant()
                self.doFit(flag)
            else:
                pass
    

    def manageFit(self, ws):
        """
        """
        self.checkExistingFit(ws)

        if self.existingFit is False:
            [self.importAxis(flag, ws) for flag in ["pass", "fail"]]
            self.initDatasets(ws)
            self.initNorm()
            self.initFitPdfs(ws)

            for flag in ["pass", "fail"]: 
                self.doFit(flag)
                if self.results[flag]["status"] is False and self.settings["refit_nobkg"] is True:
                    self.attempt_noBkgFit(flag)
        else:
            pass


    def importFitObjects(self, ws):
        """
        """
        status = bool(self.results["pass"]["status"]*self.results["fail"]["status"])
        if status and (self.existingFit is False):
            ws.Import(self.pdfs["fit_model"]["pass"]), ws.Import(self.pdfs["fit_model"]["fail"])
            ws.Import(self.results["pass"]["res_obj"]), ws.Import(self.results["fail"]["res_obj"])

        

###############################################################################
###############################################################################

class SimFitter(AbsFitter):
    """
    """
    def __init__(self, bin_key, settings):
        """
        """
        AbsFitter.__init__(self, bin_key, settings)
        self.efficiency = ROOT.RooRealVar(f"efficiency_{bin_key}", "Efficiency", 0.9, 0, 1)
        self.__one = ROOT.RooRealVar("one", "one", 1.0)
        self.__one.setConstant()
        self.__minus_one = ROOT.RooRealVar("minus_one", "minus one", -1.0)
        self.__minus_one.setConstant()
        self.one_minus_eff = ROOT.RooPolyVar(
            f"one_minus_eff_{self.bin_key}", "1 - efficiency", 
            self.efficiency, ROOT.RooArgList(self.__one, self.__minus_one))


    def initNorm_sim(self):
        """
        """
        exp_ntot_sig, exp_ntot_up, exp_ntot_low = 0, 0, 0
        for flag, norm_key in [["pass","nsig"], ["fail","nsig"], ["pass","nbkg"], ["fail","nbkg"]]:
            norm_obj = self.norm[norm_key][flag]
            for idx in range(3):
                if (type(norm_obj[idx]) is str) and ("n" in norm_obj[idx]):
                    norm_obj[idx] = float(norm_obj[idx].replace("n",""))*self.histo_data[flag].sumEntries()
                else:
                    norm_obj[idx] = float(norm_obj[idx])
            if norm_key == "nbkg":
                self.norm["nbkg"].update({ flag : ROOT.RooRealVar(
                    f"nbkg_{flag}_{self.bin_key}", f"nbkg {flag}", norm_obj[0], norm_obj[1], norm_obj[2]) })
            else:
                exp_ntot_sig += norm_obj[0]
                exp_ntot_low += norm_obj[1]
                exp_ntot_up += norm_obj[2]
                
        self.norm.update({ "ntot" : ROOT.RooRealVar(
            f"ntot_{self.bin_key}", "nsig total", exp_ntot_sig, exp_ntot_low, exp_ntot_up)})

        self.norm["nsig"].update({
            "pass" : ROOT.RooProduct(f"nsig_pass_{self.bin_key}", "nsig pass", 
                                     ROOT.RooArgList(self.efficiency, self.norm["ntot"])),
            "fail" : ROOT.RooProduct(f"nsig_fail_{self.bin_key}", "nsig fail", 
                                     ROOT.RooArgList(self.one_minus_eff, self.norm["ntot"])) })
        
            
    

    def createSimDataset(self, ws):
        """
        """
        self.sample = ROOT.RooCategory(f"sample_{self.bin_key}", "Composite sample")
        self.sample.defineType("pass")
        self.sample.defineType("fail")

        # The simultaneous dataset will be a one-dimentional histogram, filled with data referring to the
        # different categories. An axis is needed: for the sake of simplicity, is used the one retrieved
        # by the attribute "self.axis["pass"]", that is the same object in memory as "self.axis["fail"]".
        #

        if self.settings["type_analysis"] == "sim":
            self.histo_data.update({ "sim" : ROOT.RooDataHist(
                f"combData_{self.bin_key}", "Combined datasets", 
                ROOT.RooArgSet(self.axis["pass"], self.axis["fail"]), ROOT.RooFit.Index(self.sample), 
                ROOT.RooFit.Import("pass", ws.data(f"Minv_data_pass_{self.bin_key}")),
                ROOT.RooFit.Import("fail", ws.data(f"Minv_data_fail_{self.bin_key}"))) })
        elif self.settings["type_analysis"] == "sim_sf":
            self.sample.defineType("mc_pass")
            self.sample.defineType("mc_fail")
            self.histo_data.update({ "sim" : ROOT.RooDataHist(
                f"combData_{self.bin_key}", "Combined datasets", 
                ROOT.RooArgSet(self.axis["pass"], self.axis["fail"]), ROOT.RooFit.Index(self.sample), 
                ROOT.RooFit.Import("pass", ws.data(f"Minv_data_pass_{self.bin_key}")),
                ROOT.RooFit.Import("fail", ws.data(f"Minv_data_fail_{self.bin_key}")),
                ROOT.RooFit.Import("mc_pass", ws.data(f"Minv_mc_pass_{self.bin_key}")), 
                ROOT.RooFit.Import("mc_fail", ws.data(f"Minv_mc_fail_{self.bin_key}"))) })
        else:
            pass
        
        self.pdfs["fit_model"].update({ "sim" : ROOT.RooSimultaneous(
            f"fitPDF_sim_{self.bin_key}", "Simultaneous pdf", self.sample) })
        
        for flag in ["pass", "fail"]:
            self.pdfs["fit_model"]["sim"].addPdf(self.pdfs["fit_model"][flag], flag)
            if self.settings["type_analysis"] == "sim_sf":
                self.pdfs["fit_model"]["sim"].addPdf(self.pdfs["fit_model"][flag], f"mc_{flag}")   


    def attempt_noBkgFit(self, ws):
        """
        """
        print("\nAttempting to refit with no background\n")

        if self.settings["bkg_model"]["pass"] != "expo" or self.settings["bkg_model"]["fail"] != "expo":
            # Need to be implemented for, e.g., cmsshape
            pass
        else:
            refit_flags = []

            nsig_fitted_pass = self.norm["nsig"]["pass"].getVal()
            nsig_fitted_fail = self.norm["nsig"]["fail"].getVal()

            nbkg_fitted_pass = self.norm["nbkg"]["pass"].getVal()
            nbkg_fitted_fail = self.norm["nbkg"]["fail"].getVal()

            if (nbkg_fitted_pass < 0.005*nsig_fitted_pass): refit_flags.append("pass")  
            if (nbkg_fitted_fail < 0.005*nsig_fitted_fail): refit_flags.append("fail")

            for flag in refit_flags:
                self.pars["tau"][flag].setVal(1), self.pars["tau"][flag].setConstant()
                self.norm["nbkg"][flag].setVal(0), self.norm["nbkg"][flag].setConstant()

            if len(refit_flags) != 0: 
                self.doFit("sim")
            else:
                pass


    def manageFit(self, ws):
        """
        """
        self.checkExistingFit(ws)

        if self.existingFit is False:
            [self.importAxis(flag, ws) for flag in ["pass", "fail"]]
            import_mc = True if self.settings["type_analysis"] == "sim_sf" else False
            self.initDatasets(ws, import_mc=import_mc)
            self.initNorm_sim()
            self.initFitPdfs(ws)
            self.createSimDataset(ws)

            self.doFit("sim")

            if self.results["sim"]["status"] is False and self.settings["refit_nobkg"] is True:
                self.attempt_noBkgFit(ws)
            print(f"\nFinishing fit on bin {self.bin_key}")
        else:
            pass


    def importFitObjects(self, ws):
        """
        """
        status = bool(self.results["sim"]["status"])
        if status and (self.existingFit is False):
            ws.Import(self.pdfs["fit_model"]["sim"]), ws.Import(self.results["sim"]["res_obj"])






            
    