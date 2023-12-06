"""
"""
import ROOT
import sys
import copy
from utilities.dataset_utils import get_totbkg_roohist
from utilities.results_utils import efficiency_from_res
from utilities.plot_utils import plot_fitted_pass_fail
from utilities.base_library import eval_efficiency, sumw2_error, eval_norm_corrected
from utilities.fit_utils import getSidebands, fit_quality, llr_test_bkg


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
        
        self.histo_data, self.axis, self.constr_set, self.res_obj  = {}, {}, {}, {}
        self.status = {"pass" : True, "fail" : True}
        self.bin_status = False

        self.pdfs = { "smearing" : {}, "mc_pdf" : {}, "conv_pdf" : {}, 
                      "bkg_pdf" : {}, "sum_pdf" : {}, "fit_model" : {} }

        if self.settings["bkg_model"]["pass"] == "mc_raw" or self.settings["bkg_model"]["fail"] == "mc_raw":
            self.bkg_categories = self.settings["bkg_categories"]
            self.bkg_tothist = {}


    def _initParams(self, flag):
        """
        Transforms the list of numbers, representing each parameter (initial
        value, lower bound, upper bound), into a RooRealVar object. 
        """
        for pars_key in self.pars.keys():
            if flag in self.pars[pars_key].keys(): 
                pars_obj = self.pars[pars_key][flag]
                p0, pMin, pMax = pars_obj[0], pars_obj[1], pars_obj[2]
                self.pars[pars_key].update({ 
                    flag : ROOT.RooRealVar(f"{pars_key}_{flag}_{self.bin_key}", f"{pars_key} {flag}", p0, pMin, pMax) })
                if p0==pMin and p0==pMax: 
                    self.pars[pars_key][flag].setConstant()
            else:
                continue


    def _createSigPdf(self, flag, ws):
        """
        Creates the signal pdf and loads it into the pdfs dictionary. The lower
        level pdfs (i.e. the smearing and MC pds) are loaded as well.
        """
        self.pdfs["smearing"].update({ 
            flag : ROOT.RooGaussian(f"smearing_{flag}_{self.bin_key}", "Gaussian smearing", 
                                    self.axis[flag], self.pars["mu"][flag], self.pars["sigma"][flag]) 
            })

        self.pdfs["mc_pdf"].update({
            flag : ROOT.RooHistPdf(f"mc_pdf_{flag}_{self.bin_key}", "MC pdf",
                                   self.axis[flag], ws.data(f"Minv_mc_{flag}_{self.bin_key}"), 3)
            })
        
        if self.settings["type_analysis"] == "sim_sf":
            self.pdfs["mc_pdf"].update({ 
                flag : ROOT.RooExtendPdf(f"conv_{flag}_{self.bin_key}", "pdf MC", 
                                         self.pdfs["mc_pdf"][flag], ws.data(f"Minv_mc_{flag}_{self.bin_key}").sumEntries())
                })
        
        self.pdfs["conv_pdf"].update({
            flag : ROOT.RooFFTConvPdf(f"conv_{flag}_{self.bin_key}", "Convolution pdf", 
                                      self.axis[flag], self.pdfs["mc_pdf"][flag], self.pdfs["smearing"][flag])
            })
        self.pdfs["conv_pdf"][flag].setBufferFraction(self.settings["bufFraction"][flag])

    
    def _createFixedPdfs(self, flag, ws):
        """
        """
        Ndata = self.histo_data[flag].sumEntries()
        Nbkg_raw = ws.data(f"Minv_bkg_{flag}_{self.bin_key}_total").sumEntries()

        #f_bad_charge_assign = 0
        f_ch_assign, f_ch_assign_error = eval_efficiency(ws.data(f"Minv_mc_{flag}_{self.bin_key}_SS").sumEntries(),
                                                         ws.data(f"Minv_mc_{flag}_{self.bin_key}").sumEntries(),
                                                         sumw2_error(ws.data(f"Minv_mc_{flag}_{self.bin_key}_SS")),
                                                         sumw2_error(ws.data(f"Minv_mc_{flag}_{self.bin_key}")))
    
        Nsig_corr, dNsig_corr, Nbkg_corr, dNbkg_corr = eval_norm_corrected(Ndata, Nbkg_raw, f_ch_assign, f_ch_assign_error)

        self.pars["f_ch_assign"] = { 
            flag : ROOT.RooRealVar(f"fraction_ch_{flag}_{self.bin_key}", "fraction charge misassignment",
                                   f_ch_assign, 0.99*f_ch_assign, 1.01*f_ch_assign) }
        self.pars["f_ch_assign"][flag].setError(f_ch_assign_error)

        self.norm["nsig"].update({ 
            flag : ROOT.RooRealVar(f"nsig_{flag}_{self.bin_key}", f"nsig {flag}", Nsig_corr, 0.99*Nsig_corr, 1.01*Nsig_corr) })
        self.norm["nsig"][flag].setError(dNsig_corr), 
        self.norm["nsig"][flag].setConstant()
        
        self.norm["nbkg"].update({ 
            flag : ROOT.RooRealVar(f"nbkg_{flag}_{self.bin_key}", f"nbkg {flag}", Nbkg_corr, 0.99*Nbkg_corr, 1.01*Nbkg_corr) })
        self.norm["nbkg"][flag].setError(dNbkg_corr),
        self.norm["nbkg"][flag].setConstant()

        self.pdfs["conv_pdf"].update({
            flag : ROOT.RooHistPdf(f"conv_{flag}_{self.bin_key}", "Convolution pdf", 
                                   self.axis[flag], ws.data(f"Minv_mc_{flag}_{self.bin_key}"), 0)
            })

        self.pdfs["bkg_pdf"].update({
            flag : ROOT.RooHistPdf(f"bkg_{flag}_{self.bin_key}", "Background pdf", 
                                   self.axis[flag], ws.data(f"Minv_bkg_{flag}_{self.bin_key}_total"), 0)
            })


    def _createBkgPdf(self, flag, ws):
        """
        Creates the background pdf and loads it into the pdfs dictionary. In 
        case the background is modeled with MC, the histogram of the sum of the
        backgrounds is created and loaded in the apposite attribute.
        """
        bkg_model = self.settings["bkg_model"][flag]
        if bkg_model == "expo":
            self.pdfs["bkg_pdf"].update({
                flag : ROOT.RooExponential(f"expo_bkg_{flag}_{self.bin_key}", "Exponential background",
                                           self.axis[flag], self.pars["tau"][flag])
                })            

        elif "cmsshape" in bkg_model:
            self.pdfs["bkg_pdf"].update({
                flag : ROOT.RooCMSShape(f"cmsshape_bkg_{flag}_{self.bin_key}", "CMSShape background",
                                        self.axis[flag], self.pars["alpha"][flag], self.pars["beta"][flag],
                                        self.pars["gamma"][flag], self.pars["peak"][flag])
                })
        
        elif bkg_model == "mc_raw":  #WARNING: the total bkg histo HAS TO BE APPROPRIATELY CREATED BEFORE THE CREATION OF FITTER OBJETCT
            self.pdfs["bkg_pdf"].update({ 
                flag : ROOT.RooHistPdf(f"mcbkg_{flag}_{self.bin_key}_total", "MC-driven background", 
                                       ROOT.RooArgSet(self.axis[flag]),  ws.data(f"Minv_bkg_{flag}_{self.bin_key}_total"))
                })
    
        elif bkg_model == "mc_BBlight":
            # pdf created via barlow beeston strategy
            pass
        
        else:
            sys.exit("ERROR: background model not recognized")
            

    def _createPseudodata(self, flag, ws):
        """
        Creates the pseudodata histogram, which is the sum of the signal and
        background MC histograms. The histogram is loaded into the histo_data
        attribute.
        """
        roohist_pseudodata = ROOT.RooDataHist(f"Minv_pseudodata_{flag}_{self.bin_key}", "pseudodata_histo",
                                              ROOT.RooArgSet(self.axis[flag]), "x_binning")

        roohist_pseudodata.add(ws.data(f"Minv_mc_{flag}_{self.bin_key}"))
        [roohist_pseudodata.add(ws.data(f"Minv_bkg_{flag}_{self.bin_key}_{cat}")) for cat in self.settings["bkg_categories"]]

        self.histo_data.update({flag : roohist_pseudodata})

    def setConstraints(self, flag):
        """
        """
        constr_dict = self.settings["constr"]
        constraints = ROOT.RooArgSet()
        for constr_key in constr_dict.keys():
            if flag in constr_dict[constr_key].keys():
                constr_obj = constr_dict[constr_key][flag]
                constr_dict[constr_key][flag] = {}
                cnt = 0
                for par in ["mu", "sig"]:
                    constr_dict[constr_key][flag].update({ 
                        par : ROOT.RooRealVar(f"{par}Constr_{constr_key}_{flag}_{self.bin_key}", f"{par}Constr {constr_key} {flag}",
                                              constr_obj[cnt], 0.1*constr_obj[cnt], 10*constr_obj[cnt])})
                    constr_dict[constr_key][flag][par].setConstant()
                    cnt+=1

                constr_dict[constr_key][flag].update({
                    "constr_pdf" : ROOT.RooGaussian(f"constr_{constr_key}_{flag}_{self.bin_key}", 
                                                    f"Constraint on {constr_key} {flag}", self.pars[constr_key][flag], 
                                                    constr_dict[constr_key][flag]["mu"], constr_dict[constr_key][flag]["sig"])
                    })
                constraints.add(constr_dict[constr_key][flag]["constr_pdf"])
            else:
                continue
        self.constr_set.update({flag : constraints})


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
                self.res_obj = {
                    "pass": ws.obj(f'results_pass_{self.bin_key}'),
                    "fail": ws.obj(f'results_fail_{self.bin_key}')
                }
                self.existingFit = True
                self.bin_status = True
            else:
                self.existingFit = False
    
        elif self.settings["type_analysis"] == 'sim':
            if type(ws.obj(f'fitPDF_sim_{self.bin_key}')) is ROOT.RooSimultaneous:
                self.res_obj["sim"] = ws.obj(f'results_sim_{self.bin_key}')
                self.existingFit = True
                self.bin_status = True
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
        self.axis[flag].setBinning( ROOT.RooUniformBinning(self.axis[flag].getRange("x_binning")[0], 
                                                           self.axis[flag].getRange("x_binning")[1], 
                                                           self.settings["Nbins"][flag]), "cache" )


    def initDatasets(self, ws):
        """
        """
        for flag in ["pass", "fail"]:

            import_mc = True if self.settings["type_analysis"] == "sim_sf" else False
            import_bkgSS = True if self.settings["bkg_model"][flag] in ["num_estimation", "cmsshape_w_prefitSS"] else False
            import_mcSS = True if self.settings["bkg_model"][flag] == "cmsshape_w_prefitSS" else False

            if self.settings["fit_on_pseudodata"] is False:

                self.histo_data.update({flag : ws.data(f"Minv_data_{flag}_{self.bin_key}")})

                if import_mc is True: 
                    self.histo_data.update({f"mc_{flag}" : ws.data(f"Minv_mc_{flag}_{self.bin_key}")})
                
                if import_mcSS is True:
                    self.histo_data.update({f"mcSS_{flag}" : ws.data(f"Minv_mc_{flag}_{self.bin_key}_SS")})

                if import_bkgSS is True: 
                    self.histo_data.update({f"bkgSS_{flag}" : ws.data(f"Minv_bkg_{flag}_{self.bin_key}_total")})
                    self.histo_data[f"bkgSS_{flag}"].Print("v")
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
            if self.settings["bkg_model"][flag] != "num_estimation":
                self._createSigPdf(flag, ws), self._createBkgPdf(flag, ws)
            else:
                self._createFixedPdfs(flag, ws)
            self.pdfs["sum_pdf"].update({ 
                flag : ROOT.RooAddPdf(f"sum_{flag}_{self.bin_key}", "Signal+Bkg model", 
                                      ROOT.RooArgList(self.pdfs["conv_pdf"][flag], self.pdfs["bkg_pdf"][flag]), 
                                      ROOT.RooArgList(self.norm["nsig"][flag], self.norm["nbkg"][flag]))
                })

            # Creating a new pdf for plotting purposes, I guess there are other ways but this is the simplest
            self.pdfs["fit_model"].update({
                flag : ROOT.RooAddPdf(self.pdfs["sum_pdf"][flag], f'fitPDF_{flag}_{self.bin_key}')
                })


    def doFit(self, flag, fit_bkgSS=False, sidebands_lims=[75,105]):
        """
        """
        par_minos_set = ROOT.RooArgSet()
        if self.settings["type_analysis"] == "sim" and self.settings["useMinos"] == "eff":
            par_minos_set.add(self.efficiency)
        elif self.settings["type_analysis"] == "sim" and self.settings["useMinos"] == "all": 
            par_minos_set.add(self.pdfs["fit_model"][flag].getParameters(self.histo_data[flag]))

        doRegularFit = not (self.settings["bkg_model"][flag]=="num_estimation" or fit_bkgSS)
        
        if doRegularFit:
            pdf_fit = self.pdfs["fit_model"][flag]
            data_fit = self.histo_data[flag]
            range_fit = ""
            subs_res = ""
            checks = self.settings["fit_checks"]

        elif fit_bkgSS:
            pdf_fit = self.pdfs["bkg_pdf"][flag]
            data_fit = self.histo_data[f"bkgSS_{flag}"]
            sidebands_lims = getSidebands(self.histo_data[f"mcSS_{flag}"], self.axis[flag])
            self.axis[flag].setRange("sideband_under", self.axis[flag].getMin(), sidebands_lims[0])
            self.axis[flag].setRange("sideband_over", sidebands_lims[1], self.axis[flag].getMax())
            range_fit = "sideband_under,sideband_over"
            subs_res = "_SS"
            checks = "pseudodata"  # Don't do the chi2 check if it is a prefit on SS data

        elif self.settings["bkg_model"][flag] == "num_estimation":
            pdf_fit = self.pdfs["bkg_pdf"][flag]
            data_fit = self.histo_data[f"bkgSS_{flag}"]
            range_fit = ""
            subs_res = ""
            checks = "pseudodata"   # Don't do the chi2 check if it is a prefit on SS data
        
        else: 
            sys.exit("ERROR: Control applied on fit setting is not consistent!")

        print("")
        print("Data fit", data_fit.GetName())
        print("Pdf fit", pdf_fit.GetName())
        print("Range fit", range_fit)
        print("Subs res", subs_res)
        print("Checks", checks)
        print("")

        res = pdf_fit.fitTo(data_fit,
                            ROOT.RooFit.Range(range_fit),
                            ROOT.RooFit.Minimizer("Minuit2", "Migrad"),
                            # ROOT.RooFit.Minos(par_minos_set),
                            ROOT.RooFit.Strategy(2),
                            ROOT.RooFit.SumW2Error(False),
                            ROOT.RooFit.ExternalConstraints(self.constr_set[flag]),
                            ROOT.RooFit.Save(1), 
                            ROOT.RooFit.PrintLevel(self.settings["fit_verb"]))
        
        res.SetName(f"results_{flag}_{self.bin_key}{subs_res}")

        self.res_obj.update({ flag : res }) 

        if self.settings["bkg_model"][flag] == "num_estimation": self.fixFitObj(flag)

        fit_obj = { "axis" : self.axis[flag], "histo" : data_fit, 
                    "pdf" : pdf_fit, "res" : res }
        
        self.status.update({ flag : fit_quality(fit_obj, type_checks=checks) }) 
        print(res.status(), res.covQual(), res.edm())
        print(f"Fit status for {flag} category: {self.status[flag]}")

        if fit_bkgSS:
            if self.status[flag]:
                for par in self.pdfs["bkg_pdf"][flag].getParameters(data_fit): par.setConstant()
            else:
                print(f"WARNING: fit on SS data failed, the reported efficiency value is calculated using the number of events in the {flag} category")
                self.norm["nsig"][flag].setVal(self.histo_data[flag].sumEntries())
                self.norm["nsig"][flag].setError(self.histo_data[flag].sumEntries()**0.5)
                self.res_obj[flag].setFinalParList(ROOT.RooArgList(self.norm["nsig"][flag]))

        
    def fixFitObj(self, flag):
        """
        """
        # self.res_obj[flag] = ROOT.RooFitResult(f"results_{flag}_{self.bin_key}", f"results_{flag}_{self.bin_key}")
        self.res_obj[flag].setStatus(0)
        self.res_obj[flag].setCovQual(3)
        self.res_obj[flag].setEDM(0)
        self.res_obj[flag].setFinalParList(
            ROOT.RooArgList(self.pars["f_ch_assign"][flag], self.norm["nbkg"][flag], self.norm["nsig"][flag]))


    def saveFig(self, ws, figpath):
        """
        """
        if self.settings["type_analysis"] == "indep":
            eff, d_eff = efficiency_from_res(self.res_obj["pass"], self.res_obj["fail"])
            fig_status = self.bin_status
        if self.settings["type_analysis"] == "sim":
            eff, d_eff = self.efficiency.getVal(), self.efficiency.getError()
            fig_status = bool(self.status["sim"])
            self.res_obj.update({
                "pass" : self.res_obj["sim"], "fail" : self.res_obj["sim"]["res_obj"]})
        else:
            pass
        # ATTENZIONE: CONTROLLARE CHE SUCCEDE SE SI FA IL FIT CON L'OPZIONE "SIM_SF"

        eff_mc, deff_mc = eval_efficiency(ws.data(f"Minv_mc_pass_{self.bin_key}").sumEntries(), 
                                          ws.data(f"Minv_mc_fail_{self.bin_key}").sumEntries(),
                                          sumw2_error(ws.data(f"Minv_mc_pass_{self.bin_key}")),
                                          sumw2_error(ws.data(f"Minv_mc_fail_{self.bin_key}")))
        
        scale_factor = eff/eff_mc
        d_scale_factor = scale_factor*((d_eff/eff)**2 + (deff_mc/eff_mc)**2)**0.5

        plot_objects={}
        fig_folder = figpath["good"] if fig_status is True else figpath["check"]
        for flag in ["pass", "fail"]:
            plot_objects.update({ 
                flag : {"axis":self.axis[flag], "data":self.histo_data[flag],
                        "model":self.pdfs["fit_model"][flag], "res":self.res_obj[flag]}
                })

        plot_objects.update({"efficiency" : [eff, d_eff],
                             "efficiency_mc" : [eff_mc, deff_mc],
                             "scale_factor" : [scale_factor, d_scale_factor]})
        
        plot_fitted_pass_fail(self.settings["type_analysis"], plot_objects, self.bin_key, figpath=fig_folder)


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
        print("\nAttempting to refit with no background\n")
        # Implementare una funzione che gestisce queste due possibilità
        if self.settings["type_refit"] == "benchmark":
            nsig_fitted, nbkg_fitted = self.norm["nsig"][flag], self.norm["nbkg"][flag]
            refit_ctrl = nbkg_fitted.getVal() < 0.005*nsig_fitted.getVal()
        elif self.settings["type_refit"] == "llr":
            # NOT TESTED YET
            # refit_ctrl = llr_test_bkg(self.histo_data[flag], self.pdfs["fit_model"][flag], alpha=0.05)
            pass
        else:
            sys.exit("ERROR: refit type not recognized")

        if refit_ctrl is True:
            self.norm["nbkg"][flag].setVal(0)
            self.norm["nbkg"][flag].setConstant()
            for par in self.pdfs["bkg_pdf"][flag].getParameters(self.histo_data[flag]): par.setConstant()
            self.doFit(flag)
        else: 
            print("Refit not triggered\n")


    def manageFit(self, ws):
        """
        """
        self.checkExistingFit(ws)  

        if self.existingFit is False:

            [self.importAxis(flag, ws) for flag in ["pass", "fail"]], print("Imported axis")
            self.initDatasets(ws), print("Imported datasets")
            self.initNorm(),       print("Imported normalization parameters")
            self.initFitPdfs(ws),  print("Imported fit pdfs")
            
            for flag in ["pass", "fail"]: 

                if "constr" in self.settings.keys(): self.setConstraints(flag)

                if self.settings["bkg_model"][flag] == "cmsshape_w_prefitSS":
                    print("Fixing bkg parameters by operating a fit on SS data with CMSShape")
                    self.doFit(flag, fit_bkgSS=True) 
                    for par in self.pdfs["bkg_pdf"][flag].getParameters(self.histo_data[flag]):
                        par.Print()
        
                if self.status[flag] is True: 
                    print(f"Starting fit for {flag} category")
                    self.doFit(flag)
                    if self.status[flag] is False and self.settings["refit_nobkg"]:
                        self.attempt_noBkgFit(flag)

            self.bin_status = bool(self.status["pass"]*self.status["fail"])
            print(f"Fitted bin {self.bin_key} with status {self.bin_status}\n")
        else:
            print("Not possible to refit an existing PDF! \nUsing the results obtained previously")


    def importFitObjects(self, ws):
        """
        """
        if self.bin_status and (self.existingFit is False):
            ws.Import(self.pdfs["fit_model"]["pass"]), ws.Import(self.pdfs["fit_model"]["fail"])
            ws.Import(self.res_obj["pass"]), ws.Import(self.res_obj["fail"])


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
        self.one_minus_eff = ROOT.RooPolyVar(f"one_minus_eff_{self.bin_key}", "1 - efficiency", 
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
        self.sample.defineType("pass"), self.sample.defineType("fail")

        # The simultaneous dataset will be a one-dimentional histogram, filled with data referring to the
        # different categories. An axis is needed: for the sake of simplicity, is used the one retrieved
        # by the attribute "self.axis["pass"]", that is the same object in memory as "self.axis["fail"]".
        if self.settings["type_analysis"] == "sim":
            self.histo_data.update({ 
                "sim" : ROOT.RooDataHist(f"combData_{self.bin_key}", "Combined datasets", 
                                         ROOT.RooArgSet(self.axis["pass"], self.axis["fail"]),
                                         ROOT.RooFit.Index(self.sample), 
                                         ROOT.RooFit.Import("pass", ws.data(f"Minv_data_pass_{self.bin_key}")),
                                         ROOT.RooFit.Import("fail", ws.data(f"Minv_data_fail_{self.bin_key}"))) 
                })
        elif self.settings["type_analysis"] == "sim_sf":
            self.sample.defineType("mc_pass"), self.sample.defineType("mc_fail")
            self.histo_data.update({ 
                "sim" : ROOT.RooDataHist(f"combData_{self.bin_key}", "Combined datasets", 
                                         ROOT.RooArgSet(self.axis["pass"], self.axis["fail"]), 
                                         ROOT.RooFit.Index(self.sample), 
                                         ROOT.RooFit.Import("pass", ws.data(f"Minv_data_pass_{self.bin_key}")),
                                         ROOT.RooFit.Import("fail", ws.data(f"Minv_data_fail_{self.bin_key}")),
                                         ROOT.RooFit.Import("mc_pass", ws.data(f"Minv_mc_pass_{self.bin_key}")), 
                                         ROOT.RooFit.Import("mc_fail", ws.data(f"Minv_mc_fail_{self.bin_key}")))
                })
        else:
            sys.exit("ERROR: type_analysis not recognized, retry with 'sim' or 'sim_sf'")

        
        self.pdfs["fit_model"].update({ 
            "sim" : ROOT.RooSimultaneous(f"fitPDF_sim_{self.bin_key}", "Simultaneous pdf", self.sample)
            })
        
        typedata_flag = ["pass", "fail", "mc_pass", "mc_fail"] if "_sf" in self.settings["type_analysis"] else ["pass", "fail"]

        for flag in typedata_flag: self.pdfs["fit_model"]["sim"].addPdf(self.pdfs["fit_model"][flag], flag)
            

    def attempt_noBkgFit(self, ws):
        """
        """
        print("\nAttempting to refit with no background\n")

        if self.settings["bkg_model"]["pass"] == "expo" and self.settings["bkg_model"]["fail"] == "expo":
            
            refit_flags = []

            nsig_fitted_pass, nbkg_fitted_pass = self.norm["nsig"]["pass"], self.norm["nbkg"]["pass"]
            nsig_fitted_fail, nbkg_fitted_fail = self.norm["nsig"]["fail"].getVal(), self.norm["nbkg"]["fail"] 

            # Rifare come sopra per il caso di fit indipendenti
            if (nbkg_fitted_pass.getVal() < 0.005*nsig_fitted_pass.getVal()): refit_flags.append("pass")  
            if (nbkg_fitted_fail.getVal() < 0.005*nsig_fitted_fail.getVal()): refit_flags.append("fail")

            for flag in refit_flags:  #Sistemare in maniera più elegante
                self.pars["tau"][flag].setVal(1), self.pars["tau"][flag].setConstant()
                self.norm["nbkg"][flag].setVal(0), self.norm["nbkg"][flag].setConstant()

            if len(refit_flags) != 0: 
                self.doFit("sim")
        else:
            sys.exit("ERROR: refit with no background is not possible for this background model")


    def manageFit(self, ws):
        """
        """
        self.checkExistingFit(ws)

        if self.existingFit is False:
            [self.importAxis(flag, ws) for flag in ["pass", "fail"]]
            self.initDatasets(ws)
            self.initNorm_sim()
            self.initFitPdfs(ws)
            self.createSimDataset(ws)

            self.doFit("sim")

            if self.status["sim"] is False and self.settings["refit_nobkg"]: 
                self.attempt_noBkgFit(ws)

            print(f"Fitted bin {self.bin_key} with status {self.status['sim']}\n")
        else:
            print("Not possible to refit an existing PDF! \nUsing the results obtained previously")


    def importFitObjects(self, ws):
        """
        """
        if self.status["sim"] and (self.existingFit is False):
            ws.Import(self.pdfs["fit_model"]["sim"]), ws.Import(self.res_obj["sim"])






            
    