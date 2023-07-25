"""
"""
import ROOT
import os
import pickle
from array import array
from utilities.fit_utils import fit_quality, pearson_chi2_eval
from utilities.base_library import eval_efficiency, binning, bin_dictionary


def efficiency_from_res(res_pass, res_fail):
    """
    """
    pars_pass, pars_fail = res_pass.floatParsFinal(), res_fail.floatParsFinal()

    for par in pars_pass:
        if "nsig" in par.GetName():
            npass = par
    for par in pars_fail:
        if "nsig" in par.GetName():
            nfail = par

    eff, d_eff = eval_efficiency(npass.getVal(), nfail.getVal(), npass.getError(), nfail.getError())

    return eff, d_eff

###############################################################################

class results_manager:
    """
    """

    def __init__(self, type_analysis, binning_pt, binning_eta, import_ws = ""):
        """
        """
        self._dict_results = {}
        self._analysis = type_analysis if type_analysis in ["indep", "sim"] else ""
        
        if type(import_ws) is ROOT.RooWorkspace:
            bin_dict = bin_dictionary(binning_pt, binning_eta)
            for bin_key in bin_dict.keys():

                if self._analysis == 'indep':
                    res_pass = import_ws.obj(f"results_pass_{bin_key}")
                    res_fail = import_ws.obj(f"results_fail_{bin_key}")
                    self.add_result({"pass":res_pass, "fail":res_fail}, bin_key)
                elif self._analysis == 'sim':
                    self.add_result({"sim":import_ws.obj(f"results_{bin_key}")}, bin_key)
                else:
                    print("ERROR: analysis type not recognized")
                

    def add_result(self, res, bin_key):
        """
        """
        Npass, sigma_Npass = 0, 0
        Nfail, sigma_Nfail = 0, 0

        if self._analysis == 'indep':
            # res_pass = ws.obj(f"results_pass_({bin_pt}|{bin_eta})")
            # res_fail = ws.obj(f"results_fail_({bin_pt}|{bin_eta})")

            res_pass, res_fail = res["pass"], res["fail"]
            
            # histo_pass = self._ws.data(f"Minv_data_pass_({bin_pt}|{bin_eta})")
            # histo_fail = self._ws.data(f"Minv_data_fail_({bin_pt}|{bin_eta})")

            if type(res_pass) is ROOT.RooFitResult and type(res_fail) is ROOT.RooFitResult:

                new_res = {
                    bin_key : {
                        "efficiency" : efficiency_from_res(res_pass, res_fail),
                        "pars_pass" : res_pass.floatParsFinal(),
                        "corrmatrix_pass" : res_pass.correlationMatrix(),
                        "status_pass" : (res_pass.status(), res_pass.covQual(), res_pass.edm()),
                        "pars_fail": res_fail.floatParsFinal(),
                        "corrmatrix_fail": res_fail.correlationMatrix(),
                        "status_fail" : (res_fail.status(), res_fail.covQual(), res_fail.edm())
                        }
                    }
                self._dict_results.update(new_res)

        elif self._analysis == 'sim':
            # res = ws.obj(f"results_({bin_pt}|{bin_eta})")
            results = res["sim"]
            #new_res = {f"{bin_pt},{bin_eta}": {
            new_res = {f"{bin_key}": {
                # "efficiency" : (res.floatParsFinal().find(f"efficiency_({bin_pt}|{bin_eta})").getVal(),
                #                 res.floatParsFinal().find(f"efficiency_({bin_pt}|{bin_eta})").getError()),
                "efficiency" : (results.floatParsFinal().find(f"efficiency_{bin_key}").getVal(),
                                results.floatParsFinal().find(f"efficiency_{bin_key}").getError()),
                "parameters": results.floatParsFinal(),
                "corrmatrix": results.covarianceMatrix(),
                "status": (results.status(), results.covQual(), results.edm())
                }
            }
            self._dict_results.update(new_res)
        else:
            pass

    def Open(self, filename):
        with open(filename, "rb") as file:
            self._dict_results = pickle.load(file)

    def Write(self, filename):
        with open(filename, "wb") as file:
            pickle.dump(self._dict_results, file)
            # file.close()
    
    def dictionary(self):
        return self._dict_results

    def getEff(self, bin_key='', bin_pt=0, bin_eta=0):
        """
        """
        if bin_key == "all":
            eff = [] 
            [eff.append(self._dict_results[key]["efficiency"]) for key in self._dict_results.keys()]
        elif bin_key=='' and bin_pt!=0 and bin_eta != 0:
            pass
        else:
            eff = self._dict_results[bin_key]["efficiency"]

        return eff
    
    def getStatus(self, bin_key=''):
        """
        """
        b_keys = [bin_key] if bin_key!='' else self._dict_results.keys()
        status = {"pass":[], "fail":[]} if self._analysis == 'indep' else []
        for key in b_keys:
            if self._analysis == 'indep':
                status["pass"].append(self._dict_results[key]["status_pass"][0])
                status["fail"].append(self._dict_results[key]["status_fail"][0])
            elif self._analysis == 'sim':
                status.append(self._dict_results[key]["status"][0])

        return status
    
    def getCovQual(self, bin_key=''):
        """
        """
        b_keys = [bin_key] if bin_key!='' else self._dict_results.keys()
        covq = {"pass":[], "fail":[]} if self._analysis == 'indep' else []
        for key in b_keys:
            if self._analysis == 'indep':
                covq["pass"].append(self._dict_results[key]["status_pass"][1])
                covq["fail"].append(self._dict_results[key]["status_fail"][1])
            elif self._analysis == 'sim':
                covq.append(self._dict_results[key]["status"][1])

        return covq
    
    def getEDM(self, bin_key=''):
        """
        """
        b_keys = [bin_key] if bin_key!='' else self._dict_results.keys()
        edm = {"pass":[], "fail":[]} if self._analysis == 'indep' else []
        for key in b_keys:
            if self._analysis == 'indep':
                edm["pass"].append(self._dict_results[key]["status_pass"][2])
                edm["fail"].append(self._dict_results[key]["status_fail"][2])
            elif self._analysis == 'sim':
                edm.append(self._dict_results[key]["status"][2])

        return edm
    
    def getPValue(self, bin_key=""):
        """
        TO BE IMPLEMENTED
        """
        b_keys = [bin_key] if bin_key!='' else self._dict_results.keys()
        pvalue = {"pass":[], "fail":[]} if self._analysis == 'indep' else []
        for key in b_keys:
            if self._analysis == 'indep':
                pass
            elif self._analysis == 'sim':
                pass
        return pvalue

###############################################################################

def init_pass_fail_histos(histo_name, histo_title, bins_var, bins_pt, bins_eta):

    histos = {}

    for flag in ["pass", "fail"]:
        h1d = ROOT.TH1D(f"{histo_name}_{flag}", f"{histo_title} - {flag}ing probes", 
                        len(bins_var)-1, bins_var)
        h2d = ROOT.TH2D(f"{histo_name}_{flag}_2d", f"{histo_title} - {flag}ing probes", 
                        len(bins_pt)-1, bins_pt, len(bins_eta)-1, bins_eta)
        histos.update({f"{histo_name}_{flag}" : h1d})
        histos.update({f"{histo_name}_{flag}_2d" : h2d})

    return histos
    
###############################################################################

def init_results_histos(histo_name, histo_title, bins_var, bins_pt, bins_eta):

    histos = {}

    h1d = ROOT.TH1D(f"h_{histo_name}", histo_title, len(bins_var)-1, bins_var)
    h2d = ROOT.TH2D(f"h_{histo_name}_2d", f"{histo_title} ",
                    len(bins_pt)-1, bins_pt, len(bins_eta)-1, bins_eta)
    histos.update({f"{histo_name}" : h1d})
    histos.update({f"{histo_name}_2d" : h2d})

    return histos

###############################################################################

def save_eff_results(ws_name, type_analysis, binning_pt, binning_eta):
    """
    """

    file_in = ROOT.TFile.Open(ws_name, "UPDATE")
    ws = file_in.Get("w")

    bins_pt, bins_eta = binning(binning_pt), binning(binning_eta)

    bin_dict = bin_dictionary(binning_pt, binning_eta)
    
    
    results = results_manager(type_analysis, binning_pt, binning_eta, import_ws=ws)

    bins_eff = array("d", [round(0.85 + 0.003*i, 3) for i in range(51)])
    bins_rel_err = array("d", [round(0 + 0.0002*i, 4) for i in range(51)])

    histos ={}
    histos.update(init_results_histos("efficiency", "Efficiency", 
                                      bins_eff, bins_pt, bins_eta))
    histos.update(init_results_histos("eff_rel_error", "Relative error on efficiency",
                                      bins_rel_err, bins_pt, bins_eta))
    
    for bin_key in bin_dict.keys():
        
        _, bin_pt, bin_eta = bin_dict[bin_key]

        eff, d_eff = results.getEff(bin_key)
        # print(eff[0], eff[1])
        histos["efficiency"].Fill(eff)
        histos["efficiency_2d"].SetBinContent(bin_pt, bin_eta, eff)
        histos["eff_rel_error"].Fill(d_eff/eff)
        histos["eff_rel_error_2d"].SetBinContent(bin_pt, bin_eta, d_eff/eff)
    

    file_out = ROOT.TFile.Open(ws_name.replace("ws", "hres"), "RECREATE")

    [histo.Write() for histo in histos.values()]

    file_out.Close()

    
        
###############################################################################
###############################################################################

if __name__ == '__main__':


    ws_name = "../root_files/ws_iso_indep_bmark_2gev.root"
    
    save_eff_results(ws_name, "indep", "pt", "eta")
    
