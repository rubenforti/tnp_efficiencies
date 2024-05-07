"""
"""
import sys
import ROOT
import pickle
from utilities.binning_utils import bin_dictionary 
from utilities.base_lib import efficiency_from_res, eval_efficiency, sumw2_error


class results_manager:
    """
    """

    def __init__(self, type_analysis, binning_pt, binning_eta, import_ws="", import_txt="", altSig_check=False):
        """
        """
        self._dict_results = {}
        self._analysis = type_analysis if type_analysis in ["indep", "sim"] else ""
        if self._analysis == "":
            print("ERROR: analysis type not recognized")
            return None
        
        bin_dict = bin_dictionary(binning_pt, binning_eta)
        
        idx_list = 3  # Number of first useful row in the txt file (first 3 are comments) 

        for bin_key in bin_dict.keys():
        
            if type(import_ws) is ROOT.RooWorkspace:
                if self._analysis == 'indep':
                    res_pass = import_ws.obj(f"results_pass_{bin_key}")
                    res_fail = import_ws.obj(f"results_fail_{bin_key}")
                    self.add_result({"pass":res_pass, "fail":res_fail}, bin_key, import_ws)
                elif self._analysis == 'sim':
                    self.add_result({"sim":import_ws.obj(f"results_sim_{bin_key}")}, bin_key, import_ws)
                else:
                    print("ERROR: analysis type not recognized")
        
            elif type(import_txt) is list:
                if altSig_check is False:
                    self.add_result_from_txt(import_txt, idx_list, bin_key)
                else:
                    self.add_result_from_txt_altSig(import_txt, idx_list, bin_key)
                idx_list += 1


    def add_result(self, res, bin_key, ws):
        """
        """
        Npass, sigma_Npass = 0, 0
        Nfail, sigma_Nfail = 0, 0

        if self._analysis == 'indep':
        
            res_pass, res_fail = res["pass"], res["fail"]
            if type(res_pass) is ROOT.RooFitResult and type(res_fail) is ROOT.RooFitResult:

                new_res = {
                    bin_key : {
                        "efficiency" : efficiency_from_res(res_pass, res_fail),
                        "efficiency_MC" : eval_efficiency(ws.data(f"Minv_mc_pass_{bin_key}").sumEntries(),
                                                          ws.data(f"Minv_mc_fail_{bin_key}").sumEntries(),
                                                          sumw2_error(ws.data(f"Minv_mc_pass_{bin_key}")),
                                                          sumw2_error(ws.data(f"Minv_mc_fail_{bin_key}"))),
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
            results = res["sim"]
            if type(results) is not ROOT.RooFitResult:
                print("ERROR: result object not recognized")
                print(bin_key)
                sys.exit()
            pars = results.floatParsFinal()
            new_res = {f"{bin_key}": {
                "efficiency" : (pars.find(f"efficiency_{bin_key}").getVal(),
                                pars.find(f"efficiency_{bin_key}").getError()),
                "pars_sim": results.floatParsFinal(),
                "corrmatrix_sim": results.covarianceMatrix(),
                "status_sim": (results.status(), results.covQual(), results.edm())
                }
            }
            self._dict_results.update(new_res)
        else:
            pass
    
    def add_result_from_txt(self, row_list, idx_list, bin_key):
        elements = row_list[idx_list].split('\t')
        eff, deff = float(elements[4]), float(elements[5])
        effMC, deffMC = float(elements[6]), float(elements[7])
        self._dict_results.update({ 
            bin_key : {"efficiency" : (eff, deff), "efficiency_MC" : (effMC, deffMC)}})


    def add_result_from_txt_altSig(self, row_list, idx_list, bin_key):
        elements = row_list[idx_list].split('\t')
        eff, deff = float(elements[8]), float(elements[9])
        self._dict_results.update({bin_key : {"efficiency" : (eff, deff)}})

    def Open(self, filename):
        with open(filename, "rb") as file:
            self._dict_results = pickle.load(file)

    def Write(self, filename):
        with open(filename, "wb") as file:
            pickle.dump(self._dict_results, file)
            # file.close()
    
    def dictionary(self):
        return self._dict_results

    def getEff(self, bin_key=''):
        """
        """
        if bin_key == "all":
            eff = [] 
            [eff.append(self._dict_results[key]["efficiency"]) for key in self._dict_results.keys()]
        elif bin_key in self._dict_results.keys():
            eff = self._dict_results[bin_key]["efficiency"]
        else:
            print("Bin key not present in the results object dictionary")
            sys.exit()

        return eff
    
    def getEffMC(self, bin_key=''):
        """
        """
        if bin_key == "all":
            eff = [] 
            [eff.append(self._dict_results[key]["efficiency_MC"]) for key in self._dict_results.keys()]
        elif bin_key in self._dict_results.keys():
            eff = self._dict_results[bin_key]["efficiency_MC"]
        else:
            print("Bin key not present in the results object dictionary")
            sys.exit()

        return eff
            
    def getPars(self, flag, bin_key=''):
        """
        """
        if bin_key == "all":
            pars = [] 
            [pars.append(self._dict_results[key][f"pars_{flag}"]) for key in self._dict_results.keys()]
        elif bin_key in self._dict_results.keys():
            pars = self._dict_results[bin_key][f"pars_{flag}"]
        else:
            print("Bin key not present in the results object dictionary")
            sys.exit()

        return pars
    
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
