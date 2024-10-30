"""
"""
import sys
import ROOT
import pickle
from utilities.binning_utils import bin_dictionary 
from utilities.base_lib import efficiency_from_res, eval_efficiency, sumw2_error, import_pdf_library


class results_manager:
    """
    """

    def __init__(self, input_file, type_analysis, binning_pt, binning_eta, altModel_check=""):
        """
        """
        self._dict_results = {}
        if type_analysis in ["indep", "sim"]:
            self._analysis = type_analysis 
        else:
            sys.exit("ERROR: analysis type not recognized")
        
        if input_file.endswith(".root"):
            import_pdf_library("RooCMSShape")
            f = ROOT.TFile(input_file, "READ")
            input_obj = f.Get("w")
            type_input = ".root"
        elif input_file.endswith(".txt"):
            with open(input_file, "r") as f:
                input_obj = f.readlines()
            type_input = ".txt"
        else:
            sys.exit("ERROR: input file format not recognized")
        
        self.add_results(input_obj, type_input, binning_pt, binning_eta, altModel_check=altModel_check)



    def add_results(self, input_obj, type_input, binning_pt, binning_eta, altModel_check=""):
        """
        """
        idx_list = 3  # Number of first useful row in the txt file (first 3 are comments) 

        for bin_key in bin_dictionary(binning_pt, binning_eta).keys():
            if type_input == ".root":
                self.add_detailed_result(input_obj, bin_key)
            elif type_input == ".txt":
                self.add_detailed_result_from_txt(input_obj, idx_list, bin_key, altModel_check)
                idx_list += 1
        

        
    def add_detailed_result(self, ws, bin_key):
        """
        """
        if self._analysis == 'indep':
            res_pass, res_fail = ws.obj(f"results_pass_{bin_key}"), ws.obj(f"results_fail_{bin_key}")
            self._dict_results.update( {
                bin_key : { "efficiency" : efficiency_from_res(res_pass, res_fail),
                            "efficiency_MC" : eval_efficiency(
                                ws.data(f"Minv_mc_pass_{bin_key}").sumEntries(), ws.data(f"Minv_mc_fail_{bin_key}").sumEntries(),
                                sumw2_error(ws.data(f"Minv_mc_pass_{bin_key}")), sumw2_error(ws.data(f"Minv_mc_fail_{bin_key}"))),
                            "pars_pass" : res_pass.floatParsFinal(),
                            "corrmatrix_pass" : res_pass.correlationMatrix(),
                            "status_pass" : (res_pass.status(), res_pass.covQual(), res_pass.edm()),
                            "pars_fail": res_fail.floatParsFinal(),
                            "corrmatrix_fail": res_fail.correlationMatrix(),
                            "status_fail" : (res_fail.status(), res_fail.covQual(), res_fail.edm()) } })
            
        elif self._analysis == 'sim':
            res = ws.obj(f"results_sim_{bin_key}")
            pars = res.floatParsFinal()
            self._dict_results.update( {
                bin_key: {  "efficiency" : (pars.find(f"efficiency_{bin_key}").getVal(), pars.find(f"efficiency_{bin_key}").getError()),
                            "pars_sim": res.floatParsFinal(),
                            "corrmatrix_sim": res.covarianceMatrix(),
                            "status_sim": (res.status(), res.covQual(), res.edm()) } })
        
        else:
            pass


    def add_detailed_result_from_txt(self, row_list, idx_list, bin_key, altModel_check=""):
        """
        """
        elements = row_list[idx_list].split('\t')
        if altModel_check=="":
            eff, deff = float(elements[4]), float(elements[5])
        elif altModel_check=="altSig":
            eff, deff = float(elements[8]), float(elements[9])
        elif altModel_check=="altBkg":
            eff, deff = float(elements[10]), float(elements[11])
        else:
            sys.exit("ERROR: altModel_check not recognized")

        effMC, deffMC = float(elements[6]), float(elements[7])

        self._dict_results.update({ 
             bin_key : {"efficiency" : (eff, deff), "efficiency_MC" : (effMC, deffMC)}})


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
    

    '''
    ### TO BE APPLIED AND TESTED SOMEWHERE
    def getComparison(self, bin_key='', type="delta", bmark_val=None, bmark_err=None, useTestErrForPull=False):
        """
        """
        if bin_key not in self._dict_results.keys(): sys.exit("Bin key not present in the results object dictionary")
        if bmark_val is None: sys.exit("Benchmark value not provided")

        if useTestErrForPull and type=="pull": bmark_err=self.getEff(bin_key)[1]

        if bmark_err is None and type in ["delta_err", "ratio_err"]: sys.exit("Benchmark error not provided")

        eff, d_eff = self.getEff(bin_key)

        if type == "delta":
            return eff - bmark_val
        elif type == "delta_err":
            return d_eff - bmark_err
        elif type == "pull":
            return (eff - bmark_val) / (d_eff**2 + bmark_err**2)**0.5
        elif type == "pull_ref":
            return (eff - bmark_val) / bmark_err
        elif type == "rm1":
            return (eff/bmark_val) - 1
        elif type == "ratio_err":
            return d_eff / bmark_err
        else:
            sys.exit("Type of comparison not recognized")
    '''
