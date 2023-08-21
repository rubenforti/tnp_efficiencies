"""
"""
import ROOT
import os
import sys
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

    def __init__(self, type_analysis, binning_pt, binning_eta, import_ws="", import_txt=""):
        """
        """
        self._dict_results = {}
        self._analysis = type_analysis if type_analysis in ["indep", "sim"] else ""
        if self._analysis == "":
            print("ERROR: analysis type not recognized")
            return None
        
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
        
        elif type(import_txt) is list:
            idx_list = 3
            print(import_txt[0])
            print(import_txt[1])
            print(import_txt[2])
            print(import_txt[3])
            bin_dict = bin_dictionary(binning_pt, binning_eta)
            for bin_key in bin_dict.keys():
                self.add_result_from_txt(import_txt, idx_list, bin_key)
                idx_list += 1


    def add_result(self, res, bin_key):
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
    
    def add_result_from_txt(self, row_list, idx_list, bin_key):
        elements = row_list[idx_list].split('\t')
        eff, deff = float(elements[4]), float(elements[5])
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

def fill_res_histograms(res_bmark, res_new, hist_dict, bin_dict, nbins_eta):
    """
    """

    # Counter used in case pt bins are merged
    cnt_mergedpt = 0

    for bin_key in bin_dict.keys():

        _, bin_pt, bin_eta = bin_dict[bin_key]

        # Bin transformation needed in case the bins are merged
        if type(bin_eta) is list:
                bin_eta = int(1+(nbins_eta*(bin_eta[0]-1)/48.))
        if type(bin_pt) is list:
            bin_pt_list = bin_pt
            bin_pt = int(bin_pt_list[0] - cnt_mergedpt)
            cnt_mergedpt += bin_pt_list[-1]-bin_pt_list[0] if bin_eta==nbins_eta else 0

        eff_1, deff_1 = res_bmark.getEff(bin_key)
        eff_2, deff_2 = res_new.getEff(bin_key)

        if "delta" in hist_dict.keys():
            hist_dict["delta"].Fill(eff_2-eff_1)
            hist_dict["delta_2d"].SetBinContent(bin_pt, bin_eta, eff_2-eff_1)
        if "delta_error" in hist_dict.keys():
            hist_dict["delta_error"].Fill(deff_2-deff_1)
            hist_dict["delta_error_2d"].SetBinContent(bin_pt, bin_eta, deff_2-deff_1)
        if "pull" in hist_dict.keys():
            hist_dict["pull"].Fill((eff_2-eff_1)/deff_2)
            hist_dict["pull_2d"].SetBinContent(bin_pt, bin_eta, (eff_2-eff_1)/deff_2)
        if "rm1" in hist_dict.keys():
            hist_dict["rm1"].Fill((eff_2/eff_1)-1)
            hist_dict["rm1_2d"].SetBinContent(bin_pt, bin_eta, (eff_2/eff_1)-1)
        if "ratio_error" in hist_dict.keys():
            hist_dict["ratio_error"].Fill(deff_2/deff_1)
            hist_dict["ratio_error_2d"].SetBinContent(bin_pt, bin_eta, deff_2/deff_1)
    
        
###############################################################################
###############################################################################

if __name__ == '__main__':

    from fit_utils import fit_quality, pearson_chi2_eval
    from base_library import eval_efficiency, binning, bin_dictionary

    ws_name = "../results/benchmark_iso/ws_iso_indep_benchmark.root"
    
    save_eff_results(ws_name, "indep", "pt", "eta")
    
