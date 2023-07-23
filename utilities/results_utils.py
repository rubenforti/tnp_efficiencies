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
    
    
    def printStatus(self, bin_key=''):
        """
        """
        if bin_key=='':
            bin_keys = self._dict_results.keys()
        elif type(bin_key) is str and bin_key!='':
            bin_keys = [bin_key]
        else:
            print("ERROR: bin_key not recognized")
            sys.exit()

        if self._analysis == 'indep':
            [print(f"{b_key}: {self._dict_results[b_key]['status_pass']}, {self._dict_results[b_key]['status_fail']}")
                for b_key in bin_keys]
        elif self._analysis == 'sim':
            [print(f"{b_key}: {self._dict_results[b_key]['status']}") for b_key in bin_keys]
        else:
            pass





    '''
    # THE FOLLOWING METHODS ARE NOT NECCESSARY ANYMORE, DUE TO THE DIFFERENT ROLE ASSIGNED TO THE CLASS
    # IT'S BETTER TO KEEP THEM HERE AS REFERENCES FOR THE FUTURE

    def getCorrmatrix(self, bin_pt=0, bin_eta=0):
        """
        """
        if bin_pt == 0 and bin_eta == 0:
            corr = [] 
            if self._analysis == 'indep':
                [corr.append(self._dict_results[key]["corrmatrix_fail"], 
                             self._dict_results[key]["corrmatrix_pass"]) for key in self._dict_results.keys()]
            elif self._analysis == 'sim':
                [corr.append(self._dict_results[key]["corrmatrix"]) for key in self._dict_results.keys()]
        else:
            if self._analysis == 'indep':
                corr = (self._dict_results[f"{bin_pt},{bin_eta}"]["corrmatrix_fail"], 
                        self._dict_results[f"{bin_pt},{bin_eta}"]["corrmatrix_pass"])
            elif self._analysis == 'sim':
                corr = self._dict_results[f"{bin_pt},{bin_eta}"]["corrmatix"]

        return corr

    def getStatus(self, bin_pt, bin_eta, conditions='all'):
        """
        HAS TO BE EXTENDED FOR SIMULTANEOUS FITS
        """
        res = self._dict_results
        key = f"{bin_pt},{bin_eta}"

        eff = res[key]["efficiency"]
        migr_p = res[key]["fit_stat_pass"]["migrad_status"]
        covq_p = res[key]["fit_stat_pass"]["cov_matrix_quality"]
        edm_p = res[key]["fit_stat_pass"]["EDM"]
        migr_f = res[key]["fit_stat_fail"]["migrad_status"]
        covq_f = res[key]["fit_stat_fail"]["cov_matrix_quality"]
        edm_f = res[key]["fit_stat_fail"]["EDM"]

        if conditions == 'all':
            conditions = ['eff', 'migrad', 'covqual', 'edm']

        cond_list = []

        if 'eff' in conditions:
            cond_list.append(eff[0] > 0 and eff[0] < 1)
        if 'migrad' in conditions:
            cond_list.append(migr_p == 0 and migr_f == 0)
        if 'covqual' in conditions:
            cond_list.append(covq_p == 3 and covq_f == 3)
        if 'edm' in conditions:
            cond_list.append(edm_p < 1e-4 and edm_f < 1e-4)

        nfails = 0
        for cond in cond_list:
            if cond is False:
                nfails = nfails + 1

        return nfails

    def get_problematic_bins(self):
        """
        HAS TO BE EXTENDED FOR SIMULTANEOUS FITS
        """
        res = self._dict_results
        bins = []
        for key in res.keys():
            bin_pt, bin_eta = key.split(',')
            status = self.check_fit_status(int(bin_pt), int(bin_eta))
            if status != 0:
                print(f'Bin {key} has {status} problems')
                bins.append(key)
        return bins

    '''
    
###############################################################################

def init_pass_fail_histos(histo_name, histo_title, binning):

    histos = {}

    for flag in ["pass", "fail"]:

        h1d = ROOT.TH1D(f"{histo_name}_{flag}", f"{histo_title} - {flag}ing probes", len(binning)-1, binning)
        histos.update({f"{histo_name}_{flag}" : h1d})

    return histos

###############################################################################

def init_pass_fail_h2d(histo_name, histo_title, binning_pt, binning_eta):

    histos = {}

    for flag in ["pass", "fail"]:

        h2d = ROOT.TH2D(f"{histo_name}_{flag}", f"{histo_title} - {flag}ing probes", 
                        len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta)
        histos.update({f"{histo_name}_{flag}" : h2d})

    return histos
    
###############################################################################

def compare_with_benchmark(results, ref_txt):

    with open(ref_txt, "r") as file:
        row_list = file.readlines()
    

    file = ROOT.TFile("prova_histo2d.root", "RECREATE")

    
    binning_pt = array('d', [24., 26., 28., 30., 32., 34.,
                       36., 38., 40., 42., 44., 47., 50., 55., 60., 65.])
    binning_eta = array('d', [round(-2.4 + i*0.1, 2) for i in range(49)])
    
    h_delta_eff = ROOT.TH1D("delta_eff", "delta eff", 50, -1e-4, 1e-4)
    h2d_delta_eff = ROOT.TH2D("delta_eff_2d", "delta eff", len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta)
    h_delta_deff = ROOT.TH1D("delta_deff", "delta deff", 50, -1e-5, 1e-5)
    h2d_error_delta_eff = ROOT.TH2D("error_delta_eff_2d", "error delta eff", len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta)
    h2d_pull = ROOT.TH2D("pull_delta_eff_2d", "pull delta eff", len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta)

    idx_list = 3

    print(row_list[0])
    print(row_list[1])
    print(row_list[2])
    print(row_list[3])

    for i in range(1, 16):
        for j in range(1, 49):
            eff, deff = results.getEfficiency(i,j)
            elements = row_list[idx_list].split('\t')

            h_delta_eff.Fill(eff-float(elements[4]))
            h_delta_deff.Fill(deff-float(elements[5]))

            h2d_delta_eff.SetBinContent(i, j, eff-float(elements[4]))
            h2d_error_delta_eff.SetBinContent(i, j, float(elements[5]))

            h2d_pull.SetBinContent(i, j, abs(eff-float(elements[4]))/float(elements[5]))

            idx_list += 1

    c0 = ROOT.TCanvas("", "", 1200, 900)
    c0.Divide(2)
    c0.cd(1)
    ROOT.gStyle.SetOptStat("menr")
    ROOT.gPad.SetLogy()
    h_delta_eff.Draw()
    c0.cd(2)
    ROOT.gPad.SetLogy()
    h_delta_deff.Draw()
    c0.SaveAs("figs/delta_eff.pdf")

    c1 = ROOT.TCanvas("delta_eff", "delta_eff", 1200, 900)
    c1.cd()
    ROOT.gStyle.SetOptStat("en")
    ROOT.gPad.SetRightMargin(0.15)
    ROOT.gStyle.SetPalette(57)
    h2d_delta_eff.Draw("COLZ")
    h2d_delta_eff.SetContour(25)
    c1.SaveAs("figs/delta_eff_2d.pdf")

    c2 = ROOT.TCanvas("pull_delta_eff", "pull_delta_eff", 1200, 900)
    c2.cd()
    ROOT.gStyle.SetOptStat("en")
    ROOT.gPad.SetRightMargin(0.15)
    ROOT.gStyle.SetPalette(57)
    h2d_pull.Draw("COLZ")
    h2d_pull.SetContour(25)
    c2.SaveAs("figs/pull_delta_eff_2d.pdf")


    h2d_delta_eff.Write()
    h2d_error_delta_eff.Write()
    h2d_pull.Write()

    file.Close()

###############################################################################

def compare_eff_results(ws_bmark, ws_new, binning_pt, binning_eta, file_output, aux_filename="", aux_dict={}):
    """
    Compare the efficiencies and their error between two results files. The first file is 
    considered as benchmark
    """
    
    for t_an in ["indep", "sim"]:
        an_bmark = t_an if t_an in ws_bmark.GetName() else ""
        an_new = t_an if t_an in ws_new.GetName() else ""

    
    res_benchmark = results_manager("indep", binning_pt, binning_eta, import_ws=ws_bmark)
    res_new = results_manager("indep", binning_pt, binning_eta, import_ws=ws_new)

    if aux_filename!="":
        aux_file = ROOT.TFile(aux_filename, "READ")
        aux_ws = aux_file.Get("w")


    for key_aux in aux_dict.keys():

        res_pass = aux_ws.obj(f"results_pass_{aux_dict[key_aux]}")
        res_fail = aux_ws.obj(f"results_fail_{aux_dict[key_aux]}")
        new_res = {"pass":res_pass, "fail":res_fail}

        res_new.add_result(new_res, key_aux)

    file_out = ROOT.TFile(file_output, "RECREATE")

    bins_pt, bins_eta = binning(binning_pt), binning(binning_eta)
    nbins_pt, nbins_eta = len(bins_pt)-1, len(bins_eta)-1

    bin_dict = bin_dictionary(binning_pt, binning_eta)
    
    bmark_dict = res_benchmark.dictionary()
 


    h_delta_eff = ROOT.TH1D("delta_eff", "delta eff", 50, -1e-2, 1e-2)
    h2d_delta_eff = ROOT.TH2D("delta_eff_2d", "delta eff", nbins_pt, bins_pt, nbins_eta, bins_eta)
    h_delta_deff = ROOT.TH1D("delta_error_eff", "delta error eff", 50, -1e-3, 1e-2)
    h2d_delta_deff = ROOT.TH2D("delta_error_eff_2d", "delta error eff", nbins_pt, bins_pt, nbins_eta, bins_eta)
    h_pull = ROOT.TH1D("pull_eff", "pull eff", 50, -2, 2)
    h2d_pull = ROOT.TH2D("pull_eff_2d", "pull eff", nbins_pt, bins_pt, nbins_eta, bins_eta)


    for bin_key in bmark_dict.keys():

        _, bin_pt, bin_eta = bin_dict[bin_key]

        print(bin_key, bin_pt, bin_eta)

        # Bin transformation needed in case the bins are merged
        if type(bin_pt) is list:
            bin_pt = int(1+(nbins_pt*(bin_pt[0]-1)/15.))
        if type(bin_eta) is list:
            bin_eta = int(1+(nbins_eta*(bin_eta[0]-1)/48.))

        eff_1, deff_1 = res_benchmark.getEff(bin_key)
        eff_2, deff_2 = res_new.getEff(bin_key)

        delta_eff = eff_2-eff_1
        delta_deff = deff_2-deff_1

        h_delta_eff.Fill(delta_eff)
        h2d_delta_eff.SetBinContent(bin_pt, bin_eta, delta_eff)

        h_delta_deff.Fill(delta_deff)
        h2d_delta_deff.SetBinContent(bin_pt, bin_eta, delta_deff)

        h_pull.Fill(delta_eff/deff_2)
        h2d_pull.SetBinContent(bin_pt, bin_eta, delta_eff/deff_1)

    file_out.cd()
    h_delta_eff.Write()
    h2d_delta_eff.Write()
    h_delta_deff.Write()
    h2d_delta_deff.Write()
    h_pull.Write()
    h2d_pull.Write()

    file_out.Close()
    
        
###############################################################################
###############################################################################

if __name__ == '__main__':


    res_1 = results_manager("indep")
    res_2 = results_manager("sim")
    res_1.Open("results/benchmark_iso/results_iso_indep_benchmark.pkl")
    res_2.Open("results/benchmark_iso_sim/results_iso_sim.pkl")



    # res.write("results/benchmark_iso/new_results.pkl")


    compare_analysis(res_1, res_2)
    
