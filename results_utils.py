"""
"""
import ROOT
import os
import pickle
from array import array
from utilities import fit_quality, eval_efficiency, pearson_chi2_eval


class results_manager:
    """
    """

    def __init__(self, type_estimate):
        """
        """
        self._dict_results = {}
        self._analysis = type_estimate
        # self._ws = workspace

    def Open(self, filename):
        with open(filename, "rb") as file:
            self._dict_results = pickle.load(file)

    def Write(self, filename):
        with open(filename, "wb") as file:
            pickle.dump(self._dict_results, file)
            # file.close()
    

    def dictionary(self):
        return self._dict_results

    def add_result(self, ws, bin_pt, bin_eta, check=False, old_checks=False):
        """
        """
        Npass, sigma_Npass = 0, 0
        Nfail, sigma_Nfail = 0, 0

        if check is True:
            # SISTEMARE, COSÌ NON FUNZIONA !!!!!!!!
            quality = 1
            for result in res:
                quality = quality*fit_quality(result, old_checks=old_checks)
            goodfit = bool(quality)
        else:
            goodfit = True
        

        if goodfit and (f"{bin_pt},{bin_eta}" not in self._dict_results) and self._analysis == 'indep':
            res_pass = ws.obj(f"results_pass_({bin_pt}|{bin_eta})")
            res_fail = ws.obj(f"results_fail_({bin_pt}|{bin_eta})")
            print(res_pass.GetName(), res_fail.GetName())
            for par in res_pass.floatParsFinal():
                if par.GetName() == f'nsig_pass_({bin_pt}|{bin_eta})':
                    Npass = par.getVal()
                    sigma_Npass = par.getError()
            for par in res_fail.floatParsFinal():
                if par.GetName() == f'nsig_fail_({bin_pt}|{bin_eta})':
                    Nfail = par.getVal()
                    sigma_Nfail = par.getError()

            # histo_pass = self._ws.data(f"Minv_data_pass_({bin_pt}|{bin_eta})")
            # histo_fail = self._ws.data(f"Minv_data_fail_({bin_pt}|{bin_eta})")

            new_res = {
                f"{bin_pt},{bin_eta}": {
                    "efficiency" : eval_efficiency(Npass, Nfail, sigma_Npass, sigma_Nfail),
                    "pars_pass" : res_pass.floatParsFinal(),
                    "corrmatrix_pass" : res_pass.correlationMatrix(),
                    "status_pass" : (res_pass.status(), res_pass.covQual(), res_pass.edm()),
                    "pars_fail": res_fail.floatParsFinal(),
                    "corrmatrix_fail": res_fail.correlationMatrix(),
                    "status_fail" : (res_fail.status(), res_fail.covQual(), res_fail.edm())
                    }
                }

            self._dict_results.update(new_res)

        elif goodfit and self._analysis == 'sim':
            res = ws.obj(f"results_({bin_pt}|{bin_eta})")
            new_res = {f"{bin_pt},{bin_eta}": {
                "efficiency" : (res.floatParsFinal().find(f"efficiency_({bin_pt}|{bin_eta})").getVal(),
                                res.floatParsFinal().find(f"efficiency_({bin_pt}|{bin_eta})").getError()),
                "parameters": res.floatParsFinal(),
                "corrmatrix": res.covarianceMatrix(),
                "status": (res.status(), res.covQual(), res.edm())
                }
            }
            self._dict_results.update(new_res)
        else:
            pass

    def getEfficiency(self, bin_pt=0, bin_eta=0):

        if bin_pt == 0 and bin_eta == 0:
            eff = [] 
            [eff.append(self._dict_results[key]["efficiency"]) for key in self._dict_results.keys()]
                # d_eff.append(self._dict_results[key]["efficiency"][1])
        else:
            eff = self._dict_results[f"{bin_pt},{bin_eta}"]["efficiency"]
            # d_eff = self._dict_results[f"{bin_pt},{bin_eta}"]["efficiency"][1]
        return eff

    def getCorr(self, bin_pt=0, bin_eta=0):

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
    


def compare_analysis(res_1, res_2):
    """
    Compare the efficiencies and their error between two results files. The first file is 
    considered as benchmark
    """
    

    file = ROOT.TFile("prova_histo2d.root", "RECREATE")

    
    binning_pt = array('d', [24., 26., 28., 30., 32., 34.,
                       36., 38., 40., 42., 44., 47., 50., 55., 60., 65.])
    binning_eta = array('d', [round(-2.4 + i*0.1, 2) for i in range(49)])
    
    h_delta_eff = ROOT.TH1D("delta_eff", "delta eff", 50, -1e-4, 1e-4)
    h2d_delta_eff = ROOT.TH2D("delta_eff_2d", "delta eff", len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta)
    h_delta_deff = ROOT.TH1D("delta_deff", "delta deff", 50, -1e-5, 1e-5)
    h2d_pull = ROOT.TH2D("pull_delta_eff_2d", "pull delta eff", len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta)

    for i in range(1, 2):
        for j in range(1, 49):
            eff_1, deff_1 = res_1.getEfficiency(i,j)
            eff_2, deff_2 = res_2.getEfficiency(i,j)
            
            h_delta_eff.Fill(eff_1-eff_2)
            h_delta_deff.Fill(deff_1-deff_2)

            h2d_delta_eff.SetBinContent(i, j, eff_1-eff_2)
            h2d_pull.SetBinContent(i, j, abs(eff_1-eff_2)/eff_1)

    c0 = ROOT.TCanvas("", "", 1200, 900)
    c0.Divide(2)
    c0.cd(1)
    ROOT.gStyle.SetOptStat("men")
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
    h2d_pull.Draw("COLZ")
    h2d_pull.SetContour(25)
    c2.SaveAs("figs/pull_delta_eff_2d.pdf")

    h_delta_eff.Write()
    h_delta_deff.Write()
    h2d_delta_eff.Write()
    h2d_pull.Write()

    file.Close()
    
        
    
    









if __name__ == '__main__':


    res_1 = results_manager("indep")
    res_2 = results_manager("sim")
    res_1.Open("results/benchmark_iso/new_results.pkl")
    res_2.Open("results/results_iso_sim.pkl")



    # res.write("results/benchmark_iso/new_results.pkl")


    compare_analysis(res_1, res_2)
    
