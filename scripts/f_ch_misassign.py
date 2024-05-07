import ROOT
from utilities.base_lib import bin_dictionary, binning, eval_efficiency
from copy import copy
import sys

import numpy as np
import matplotlib.pyplot as plt


histos_OS = {
    "pass": {
        "plus" : 0,
        "minus" : 0
    },
    "fail": {
        "plus" : 0,
        "minus" : 0
    }
} 
histos_SS = {
    "pass": {
        "plus" : 0,
        "minus" : 0
    },
    "fail": {
        "plus" : 0,
        "minus" : 0
    }
}



res_file = ROOT.TFile.Open("/scratch/rforti/tnp_efficiencies_results/tracking/benchmark_opt_with_legacy_checks/res_tracking_indep_benchmark.root")

h_rel_errors = res_file.Get("h_rel_err_efficiency_2d")


for ch in ["plus", "minus"]:

    file_OS = ROOT.TFile.Open(f"datasets/OS/tnp_tracking{ch}_mc_vertexWeights1_oscharge1.root")
    file_SS = ROOT.TFile.Open(f"datasets/SS/tnp_tracking{ch}_mc_vertexWeights1_oscharge0.root")

    for flag in ["pass", "fail"]:
        
        th3_tmp_OS = file_OS.Get(f"{flag}_mu_DY_postVFP")
        th3_tmp_SS = file_SS.Get(f"{flag}_mu_DY_postVFP")

        histos_OS[flag][ch] = copy(th3_tmp_OS)
        histos_SS[flag][ch] = copy(th3_tmp_SS)


ch_misassign_pass = ROOT.TH2D("ch_misassign_pass", "ch_misassign_pass",
                              4, binning("pt_tracking"), 48, binning("eta"))
ch_misassign_fail = ROOT.TH2D("ch_misassign_fail", "ch_misassign_fail",
                              4, binning("pt_tracking"), 48, binning("eta"))
ch_misassign_ratio = ROOT.TH2D("ch_misassign_ratio", "ch_misassign_ratio",
                              4, binning("pt_tracking"), 48, binning("eta"))


eff_rel_bias_OS_SS = ROOT.TH2D("eff_rel_bias_OS_SS", "eff_rel_bias_OS_SS",
                                4, binning("pt_tracking"), 48, binning("eta"))

eff_pull_OS_SS = ROOT.TH2D("eff_pull_OS_SS", "eff_pull_OS_SS",
                        4, binning("pt_tracking"), 48, binning("eta"))

delta_rel_vs_OS = ROOT.TH2D("delta_rel_vs_OS", "delta_rel_vs_OS",
                        4, binning("pt_tracking"), 48, binning("eta"))

h_delta_rel_vs_OS = ROOT.TH1D("h_delta_rel_vs_OS", "h_delta_rel_vs_OS",
                              50, -1.7, 0.3)

delta_rel_vs_SS = ROOT.TH2D("delta_rel_vs_SS", "delta_rel_vs_SS",
                        4, binning("pt_tracking"), 48, binning("eta"))


correlation_dRel_fCh = ROOT.TH2D("correlation_dRel_fCh", "correlation_dRel_fCh",
                                 50, 0, 0.27, 50, -1.7, 0.2)

manual_bias_test = ROOT.TH2D("manual_bias_test", "manual_bias_test",
                             4, binning("pt_tracking"), 48, binning("eta"))

delta_rel = []


for b_pt in range(1, 5):
    for b_eta in range(1, 49):

        n_pass_OS, n_fail_OS, n_pass_SS, n_fail_SS = 0, 0, 0, 0
        err_n_pass_OS, err_n_fail_OS, err_n_pass_SS, err_n_fail_SS = 0, 0, 0, 0
        

        for ch in ["plus", "minus"]:

            h_pass_OS = histos_OS["pass"][ch].ProjectionX("h_pass_OS", b_pt, b_pt, b_eta, b_eta, "e")
            h_fail_OS = histos_OS["fail"][ch].ProjectionX("h_fail_OS", b_pt, b_pt, b_eta, b_eta, "e")
            h_pass_SS = histos_SS["pass"][ch].ProjectionX("h_pass_SS", b_pt, b_pt, b_eta, b_eta, "e")
            h_fail_SS = histos_SS["fail"][ch].ProjectionX("h_fail_SS", b_pt, b_pt, b_eta, b_eta, "e")

            n_pass_OS += h_pass_OS.Integral()
            n_fail_OS += h_fail_OS.Integral()
            n_pass_SS += h_pass_SS.Integral()
            n_fail_SS += h_fail_SS.Integral()

            for i in range (1, 81):
                err_n_pass_OS += h_pass_OS.GetBinError(i)**2
                err_n_fail_OS += h_fail_OS.GetBinError(i)**2
                err_n_pass_SS += h_pass_SS.GetBinError(i)**2
                err_n_fail_SS += h_fail_SS.GetBinError(i)**2

        err_n_pass = err_n_pass_OS + err_n_pass_SS
        err_n_fail = err_n_fail_OS + err_n_fail_SS

        err_n_pass = err_n_pass**0.5
        err_n_fail = err_n_fail**0.5
        err_n_pass_OS = err_n_pass_OS**0.5
        err_n_fail_OS = err_n_fail_OS**0.5
        err_n_pass_SS = err_n_pass_SS**0.5
        err_n_fail_SS = err_n_fail_SS**0.5

        eff_legacy, d_eff_legacy = eval_efficiency(
            n_pass_SS+n_pass_OS, n_fail_SS+n_fail_OS, err_n_pass, err_n_fail)
        
        eff_OS, d_eff_OS = eval_efficiency(
            n_pass_OS, n_fail_OS, err_n_pass_OS, err_n_fail_OS)
        
        eff_SS, d_eff_SS = eval_efficiency(
            n_pass_SS, n_fail_SS, err_n_pass_SS, err_n_fail_SS)
        
        f_pass = n_pass_SS/(n_pass_SS+n_pass_OS)
        f_fail = n_fail_SS/(n_fail_SS+n_fail_OS)

        ch_misassign_pass.SetBinContent(b_pt, b_eta, f_pass)
        ch_misassign_fail.SetBinContent(b_pt, b_eta, f_fail)
        ch_misassign_ratio.SetBinContent(b_pt, b_eta, f_fail/f_pass)

        eff_rel_bias_OS_SS.SetBinContent(b_pt, b_eta, (eff_SS-eff_OS)/eff_OS)
        eff_pull_OS_SS.SetBinContent(b_pt, b_eta, (eff_OS-eff_SS)/(d_eff_OS**2+d_eff_SS**2)**0.5)
        
        rel_error = h_rel_errors.GetBinContent(b_pt, b_eta)
        d_OS = (eff_legacy - eff_OS)/(rel_error*eff_OS)
        d_SS = (eff_legacy - eff_SS)/(rel_error*eff_SS)

        delta_rel_vs_OS.SetBinContent(b_pt, b_eta, d_OS)
        delta_rel_vs_SS.SetBinContent(b_pt, b_eta, d_SS) 

        h_delta_rel_vs_OS.Fill(d_OS)

        correlation_dRel_fCh.Fill(f_fail, d_OS)

        #print(round(eff_legacy, 3), "", round(eff_OS,3), round(d_OS,1), "", round(eff_SS,3), round(d_SS,1))
        delta_rel.append(d_OS)

        bias_factor_fail = np.linspace((1-2*f_fail), 1, 75) #factor to be applied on the number of failing events in the evaluation of the OS+SS efficiency

        biased_eff, biased_eff_error = [], []
        for i in range(75):
            eff_leg_biased, err_eff_biased = eval_efficiency(
                n_pass_SS+n_pass_OS, (n_fail_SS+n_fail_OS)*bias_factor_fail[i], err_n_pass, err_n_fail*bias_factor_fail[i])
            biased_eff.append(eff_leg_biased)
            biased_eff_error.append(err_eff_biased)
        
        biased_eff = np.array(biased_eff) 
        biased_eff_error = np.array(biased_eff_error)
        #delta_rel_biased = (biased_eff-np.ones(75)*eff_OS)/(np.ones(75)*(eff_OS*rel_error))



        plt.figure()
        plt.xlabel("Bias factor applied to nFail")
        plt.ylabel("Efficiency")
        plt.errorbar(bias_factor_fail, biased_eff, biased_eff_error, np.zeros(75), linestyle="", marker=".", label="Biased efficiency")
        plt.plot(bias_factor_fail, np.ones(75)*eff_legacy, linestyle=":", marker="", label="Legacy")
        plt.plot(bias_factor_fail, np.ones(75)*eff_OS*(1+rel_error), linestyle="-", marker="", color="red")
        plt.plot(bias_factor_fail, np.ones(75)*eff_OS*(1-rel_error), linestyle="-", marker="", color="red")
        plt.plot(bias_factor_fail, np.ones(75)*eff_OS, linestyle="--", marker="", color="red", label="OS")
        plt.axvline(1-f_fail, linestyle="-", color="green", label="Expected fCH fail")
        plt.legend()
        plt.show()
        plt.savefig(f"bias_mc_plots/bias_{b_pt}_{b_eta}.png")
        plt.close()

        test_bias_up = biased_eff[biased_eff>(eff_OS*(1+rel_error))]
        test_bias_down = biased_eff[biased_eff<(eff_OS*(1-rel_error))]
        print(len(test_bias_up), len(test_bias_down))

        manual_bias_test.SetBinContent(b_pt, b_eta, len(test_bias_up) + len(test_bias_down))






delta_rel = np.array(delta_rel)
print("mean", np.mean(delta_rel))
print("std", np.std(delta_rel))
print("median", np.median(delta_rel))
print("max", np.max(delta_rel))
print("min", np.min(delta_rel))



file_out = ROOT.TFile.Open("ch_misassign.root", "recreate")
file_out.cd()
ch_misassign_pass.Write()
ch_misassign_fail.Write()
ch_misassign_ratio.Write()
eff_rel_bias_OS_SS.Write()
eff_pull_OS_SS.Write()
delta_rel_vs_OS.Write()
h_delta_rel_vs_OS.Write()
delta_rel_vs_SS.Write()
correlation_dRel_fCh.Write()
manual_bias_test.Write()
file_out.Close()




        
