import ROOT
from utilities.base_library import bin_dictionary
import numpy as np
from scipy.stats import ks_2samp as KS
import matplotlib.pyplot as plt
from statistics import median
import sys

file = ROOT.TFile("/scratch/rforti/tnp_efficiencies_results/tracking/bkg_figs/ws_tracking_bkg.root")
ws = file.Get("w")

file_SS = ROOT.TFile("/scratch/rforti/tnp_efficiencies_results/tracking/bkg_SS_figs/ws_tracking_bkg_SS.root")
ws_SS = file_SS.Get("w")


bkg_categories = [# "bkg_WW", "bkg_WZ", "bkg_ZZ", 
                  # "bkg_TTSemileptonic", "bkg_TTFullyleptonic", "bkg_Ztautau",
                  "bkg_WplusJets", "bkg_WminusJets"]
                  
bin_dict = bin_dictionary("pt_tracking", "eta")

ratio_samecharge_mcbkg, pval_KS_samecharge = [], []
ratio_SS_OS_wjets, pval_KS_wjets_SS = [], []

ch = ["plus", "minus"]


for b_key in bin_dict.keys():
    


    nbkg_mc_OS, nbkg_mc_SS = 0, 0
    
    array_mc_OS, array_mc_SS, array_samecharge = [0]*80, [0]*80, [0]*80

    errors_mc_OS, errors_mc_SS = [0]*80, [0]*80

    
    for cat in bkg_categories: 
        data = ws.data(f"Minv_bkg_fail_{b_key}_{cat.replace('bkg_', '')}")
        nbkg_mc_OS += data.sumEntries()
        for i in range(80): 
            data.get(i)
            array_mc_OS[i] += data.weight(i)
            errors_mc_OS[i] += data.weightError(ROOT.RooAbsData.SumW2)**2

        data_SS = ws_SS.data(f"Minv_bkg_fail_{b_key}_{cat.replace('bkg_', '')}_SS")
        nbkg_mc_SS += data_SS.sumEntries()
        for i in range(80): 
            data_SS.get(i)
            array_mc_SS[i] += data_SS.weight(i)
            errors_mc_SS[i] += data_SS.weightError(ROOT.RooAbsData.SumW2)**2


    data_samecharge = ws.data(f"Minv_bkg_fail_{b_key}_SameCharge")
    nbkg_ss = data_samecharge.sumEntries()
    for i in range(80): array_samecharge[i] += data_samecharge.weight(i)
    
    ratio_samecharge_mcbkg.append(nbkg_ss/nbkg_mc_OS)
    ratio_SS_OS_wjets.append(nbkg_mc_SS/nbkg_mc_OS)

    distr_ratio_SS_OS = np.divide(array_mc_SS, array_mc_OS)

    errors_mc_OS = np.array(errors_mc_OS)**0.5
    errors_mc_SS = np.array(errors_mc_SS)**0.5
    array_mc_OS = np.array(array_mc_OS)
    array_mc_SS = np.array(array_mc_SS)

    error_distr_ratio_SS_OS = abs(distr_ratio_SS_OS)*((np.divide(errors_mc_SS, array_mc_OS)**2 + np.divide(errors_mc_OS, array_mc_OS)**2)**0.5)
    # error_distr_ratio_SS_OS = [val for val in error_distr_ratio_SS_OS if val!=np.inf]

    for i in range(80):
        err_i = error_distr_ratio_SS_OS[i]
        val_i = distr_ratio_SS_OS[i]
        if err_i == np.inf or err_i > 3*val_i or val_i == 0: 
            error_distr_ratio_SS_OS[i] = 0
            #distr_ratio_SS_OS[i] = 0

    med_val = median([val for val in distr_ratio_SS_OS if val!=np.inf])
    mad = median([abs(val - med_val) for val in distr_ratio_SS_OS if val!=np.inf])
    
    
    array_mc_OS = np.array(array_mc_OS)/nbkg_mc_OS
    array_mc_SS = np.array(array_mc_SS)/nbkg_mc_SS
    array_samecharge = np.array(array_samecharge)/nbkg_ss
    
    pval_KS_samecharge.append(KS(array_mc_OS, array_samecharge, alternative="twosided").pvalue)
    pval_KS_wjets_SS.append(KS(array_mc_OS, array_mc_SS, alternative="twosided").pvalue)
    
    if(pval_KS_samecharge[-1] < 1./192) or (pval_KS_wjets_SS[-1] < 1./192): print(b_key)
    axis = np.linspace(50.5, 129.5, 80)
    plt.figure()
    plt.subplot(2,1,1)
    plt.plot(axis, array_mc_OS, label="OS")
    plt.plot(axis, array_mc_SS, label="SS")
    plt.legend()
    plt.subplot(2,1,2)
    plt.errorbar(axis, distr_ratio_SS_OS-1, error_distr_ratio_SS_OS, linestyle="", marker=".", color="blue")
    plt.plot(axis, np.zeros(80), linestyle="--", color="red")
    plt.show()
    plt.savefig(f"wjets_figs/{b_key}_wjets_ratio.pdf")
    

    

ratio_samecharge_mcbkg = np.array(ratio_samecharge_mcbkg)
pval_KS_samecharge = np.array(pval_KS_samecharge)
ratio_SS_OS_wjets = np.array(ratio_SS_OS_wjets)
pval_KS_wjets_SS = np.array(pval_KS_wjets_SS)


print("Ratio SameCharge/mcBkg:", ratio_samecharge_mcbkg.mean(), ratio_samecharge_mcbkg.max(), ratio_samecharge_mcbkg.min())
print("Pval KS_test SameCharge/mcBkg (OS):", pval_KS_samecharge.mean(), pval_KS_samecharge.max(), pval_KS_samecharge.min())
print("Ratio SS/OS mcBkg:", ratio_SS_OS_wjets.mean(), ratio_SS_OS_wjets.max(), ratio_SS_OS_wjets.min()) 
print("Pval KS_test SS/OS mcBkg:", pval_KS_wjets_SS.mean(), pval_KS_wjets_SS.max(), pval_KS_wjets_SS.min())
    


