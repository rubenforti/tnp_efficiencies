"""
"""
import ROOT
import statistics
from utilities.res_tools import results_manager
from utilities.dataset_utils import import_pdf_library
from utilities.base_lib import binning, bin_dictionary, eval_efficiency, sumw2_error


base_folder = "/scratch/rforti/tnp_efficiencies_results/tracking"

ws_filename = base_folder+"/legacy_fit_onlyFailSA_allMC/ws_tracking.root"

# benchmark_filename = base_folder+"/legacy_fit_onlyFailSA_allMC/ws_tracking.root"

import_pdf_library("RooCMSShape")

file = ROOT.TFile(ws_filename, "READ")
ws = file.Get("w")
results = results_manager("indep", "pt_tracking", "eta", import_ws=ws)


benchmark_filename = base_folder+"/../egm_tnp_results/tracking/allEfficiencies.txt"
with open(benchmark_filename, "r") as file_bmark:
    row_list = file_bmark.readlines()
    print(type(row_list))
    res_benchmark = results_manager("indep", "pt_tracking", "eta", import_txt=row_list)

'''
file_bmark = ROOT.TFile(benchmark_filename, "READ")
ws_bmark = file_bmark.Get("w")
res_benchmark = results_manager("indep", "pt_tracking", "eta", import_ws=ws_bmark)
'''

bin_dict = bin_dictionary("pt_tracking", "eta")


stats = {
    "delta_eff" : [],
    "delta_err_eff" : [],
    "pull" : [],
    "rm1" : [],
    "ratio_error" : [],
}


exclude_bins = ["[55.0to65.0][-0.3to-0.2]", "[55.0to65.0][-0.2to-0.1]", 
                "[55.0to65.0][0.2to0.3]", "[55.0to65.0][0.5to0.6]",   "[55.0to65.0][0.7to0.8]"]


for bin_key in bin_dict.keys():

    
    # if bin_key in exclude_bins: continue

    eff, d_eff = results.getEff(bin_key)
    eff_benchmark, d_eff_benchmark = res_benchmark.getEff(bin_key)

    #if abs((eff - eff_benchmark)/d_eff_benchmark) > 2 : continue
    
    stats["delta_eff"].append(eff - eff_benchmark)
    stats["delta_err_eff"].append(d_eff - d_eff_benchmark)
    stats["pull"].append((eff - eff_benchmark)/d_eff_benchmark)
    stats["rm1"].append((eff/eff_benchmark) - 1)
    stats["ratio_error"].append((d_eff/d_eff_benchmark) - 1)


    '''
    eff_mc, d_eff_mc = eval_efficiency(ws.data(f"Minv_mc_pass_{bin_key}").sumEntries(), 
                                       ws.data(f"Minv_mc_fail_{bin_key}").sumEntries(),
                                       sumw2_error(ws.data(f"Minv_mc_pass_{bin_key}")),
                                       sumw2_error(ws.data(f"Minv_mc_fail_{bin_key}")))


    stats["delta_eff"].append(eff - eff_mc)
    stats["delta_err_eff"].append(d_eff - eff_mc)
    stats["pull"].append((eff - eff_mc)/d_eff)
    stats["rm1"].append((eff/eff_mc) - 1)
    stats["ratio_error"].append((d_eff/d_eff_mc) - 1)

    # print(bin_key, (eff - eff_mc)/d_eff)
    '''




name_print = ws_filename.split("tracking/")[1]

print(f"\nStats {name_print}: \n")

print("Num bins = ", len(stats["delta_eff"]))
print("")
for stat in stats.keys():
    mean = sum(stats[stat])/len(stats[stat])
    std = (sum([(val - mean)**2 for val in stats[stat]])/len(stats[stat]))**0.5
    med_val = statistics.median(stats[stat])
    mad = statistics.median([abs(val - med_val) for val in stats[stat]])

    print(f"{stat} :  Mean={mean},  Std={std},  Median={med_val},  MAD={mad}\n")

print("\n\n")
    

    
    









