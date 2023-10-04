"""
"""
import ROOT
from statistics import median
from utilities.results_utils import results_manager
from utilities.base_library import binning, bin_dictionary


ws_filename = "results/benchmark_iso_r628/ws_iso_indep_bmark_new.root"

# benchmark_filename = "results/benchmark_iso_r628/ws_iso_indep_bmark.root"

benchmark_filename = "results/benchmark_iso/old_results.txt"


file = ROOT.TFile(ws_filename, "READ")
ws = file.Get("w")
results = results_manager("indep", "pt", "eta", import_ws=ws)


with open(benchmark_filename, "r") as file_bmark:
    row_list = file_bmark.readlines()
    print(type(row_list))
    res_benchmark = results_manager("indep", "pt", "eta", import_txt=row_list)

'''
file_bmark = ROOT.TFile(benchmark_filename, "READ")
ws_bmark = file_bmark.Get("w")
res_benchmark = results_manager("indep", "pt", "eta", import_ws=ws_bmark)
'''

bin_dict = bin_dictionary("pt", "eta")


stats = {
    "delta_eff" : [],
    "delta_err_eff" : [],
    "pull" : [],
    "rm1" : [],
    "ratio_error" : [],
}

for bin_key in bin_dict.keys():

    eff, d_eff = results.getEff(bin_key)
    eff_bmark, d_eff_bmark = res_benchmark.getEff(bin_key)

    stats["delta_eff"].append(eff - eff_bmark)
    stats["delta_err_eff"].append(d_eff - d_eff_bmark)
    stats["pull"].append((eff - eff_bmark)/d_eff_bmark)
    stats["rm1"].append((eff/eff_bmark) - 1)
    stats["ratio_error"].append((d_eff/d_eff_bmark) - 1)


for stat in stats.keys():
    med_val = median(stats[stat])
    mad = median([abs(val - med_val) for val in stats[stat]])

    print(f"{stat} :  Median={med_val},  MAD={mad}\n")
    









