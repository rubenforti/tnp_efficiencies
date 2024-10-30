"""
"""
import ROOT
import sys
import statistics
import argparse
from utilities.results_manager import results_manager
from utilities.binning_utils import bin_dictionary
from utilities.base_lib import default_binnings_by_eff, binnings


stats = {
    "delta_eff" : [],
    "delta_err_eff" : [],
    "pull" : [],
    "pull_ref" : [],
    "rm1" : [],
    "ratio_error" : [],
}

stats_cut = {
    "delta_eff" : 0.1,
    "delta_err_eff" : 0.1,
    "pull" : 5,
    "pull_ref" : 5,
    "rm1" : 0.1,
    "ratio_error" : 1.5,
}


def make_stat(stat_key, eff_test=None, eff_bmark=None, d_eff_test=None, d_eff_bmark=None):

    if stat_key == "delta_eff":
        return eff_test - eff_bmark
    elif stat_key == "delta_err_eff":
        return d_eff_test - d_eff_bmark
    elif stat_key == "pull":
        return (eff_test - eff_bmark)/((d_eff_bmark**2 + d_eff_test**2)**0.5)
    elif stat_key == "pull_ref":
        return (eff_test - eff_bmark)/d_eff_bmark
    elif stat_key == "rm1":
        return (eff_test/eff_bmark) - 1
    elif stat_key == "ratio_error":
        return (d_eff_test/d_eff_bmark) - 1
    else:
        sys.exit("ERROR: stat_key not recognized")
    


parser = argparse.ArgumentParser(description='Evaluate some statistics')

parser.add_argument('-i', '--input', type=str, help='Input file', required=True)
parser.add_argument('-b', '--benchmark', type=str, help='Benchmark file')
parser.add_argument('-e', '--eff', type=str, help='Efficiency type', required=True)
parser.add_argument('--binning_pt', type=str, help='Binning in pt', default="")
parser.add_argument('--binning_eta', type=str, help='Binning in eta', default="")
parser.add_argument('--altModel_txt', type=str, default="", choices=["", "altSig", "altBkg"],
                    help='Benchmark efficiencies (in a .txt file) are selected with the specified method')
parser.add_argument('--cmpMC', action='store_true', help='Compare with MC efficiencies')
parser.add_argument('--exclude_bins', type=str, nargs="+", default=[], help='List of bins to exclude')
parser.add_argument('--ptBins', type=int, nargs=2, default=[0,-1], help='Select a range of pt bins')
parser.add_argument('--etaBins', type=int, nargs=2, default=[0,-1], help='Select a range of eta bins')

args = parser.parse_args()

if args.binning_pt=="": args.binning_pt = default_binnings_by_eff[args.eff][0]
if args.binning_eta=="": args.binning_eta = default_binnings_by_eff[args.eff][1]

if args.ptBins[1]==-1: args.ptBins[1] = (len(binnings[args.binning_pt])-1) - 1
if args.etaBins[1]==-1: args.etaBins[1] = (len(binnings[args.binning_eta])-1) - 1

res_test = results_manager(args.input, "indep", args.binning_pt, args.binning_eta)
if not args.cmpMC:
    res_bmark = results_manager(args.benchmark, "indep", args.binning_pt, args.binning_eta,
                                    altModel_check=args.altModel_txt)
else:
    res_bmark = None

for bin_key, [_, idx_pt, idx_eta] in bin_dictionary(args.binning_pt, args.binning_eta).items():

    if bin_key in args.exclude_bins: continue
    if idx_pt-1 < args.ptBins[0] or idx_pt-1 > args.ptBins[1]: continue
    if idx_eta-1 < args.etaBins[0] or idx_eta-1 > args.etaBins[1]: continue

    eff, d_eff = res_test.getEff(bin_key)
    eff_bmark, d_eff_bmark = res_bmark.getEff(bin_key) if not args.cmpMC else res_test.getEffMC(bin_key)

    for stat_key in stats.keys():
        out_stat = make_stat(stat_key, eff, eff_bmark, d_eff, d_eff_bmark)
        if abs(out_stat) < stats_cut[stat_key]:
            stats[stat_key].append(out_stat)
    

test_name_print = args.input.split(f"/")[-1]
bmark_name_print = args.benchmark.split(f"/")[-1] if args.benchmark is not None else "MC"

print(f"\n\nStats {test_name_print} vs {bmark_name_print}: \n")
print("----"*20)
print("Num bins = ", len(stats["delta_eff"]))
print("----"*20)
print("----"*20)
for stat in stats.keys():
    mean = sum(stats[stat])/len(stats[stat])
    std = (sum([(val - mean)**2 for val in stats[stat]])/len(stats[stat]))**0.5
    med_val = statistics.median(stats[stat])
    mad = statistics.median([abs(val - med_val) for val in stats[stat]])

    print(f"{stat} :  Mean={mean:.5f},  Std={std:.5f},  Median={med_val:.5f},  MAD={mad:.5f}")
    print("----"*20)
print("\n\n")
    

    
    









