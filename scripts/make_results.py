"""
"""
import os
import sys
import math
import ROOT
import argparse
from array import array
from copy import copy
from utilities.base_lib import default_binnings_by_eff, safeSystem
from utilities.binning_utils import bin_dictionary
from utilities.results_manager import results_manager
from utilities.results_utils import init_histos, fill_res_histograms, fill_resCmp_histograms


NBINS = 35

eff_min = 0.895
rel_err_eff_min, rel_err_eff_max = 1e-4, 7e-3
sf_min, sf_max =  0.99, 1.03

delta_min = -2e-3
delta_error_min = -2e-3
pull_min = -1
rm1_min = -1e-3
ratio_error_min = -1


res_dict = {
    "efficiency" : { 
        "title" : "Efficiency", 
        "array" : array("d", [round(eff_min + (1-eff_min)*(i/NBINS), 4) for i in range(NBINS+1)])},
    "rel_err_efficiency" : {
        "title" : "Relative error on efficiency", 
        "array" : array("d", [round(rel_err_eff_min + (rel_err_eff_max-rel_err_eff_min)*(i/NBINS), 6) for i in range(NBINS+1)])},
    "efficiency_MC" : {
        "title" : "Efficiency MC", 
        "array" : array("d", [round(eff_min + (1-eff_min)*(i/NBINS), 4) for i in range(NBINS+1)])},
    "sf" : {
        "title" : "Scale Factor", 
        "array" : array("d", [round(sf_min + (sf_max-sf_min)*(i/NBINS), 5) for i in range(NBINS+1)])},
}
resCmp_dict = {
    "delta" : {
        "title" : "Delta efficiency", 
        "array" : array("d", [round(delta_min + (-2*delta_min/NBINS)*i, 6) for i in range(NBINS+1)])},
    "delta_err" : {
        "title" : "Delta error on efficiency", 
        "array" : array("d", [round(delta_error_min + (-2*delta_error_min/NBINS)*i, 6) for i in range(NBINS+1)])},
    "pull" : {
        "title" : "Pull", 
        "array" : array("d", [round(pull_min + (-2*pull_min/NBINS)*i, 6) for i in range(NBINS+1)])},
    "pull_ref" : {
        "title" : "Pull (referred to bmark error)", 
        "array" : array("d", [round(pull_min + (-2*pull_min/NBINS)*i, 6) for i in range(NBINS+1)])},
    "rm1" : {
        "title" : "Relative bias", 
        "array" : array("d", [round(rm1_min + (-2*rm1_min/NBINS)*i, 6) for i in range(NBINS+1)])},
    "ratio_err" : {
        "title" : "Ratio of errors minus 1", 
        "array" : array("d", [round(ratio_error_min + (-2*ratio_error_min/NBINS)*i, 6) for i in range(NBINS+1)])}
}

def zeroBkg_estimate(pars):
        for par in pars:
            if "nbkg" in par.GetName() and "fail" in par.GetName():
                if math.isclose(par.getVal(), 0.5, abs_tol=1.5): #
                    # if the fit estimates less than 2 bkg events, it's considered as zero-background
                    return True
                else:
                    return False
        return True

###############################################################################


parser = argparse.ArgumentParser(description='Produce results and compare with a benchmark')

subparsers = parser.add_subparsers(help='Sub-command help')

parser.add_argument('-i', '--test_file', type=str, help='Input file (results to be tested)', required=True)
parser.add_argument('-e', '--eff', type=str, help='Efficiency type', required=True)
parser.add_argument('--altModel_txt', type=str, default="", choices=["", "altSig", "altBkg"],
                    help='Benchmark efficiencies (in a .txt file) are selected with the specified method')
parser.add_argument('--binning_pt', type=str, help='Binning in pt', default="")
parser.add_argument('--binning_eta', type=str, help='Binning in eta', default="")
parser.add_argument('--postfix', type=str, default="", help='Postfix name for output file')
parser.add_argument('--type_analysis', type=str, default="indep", choices=["indep", "sim"], help='Type of analysis')
parser.add_argument('--pseudodata', action='store_true', help='Evaluation on pseudodata')
parser.add_argument('--figs', action='store_true', help='Produce figures')
parser.add_argument('--setLog', action='store_true', help='Set log scale on y axis (when plotting)')

parser_cmp = subparsers.add_parser('cmp', help='Compare results')
parser_cmp.add_argument('-b', '--bmark_file', type=str, help='Benchmark file', required=True)
parser_cmp.add_argument('-c', '--comparisons', type=str, nargs="+", default=list(resCmp_dict.keys()),
                        help='List of comparison statistics to be evaluated')
parser_cmp.add_argument('--projection', action='store_true', help='Save projections on pt and eta')
parser_cmp.add_argument('--isolate_effect', type=str, default="", choices=["pass", "fail"], 
                        help='Isolate effect on pass or fail')
parser_cmp.add_argument('--eval_onlyFitsWithBkg', action='store_true',
                         help='Make comparisons only on bins with Bkg (selection only on fail)')
parser_cmp.add_argument('--eval_asymm', action='store_true', help='Evaluate asymmetry in MINOS errors')
parser_cmp.add_argument('--add_plain_res', action='store_true', help='Add plain results in the root files')

args = parser.parse_args()

# Selection of which histograms to produce
if not hasattr(args, 'bmark_file'):
    list_histograms = res_dict if not args.pseudodata else resCmp_dict
elif not args.add_plain_res:
    list_histograms = {k:v for k,v in resCmp_dict.items() if k in args.comparisons}
else:
    list_histograms = {res_dict.items(), resCmp_dict.items()}

if args.binning_pt=="": args.binning_pt = default_binnings_by_eff[args.eff][0]
if args.binning_eta=="": args.binning_eta = default_binnings_by_eff[args.eff][1]

bin_dict = bin_dictionary(args.binning_pt, args.binning_eta)

binnings_1d = {k : v["array"] for k, v in list_histograms.items()}
names_titles_dict = {k : v["title"] for k, v in list_histograms.items()}

histos = init_histos(names_titles_dict.keys(), names_titles_dict.values(), 
                     binning_var=binnings_1d, binning_pt=args.binning_pt, binning_eta=args.binning_eta, flag="")


res_test = results_manager(args.test_file, args.type_analysis, args.binning_pt, args.binning_eta, 
                          altModel_check=args.altModel_txt)
res_benchmark = None


if hasattr(args, 'bmark_file'):
    res_benchmark = results_manager(args.bmark_file, args.type_analysis, args.binning_pt, args.binning_eta, 
                                    altModel_check=args.altModel_txt)

    if args.eval_onlyFitsWithBkg:
        bin_dict_original = copy(bin_dict)
        
        for b_key in bin_dict_original.keys():
            pars = res_test.getPars("fail", b_key) if args.type_analyisis=="indep" else res_test.getPars("sim", b_key)
            if zeroBkg_estimate(pars): bin_dict.pop(b_key)






fill_res_histograms(res_test, histos, bin_dict)
fill_resCmp_histograms(res_benchmark, res_test, histos, bin_dict, isPseudodata=args.pseudodata) 

'''
histos_copy = histos.copy()
for hist_key in histos_copy.keys():
    if "2d" in hist_key:
        histos[hist_key].Sumw2()
        histo_pt = histos[hist_key].ProjectionX(hist_key.replace("2d", "pt"), 1, len(bins_eta)-1)
        histo_pt.Scale(1/(len(bins_eta)-1))
        histo_eta = histos[hist_key].ProjectionY(hist_key.replace("2d", "eta"), 1, len(bins_pt)-1)
        histo_eta.Scale(1/(len(bins_pt)-1))
        histos.update({histo_pt.GetName() : histo_pt, histo_eta.GetName() : histo_eta})
'''

base_lib, outname = args.test_file.rsplit("/", 1)

if hasattr(args, 'bmark_file'):
    postfix = "-"+args.postfix if args.postfix!="" else ""
    if args.bmark_file.endswith(".txt"):
        add_flag = f"_cmp-egm-{args.altModel_txt}" if args.altModel_txt!="" else "_cmp-egm"
    else:
        add_flag = "_cmp-" + args.bmark_file.split("_")[-1].replace(".root", "") if len(args.bmark_file.split("_"))>2 else "" 

    add_flag += postfix

    if args.eval_onlyFitsWithBkg: 
        add_flag += "-onlyFitsWithBkg"

else:
    add_flag = "_"+args.postfix if args.postfix!="" else ""

if args.test_file.endswith(".root"):
    gen_name = (add_flag+"_"+outname.split("/")[-1].split("_", 1)[1]).replace(".root", "")
else:
    gen_name = outname.split("/")[-1].replace(".txt", "")

if args.figs:
    new_outFolder = base_lib+"/results"+gen_name
    if os.path.exists(new_outFolder):
        check = input(f"Directory {new_outFolder} already exists. Do you want to overwrite it? [y/n] ")
        if check=="y":
            safeSystem(f"rm -rf {new_outFolder}")
        else:
            sys.exit(0)
    os.mkdir(new_outFolder) 
    base_lib = new_outFolder       

outname = "resHist"+gen_name+".root"

file_out = ROOT.TFile(base_lib+"/"+outname, "RECREATE")
file_out.cd()
[histo.Write() for histo in histos.values()]
file_out.Close()

print("\n\nSaved file with results:", file_out.GetName())

if args.figs:
    cmd = f"python scripts/plotting/plot_hres.py -i {base_lib}/{outname} -e {args.eff}"
    if args.setLog: 
        cmd += " --log"
    safeSystem(cmd)