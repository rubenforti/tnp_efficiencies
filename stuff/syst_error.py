"""
"""

import ROOT
from utilities.results_utils import compare_eff_results
from utilities.base_library import binning


file_benchmark = ROOT.TFile("results/benchmark_iso/ws_iso_indep_benchmark.root", "READ")
ws_benchmark = file_benchmark.Get("w")

file_test = ROOT.TFile("root_files/ws_iso_indep_mcbkg.root", "READ")
ws_test = file_test.Get("w")

file_test_aux = ROOT.TFile("root_files/ws_iso_indep_mcbkg_mergedbins.root", "READ")
ws_test_aux = file_test_aux.Get("w")

prob_bins_equivalence = {"[28.0to30.0][2.0to2.1]" : "[28.0to30.0][1.8to2.4]",
                         "[28.0to30.0][2.1to2.2]" : "[28.0to30.0][1.8to2.4]",
                         "[28.0to30.0][2.3to2.4]" : "[28.0to30.0][1.8to2.4]", 
                         "[32.0to34.0][2.3to2.4]" : "[32.0to34.0][1.8to2.4]",  
                         "[32.0to34.0][-2.4to-2.3]" : "[32.0to34.0][-2.4to-1.8]", 
                         "[55.0to60.0][-2.4to-2.3]" : "[55.0to60.0][-2.4to-1.8]", 
                         "[60.0to65.0][-1.9to-1.8]" : "[60.0to65.0][-2.4to-1.8]"
                        }

compare_eff_results(ws_benchmark, ws_test, "pt", "eta", "systematic_error_eff.root", 
                    aux_dict=prob_bins_equivalence, aux_filename="root_files/ws_iso_indep_mcbkg_mergedbins.root")

