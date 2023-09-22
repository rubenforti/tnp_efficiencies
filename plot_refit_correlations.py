import ROOT
from utilities.results_utils import results_manager
from utilities.base_library import binning, bin_dictionary

ws_filename = "results/iso_sim/ws_iso_sim.root"

benchmark_filename = "results/benchmark_iso/ws_iso_sim_benchmark.root"

file = ROOT.TFile(ws_filename, "READ")
ws = file.Get("w")

benchmark_file = ROOT.TFile(benchmark_filename, "READ")
benchmark_ws = benchmark_file.Get("w")

results = results_manager("sim", "pt", "eta", import_ws=ws)
benchmark_results = results_manager("indep", "pt", "eta", import_ws=benchmark_ws)

bin_dict = 


graph_refitted = ROOT.TGraph()



c = ROOT.TCanvas("c", "c", 1600, 900)
c.cd()

pad_title = ROOT.TPad("pad_title", "pad_title", 0, 0.9, 1, 1)
pad_plot = ROOT.TPad("pad_plot", "pad_plot", 0, 0, 0.7, 1)
pad_info = ROOT.TPad("pad_info", "pad_info", 0.7, 0, 1, 1)

pad_title.SetMargin(0.1, 0.1, 0.1, 0.1), pad_title.Draw()
pad_plot.SetMargin(0.15, 0.05, 0.15, 0.05), pad_plot.Draw()
pad_info.SetMargin(0.1, 0.1, 0.1, 0.1), pad_info.Draw()

pad_title.cd()
title = ROOT.TLatex()
title.SetTextSize(0.5)
title.DrawLatexNDC(0.5, 0.5, "Simultaneous Fit Correlations")

pad_plot.cd()

