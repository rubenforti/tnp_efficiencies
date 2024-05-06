import ROOT
import os, sys
from array import array
import argparse
import ctypes
from utilities.base_library import bin_dictionary, sumw2_error, binning
from utilities.dataset_utils import gen_import_dictionary, ws_init
from matplotlib import pyplot as plt
from utilities.bkg_utils import bkg_mass_distribution
from copy import deepcopy as dcp

path = os.path.dirname(__file__)
ROOT.gSystem.cd(path)

ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True


colors = { 
    "bkg_WW" : ROOT.kGreen-1,
    "bkg_WZ" : ROOT.kGreen-3,
    "bkg_ZZ" : ROOT.kGreen+1,
    "bkg_TTSemileptonic" : ROOT.kCyan+1,
    "bkg_TTFullyleptonic" : ROOT.kCyan+4,
    "bkg_Ztautau" : ROOT.kMagenta+1,
    "bkg_WplusJets" : ROOT.kOrange+7,
    "bkg_WminusJets" : ROOT.kOrange-3,
    "bkg_SameCharge" : ROOT.kYellow+2,
    "bkg_Zjets" : ROOT.kBlue+4,
    "bkg_total" :  ROOT.kBlack, 
    "mc" : ROOT.kRed,
    "mc_SS" : ROOT.kOrange+4, 
    "pdf_bkg_fit" : ROOT.kRed-2,

    "Diboson" : ROOT.kGreen,
    "Top" : ROOT.kCyan+3,
    "Ztautau" : ROOT.kMagenta+1,
    "Wjets" : ROOT.kOrange+4,
    "Zjets" : ROOT.kBlue+2,
    "Diboson_SS" : ROOT.kGreen-2,
    "Top_SS" : ROOT.kCyan-7,
    "Ztautau_SS" : ROOT.kMagenta-1,
    "Wjets_SS" : ROOT.kOrange-2,
    "Zjets_SS" : ROOT.kBlue-7,
    "SameCharge" : ROOT.kYellow+2,
}

parser = argparse.ArgumentParser(description='Type of background to study')
parser.add_argument("-b", "--background_type", type=str, choices=["Zjets", "Wjets", "all"], default="all",
                    help='the type of background to study')
# parser.add_argument(      "--evaluateSameCharge", action="store_true")
args = parser.parse_args()

print("Studying background: " + args.background_type)


base_dir = os.path.join(path, f"{args.background_type}_tests")
if not os.path.exists(base_dir): os.makedirs(base_dir)

if args.background_type == "Zjets":
    bkg_categories = ["bkg_Zjets"]
elif args.background_type == "Wjets":
    bkg_categories = ["bkg_WplusJets", "bkg_WminusJets"]
else :
    bkg_categories = ["bkg_Zjets", "bkg_WplusJets", "bkg_WminusJets", "bkg_SameCharge"]

import_categories = bkg_categories + ["mc"]

binning_mass = "mass_50_130"
binning_pt = "pt_tracking"
binning_eta = "eta"


import_dictionary = gen_import_dictionary("/scratch/rforti/steve_histograms_2016/tracking", "tracking", import_categories,
                                          ch_set=["plus", "minus"], do_OS_tracking=True,
                                          add_SS_mc=True, add_SS_bkg=True)

bkg_cat_plot = dcp(bkg_categories)
bkg_cat_plot += [cat+"_SS" for cat in bkg_cat_plot if cat!="bkg_SameCharge"]
bkg_cat_plot += ["mc_SS"]

ws = ws_init(import_dictionary, "indep", binning_pt, binning_eta, binning_mass)

ws_filename = f"{base_dir}/{args.background_type}_bkg_OS+SS.root"
ws.writeToFile(ws_filename)

plot_mass_dir = os.path.join(base_dir, "minv_plots_w_sig")
if not os.path.exists(plot_mass_dir): os.makedirs(plot_mass_dir)

bkg_mass_distribution("tracking", ws_filename, bkg_cat_plot, binning_pt, binning_eta,
                      plot_on_signal=True, logscale="hybrid",
                      figpath=base_dir)

os.system(f"rm {plot_mass_dir}/*_pass.pdf")


bins_pt, bins_eta = binning(binning_pt), binning(binning_eta)


h2d_dict = {}

h2d_names = ["OS_SS"] if args.background_type!="all" else ["OS_SameCharge", "SS_SameCharge", "OS+SS_SameCharge"]

for h_name in h2d_names:
    h2d_dict[h_name] =  ROOT.TH2D(f"h2d_ratio_{h_name}", f"{h_name} ratio", 
                                  len(bins_pt)-1, bins_pt, len(bins_eta)-1, bins_eta)

'''

h2d_ratio_OS_SS = ROOT.TH2D("h2d_ratio_OS_SS", "SS/OS ratio", 
                            len(bins_pt)-1, bins_pt, len(bins_eta)-1, bins_eta)

h2d_pull_OS_SS = ROOT.TH2D("h2d_pull_OS_SS", "Pull of SS/OS ratio", 
                            len(bins_pt)-1, bins_pt, len(bins_eta)-1, bins_eta)
lg
'''



bin_dict = bin_dictionary(binning_pt, binning_eta)





for b_key, [_, bin_pt, bin_eta] in bin_dict.items():

    axis = ws.var(f"x_fail_{b_key}")
   
    histo_OS = ROOT.RooDataHist("hist_OS", f"hist_{args.background_type}_OS", ROOT.RooArgList(axis), "x_binning")
    histo_SS = ROOT.RooDataHist("hist_SS", f"hist_{args.background_type}_SS", ROOT.RooArgList(axis), "x_binning")
    for bkg_spec_cat in bkg_categories:
        if bkg_spec_cat == "bkg_SameCharge": continue
        
        histo_OS.add(ws.data(f"Minv_bkg_fail_{b_key}_{bkg_spec_cat.replace('bkg_', '')}"))
        histo_SS.add(ws.data(f"Minv_bkg_fail_{b_key}_{bkg_spec_cat.replace('bkg_', '')}_SS"))

    histo_mc_SS = dcp(ws.data(f"Minv_mc_fail_{b_key}_SS"))
    histo_mc_SS.SetName(f"hist_mcSS_{b_key}")

    if not args.background_type == "all":
        ratios_dict = { 
            "OS_SS" : {
                "histo_num" : histo_OS,
                "histo_den" : histo_SS,
                "array_ratio" : [0]*80,
                "errors_ratio" : [0]*80
            }
        }
    else:
        histo_SameCharge = dcp(ws.data(f"Minv_bkg_fail_{b_key}_SameCharge"))
        histo_SameCharge.SetTitle(f"hist_SameCharge")
        #histo_SameCharge.add(ws.data(f"Minv_bkg_fail_{b_key}_SameCharge"))
        histo_OS_SS = ROOT.RooDataHist("hist_OS_SS", f"hist_OS+SS", ROOT.RooArgList(axis), "x_binning")
        histo_OS_SS.add(histo_OS)
        histo_OS_SS.add(histo_SS)

        ratios_dict = {
            "OS_SameCharge" : {
                "histo_num" : histo_OS,
                "histo_den" : histo_SameCharge,
                "array_ratio" : [0]*80,
                "errors_ratio" : [0]*80
            },
            "SS_SameCharge" : {
                "histo_num" : histo_SS,
                "histo_den" : histo_SameCharge,
                "array_ratio" : [0]*80,
                "errors_ratio" : [0]*80
            },
            "OS+SS_SameCharge" : {
                "histo_num" : histo_OS_SS,
                "histo_den" : histo_SameCharge,
                "array_ratio" : [0]*80,
                "errors_ratio" : [0]*80
            }
        }


    for ratio_key, ratio_dict in ratios_dict.items():
        histo_num = ratio_dict["histo_num"]
        histo_den = ratio_dict["histo_den"]
        array_ratio = ratio_dict["array_ratio"]
        errors_ratio = ratio_dict["errors_ratio"]

        for i in range(80): 
            histo_num.get(i)
            histo_den.get(i)

            try:
                array_ratio[i] = histo_num.weight(i)/histo_den.weight(i)
            except ZeroDivisionError:
                array_ratio[i] = 0
            try:
                errors_ratio[i] = array_ratio[i]*(
                    (histo_den.weightError(ROOT.RooAbsData.SumW2)/histo_den.weight(i))**2 + 
                    (histo_num.weightError(ROOT.RooAbsData.SumW2)/histo_num.weight(i))**2 )**0.5
                if errors_ratio[i] / array_ratio[i] > 1:
                    if array_ratio[i] > 1: array_ratio[i], errors_ratio[i] = 1.5, 1
                    else: array_ratio[i], errors_ratio[i] = 0.5, 1
            except ZeroDivisionError:
                errors_ratio[i] = 0

        ratio_val = histo_num.sumEntries()/histo_den.sumEntries()
        ratio_err = ratio_val*( (sumw2_error(histo_den)/histo_den.sumEntries())**2 +
                                (sumw2_error(histo_num)/histo_num.sumEntries())**2 )**0.5
        
            
        h2d_dict[ratio_key].SetBinContent(bin_pt, bin_eta, ratio_val)
        h2d_dict[ratio_key].SetBinError(bin_pt, bin_eta, ratio_err)

    #ratio_OS_SS.append(histo_OS.sumEntries()/histo_SS.sumEntries())
   
    # h2d_ratio_OS_SS.SetBinContent(bin_pt, bin_eta, ratio_OS_SS[-1])
    # h2d_ratio_OS_SS.SetBinError(bin_pt, bin_eta, error_ratio_OS_SS_Wjets[-1])


        c = ROOT.TCanvas(f"c_{ratio_key}_{b_key}", f"c_{b_key}", 900, 900)
        c.cd()

        pad_info = ROOT.TPad("pad_info", "pad_info", 0, 0.9, 1, 1)
        pad_info.SetMargin(0.25, 0.25, 0.2, 0.1)
        pad_info.Draw()

        pad_plot = ROOT.TPad("pad_plot", "pad_plot", 0, 0.3, 1, 0.9)
        pad_plot.SetMargin(0.1, 0.05, 0.1, 0.01)
        pad_plot.Draw()

        pad_legend = ROOT.TPad("pad_legend", "pad_legend", 0.75, 0.75, 0.95, 0.9)
        pad_legend.SetMargin(0.1, 0.1, 0.1, 0.1)
        pad_legend.Draw()

        pad_ratio = ROOT.TPad("pad_ratio", "pad_ratio", 0, 0, 1, 0.3)
        pad_ratio.SetMargin(0.1, 0.05, 0.1, 0.01)
        pad_ratio.Draw()

        pad_info.cd()
        text = ROOT.TPaveText(0, 0.2, 1, 0.8, "NDC NB")
        text.AddText(f"Bin: {b_key}")
        text.AddText(f"Ratio {ratio_key} = {ratio_val:.2f} #pm {ratio_err:.2f}")
        text.Draw()
        c.Update()  

        pad_plot.cd()
        frame = axis.frame(50, 130, 80)
        frame.SetTitle("")
        frame.GetYaxis().SetTitle("Events/(1 GeV)")
        frame.GetYaxis().SetTitleOffset(1.2)

        if ratio_key not in ["OS+SS_SameCharge", "SS_SameCharge"]:
            histo_num.plotOn(frame, 
                            ROOT.RooFit.LineColor(ROOT.kBlue+2), 
                            ROOT.RooFit.MarkerColor(ROOT.kBlue+2),
                            ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))
        else:
            roopdf_num = ROOT.RooHistPdf(f"roopdf_{histo_num.GetName()}", histo_num.GetTitle(), ROOT.RooArgList(axis), histo_num)
            roopdf_num.plotOn(frame,
                             ROOT.RooFit.LineColor(ROOT.kBlue+2),
                             ROOT.RooFit.Normalization(histo_den.sumEntries()))

        histo_den.plotOn(frame, 
                        ROOT.RooFit.LineColor( ROOT.kBlue-7),
                        ROOT.RooFit.MarkerColor(ROOT.kBlue-7), 
                        ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))
        
        pdf_histo_mcSS = ROOT.RooHistPdf(f"pdf_histo_mcSS_{b_key}", f"pdf_histo_mcSS_{b_key}", 
                                         ROOT.RooArgSet(axis), histo_mc_SS)
        
        pdf_histo_mcSS.plotOn(frame,
                            ROOT.RooFit.LineColor(ROOT.kOrange+4),
                            ROOT.RooFit.Normalization(histo_mc_SS.sumEntries(), ROOT.RooAbsReal.NumEvent))
        
        pad_legend.cd()
        legend = ROOT.TLegend(0.1, 0.1, 0.9, 0.9)
        legend.SetFillColor(0)
        legend.SetTextSize(0.15)
        legend.SetTextAlign(12)
        legend.SetBorderSize(0)
        
        if ratio_key!="OS+SS_SameCharge":
            legend.AddEntry(histo_num, histo_num.GetTitle().replace("hist_", ""), "lep")
        else:
            legend.AddEntry(histo_num, "OS+SS (norm)", "l")
        legend_obj = legend.GetListOfPrimitives().Last()
        legend_obj.SetLineColor(ROOT.kBlue+2)
        legend_obj.SetLineWidth(3)
        legend_obj.SetMarkerColor(ROOT.kBlue+2)
        legend_obj.SetMarkerSize(1)


        legend.AddEntry(histo_den, histo_den.GetTitle().replace("hist_", ""), "lep")
        legend_obj = legend.GetListOfPrimitives().Last()
        legend_obj.SetLineColor(ROOT.kBlue-7)
        legend_obj.SetLineWidth(3)
        legend_obj.SetMarkerColor(ROOT.kBlue-7)
        legend_obj.SetMarkerSize(1)

        legend.AddEntry(pdf_histo_mcSS, "MC_SS", "l")
        legend_obj = legend.GetListOfPrimitives().Last()
        legend_obj.SetLineColor(colors["mc_SS"])
        legend_obj.SetLineWidth(3)
        

        if ratio_key == "OS_SS" and args.background_type=="Zjets":
            pad_plot.cd()
            aux_roohist = dcp(histo_SS)
            aux_roohist.add(histo_mc_SS)

            pdf_aux = ROOT.RooHistPdf(f"pdf_aux_{b_key}", f"pdf_aux_{b_key}", ROOT.RooArgSet(axis), aux_roohist)

            pdf_aux.plotOn(frame,
                        ROOT.RooFit.LineColor(ROOT.kRed), 
                        ROOT.RooFit.LineStyle(1),
                        ROOT.RooFit.LineWidth(1),
                        ROOT.RooFit.Normalization(aux_roohist.sumEntries(), ROOT.RooAbsReal.NumEvent))
            
            pad_legend.cd()
            legend.AddEntry(pdf_aux, "MC_SS + Zjets_SS", "l")
            legend_obj = legend.GetListOfPrimitives().Last()
            legend_obj.SetLineColor(ROOT.kRed)
            legend_obj.SetLineWidth(1)
        
        pad_plot.cd()
        frame.Draw()

        pad_legend.cd()
        legend.Draw()

        if ratio_key=="OS+SS_SameCharge":
            scale_f = histo_den.sumEntries()/histo_num.sumEntries()
            for i in range(80):
                array_ratio[i] = array_ratio[i]*scale_f
                errors_ratio[i] = errors_ratio[i]*scale_f

        pad_ratio.cd()
        ratio_graph = ROOT.TGraphErrors(80,
                                array("d", [50.5+i for i in range(80)]), array("d", array_ratio), 
                                array("d", [0]*80), array("d", errors_ratio))
        ratio_graph.GetXaxis().SetRangeUser(50, 130)
        ratio_graph.SetTitle("")
        ratio_graph.SetMarkerStyle(20)
        ratio_graph.SetMarkerSize(0.5)
        ratio_graph.GetYaxis().SetTitle(ratio_key)
        ratio_graph.GetYaxis().SetTitleSize(0.5)
        ratio_graph.Draw("ZAP")
        hline = ROOT.TLine(50, 1, 130, 1)
        hline.SetLineStyle(2)
        hline.Draw()

        path_ratio = os.path.join(path, base_dir, f"ratios_{ratio_key}")
        if not os.path.exists(path_ratio): os.makedirs(path_ratio)

        c.SaveAs(f"{path_ratio}/{b_key}_{args.background_type}_ratio.pdf")



file_out = ROOT.TFile(f"{base_dir}/{args.background_type}_features.root", "recreate")
file_out.cd()
for h2 in h2d_dict.values(): h2.Write()
file_out.Close()







'''


bin_dict = bin_dictionary(binning_pt, binning_eta)

ws_OS = file_OS.Get("w")
ws_SS = file_SS.Get("w")

sum_os, sum_ss, sum_mc_ss = 0, 0, 0

diff_list = []
ratio_list = []

for b_key in bin_dict.keys():

    histo_OS = ws_OS.data(f"Minv_bkg_fail_{b_key}_Zjets")
    histo_SS = ws_SS.data(f"Minv_bkg_fail_{b_key}_Zjets_SS")
    # histo_mc_SS = ws_SS.data(f"Minv_mc_fail_{b_key}_SS")

    diff = histo_OS.sumEntries() - histo_SS.sumEntries()
    diff_error = (sumw2_error(histo_OS)**2 + sumw2_error(histo_SS)**2)**0.5

    sum_os += histo_OS.sumEntries()
    sum_ss += histo_SS.sumEntries()
    # sum_mc_ss += histo_mc_SS.sumEntries()

    diff_list.append(diff/diff_error)
    ratio_list.append(histo_OS.sumEntries()/histo_SS.sumEntries())

    if diff/diff_error > 15:
        print(b_key)

print(sum_os, sum_ss, sum_mc_ss)

plt.figure()
plt.hist(diff_list, bins=50)
plt.show()
plt.savefig("diff_Zjets.png")

plt.figure()
plt.hist(ratio_list, bins=50)
plt.show()
plt.savefig("ratio_Zjets.png")



    
'''
