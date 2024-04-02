import ROOT
from utilities.base_library import binning, bin_dictionary, sumw2_error
from utilities.results_utils import results_manager, efficiency_from_res, eval_efficiency
from utilities.CMS_lumi import CMS_lumi
from utilities.dataset_utils import import_pdf_library
from array import array
from copy import copy

ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

import_pdf_library("RooCMSShape")


gen_folder = "/scratch/rforti/tnp_efficiencies_results/tracking"

file_bmark = ROOT.TFile(gen_folder + "/benchmark/ws_tracking.root", "READ")
ws_bmark = file_bmark.Get("w")

file_bb = ROOT.TFile(gen_folder + "/BBlight_legacySettings/ws_tracking_BBlight.root", "READ")
ws_bb = file_bb.Get("w")

binning_pt_name = "pt_tracking"
binning_eta_name = "eta"


nptBins = len(binning(binning_pt_name))-1
netaBins = len(binning(binning_eta_name))-1
ntotBins = nptBins*netaBins


lowlim, uplim = 0.972, 1.007
lowlim_sf, uplim_sf = 0.985, 1.015


res_obj = results_manager("indep", "pt_tracking", "eta", import_ws=ws_bb)


c = ROOT.TCanvas("c", "c", 3600, 2400)
c.cd()

pad_title = ROOT.TPad("pad_title", "pad_title", 0.0, 0.95, 1.0, 1.0)
pad_title.SetMargin(0.05, 0.05, 0.05, 0.05), pad_title.Draw()

pad_main = ROOT.TPad("pad_main", "pad_main", 0.0, 0.35, 1.0, 0.95)
pad_main.SetMargin(0.08, 0.025, 0.0, 0.05), pad_main.Draw()

pad_sf = ROOT.TPad("pad_sf", "pad_sf", 0.0, 0.0, 1.0, 0.35)
pad_sf.SetMargin(0.08, 0.025, 0.25, 0.0), pad_sf.Draw()


pad_title.cd()
titlebox = ROOT.TPaveText(0, 0.1, 1, 0.9, "NDC NB")
titlebox.SetFillColor(0)
titlebox.SetTextSize(0.7)
titlebox.AddText(0.5, 0.5, f"Unrolled efficiency")
titlebox.Draw()
c.Update()


eff_list, d_eff_list, bin_list, d_bin_list = [], [], [], []

eff_mc_list, d_eff_mc_list = [], []
sf_list, d_sf_list = [], []


bin_text = []


idx_plot = 0
for bin_key, [n_bin, idx_pt, idx_eta] in bin_dictionary(binning_pt_name, binning_eta_name).items():


    if idx_eta == 1:
        pt_text = bin_key.split("][")[0]
        pt_min, pt_max = pt_text.split("to")
        # print(pt_min, pt_max)
        pt_min = pt_min.replace("[", "")
        pt_min = pt_min.replace(".0", "")
        pt_max = pt_max.replace(".0", "")    
        bin_text.append(f"p_{{T}} #in  [{pt_min},{pt_max}] GeV/c")

    eff, d_eff = res_obj.getEff(bin_key)

    eff_list.append(eff)
    d_eff_list.append(d_eff)
    bin_list.append(n_bin)
    d_bin_list.append(0)

    # data_eff_graph.SetPoint(n_bin, n_bin, eff)
    # data_eff_graph.SetPointError(n_bin, d_eff, d_eff)

    mc_pass, mc_fail = ws_bb.data(f"Minv_mc_pass_{bin_key}"), ws_bb.data(f"Minv_mc_fail_{bin_key}")

    eff_mc, d_eff_mc = eval_efficiency(mc_pass.sumEntries(), mc_fail.sumEntries(), 
                                       sumw2_error(mc_pass), sumw2_error(mc_fail))

    eff_mc_list.append(eff_mc)
    d_eff_mc_list.append(d_eff_mc)

    sf_list.append(eff/eff_mc)
    d_sf_list.append(((d_eff/eff)**2 + (d_eff_mc/eff_mc)**2)**0.5)



pad_main.cd()
data_eff_graph = ROOT.TGraphErrors(ntotBins, 
                                   array("d", bin_list), array("d", eff_list), 
                                   array("d", d_bin_list), array("d", d_eff_list))

mc_eff_graph = ROOT.TGraphErrors(ntotBins,
                                 array("d", bin_list), array("d", eff_mc_list), 
                                 array("d", d_bin_list), array("d", d_eff_mc_list))


data_eff_graph.SetMarkerStyle(20)
data_eff_graph.SetMarkerSize(1.7)
data_eff_graph.SetMarkerColor(ROOT.kBlack)
data_eff_graph.SetLineWidth(2)
data_eff_graph.SetTitle("")
data_eff_graph.GetXaxis().SetTitle("")
data_eff_graph.GetXaxis().SetLabelSize(0)
data_eff_graph.GetYaxis().SetTitle("Efficiency")
data_eff_graph.GetYaxis().SetTickSize(0.00)
data_eff_graph.GetYaxis().SetLabelSize(0.04)
data_eff_graph.GetYaxis().SetTitleOffset(0.6)
data_eff_graph.GetYaxis().SetTitleSize(0.06)

data_eff_graph.GetYaxis().SetRangeUser(lowlim, uplim)
data_eff_graph.GetXaxis().SetRangeUser(-1, ntotBins+2)

data_eff_graph.Draw("ZAP")

mc_eff_graph.SetMarkerStyle(20)
mc_eff_graph.SetMarkerSize(1.7)
mc_eff_graph.SetMarkerColor(ROOT.kRed)
mc_eff_graph.SetLineWidth(2)
mc_eff_graph.SetTitle("")
mc_eff_graph.SetLineColor(ROOT.kRed)

mc_eff_graph.GetYaxis().SetRangeUser(lowlim, uplim)
mc_eff_graph.GetXaxis().SetRangeUser(-1, ntotBins+2)

mc_eff_graph.Draw("PZ same")
c.Update()

offsetXaxis = 0.5

vertline = ROOT.TLine(36, lowlim, 36, uplim)
vertline.SetLineColor(11)
vertline.SetLineStyle(4)
vertline.SetLineWidth(6)

for i in range(1, nptBins): # do not need line at canvas borders 
        vertline.DrawLine(netaBins*i-offsetXaxis, lowlim, netaBins*i-offsetXaxis, uplim)

bintext = ROOT.TLatex()
#bintext.SetNDC()
bintext.SetTextSize(0.035)  # 0.03
bintext.SetTextFont(42)
bintext.SetTextAngle(45)


sliceLabelOffset = 12.0
    
for i in range(0, nptBins):
    ytext = lowlim + 0.003
    bintext.DrawLatex(netaBins*i + netaBins/sliceLabelOffset, ytext, bin_text[i])


data_eff_legend = copy(data_eff_graph)
data_eff_legend.SetMarkerColor(ROOT.kBlack)
data_eff_legend.SetLineColor(ROOT.kBlack)
data_eff_legend.SetMarkerSize(2)
data_eff_legend.SetLineWidth(2)

mc_eff_legend = copy(mc_eff_graph)
mc_eff_legend.SetMarkerColor(ROOT.kRed)
mc_eff_legend.SetLineColor(ROOT.kRed)
mc_eff_legend.SetMarkerSize(2)
mc_eff_legend.SetLineWidth(2)

legend = ROOT.TLegend(0.09, 0.78, 0.2, 0.92)
# legend.SetBorderSize(0)
legend.SetFillStyle(1001)
legend.SetFillColor(0)
legend.SetTextFont(42)
legend.SetTextSize(0.04)
legend.AddEntry(data_eff_legend, "#varepsilon Data", "lp")
legend.AddEntry(mc_eff_legend, "#varepsilon MC", "lp")
legend.Draw()


CMS_lumi(pad_main, 5, 0, simulation=False)
pad_main.Update()



pad_sf.cd()
sf_graph = ROOT.TGraphErrors(ntotBins,
                            array("d", bin_list), array("d", sf_list), 
                            array("d", d_bin_list), array("d", d_sf_list))
sf_graph.SetMarkerStyle(20)
sf_graph.SetMarkerSize(1.7)
sf_graph.SetMarkerColor(ROOT.kBlue)
sf_graph.SetLineColor(ROOT.kBlue)
sf_graph.SetLineWidth(2)
sf_graph.SetTitle("")
sf_graph.GetYaxis().SetRangeUser(lowlim_sf, uplim_sf)
sf_graph.GetXaxis().SetRangeUser(-1, ntotBins+2)
sf_graph.GetXaxis().SetTitle("Bin number")
sf_graph.GetXaxis().SetLabelSize(0.09)
sf_graph.GetXaxis().SetTitleSize(0.11)
sf_graph.GetXaxis().SetTickSize(0.05)
sf_graph.GetXaxis().SetTitleOffset(0.9)
sf_graph.GetYaxis().SetTitle("Data / MC")
sf_graph.GetYaxis().SetTitleSize(0.06)
sf_graph.GetYaxis().SetTickSize(0.00)
sf_graph.GetYaxis().SetLabelSize(0.07)
sf_graph.GetYaxis().SetTitleOffset(0.6)
sf_graph.GetYaxis().SetNdivisions(505)

sf_graph.Draw("ZAP")


hline = ROOT.TLine(-1, 1, ntotBins+2, 1)
hline.SetLineColor(12)
hline.SetLineStyle(10)
hline.SetLineWidth(2)
hline.DrawLine(-1, 1, ntotBins+2, 1)

for i in range(1, nptBins): # do not need line at canvas borders 
        vertline.DrawLine(netaBins*i-offsetXaxis, lowlim_sf, netaBins*i-offsetXaxis, uplim_sf)


c.Update()








c.SaveAs("unrolled_eff_BBlight_dR03.png")

