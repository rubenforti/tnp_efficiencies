import ROOT
from utilities.base_library import bin_dictionary, sumw2_error
from utilities.results_utils import results_manager, efficiency_from_res, eval_efficiency
from utilities.CMS_lumi import CMS_lumi
from array import array
from copy import copy

ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

gen_folder = "/scratch/rforti/tnp_efficiencies_results/tracking"

file_bmark = ROOT.TFile(gen_folder + "/benchmark/ws_tracking.root", "READ")
ws_bmark = file_bmark.Get("w")

file_bb = ROOT.TFile(gen_folder + "/BBlight_legacySettings/ws_tracking_BBlight.root", "READ")
ws_bb = file_bb.Get("w")

file_altSig = ROOT.TFile(gen_folder + "/current_results.root", "READ")
hist_rm1_altSig = file_altSig.Get("h_bias_altSig_2d")

bin_dict = bin_dictionary()

res_obj = results_manager("indep", "pt_tracking", "eta", import_ws=ws_bb)


c = ROOT.TCanvas("c", "c", 2400, 1200)
c.cd()

pad_title = ROOT.TPad("pad_title", "pad_title", 0.0, 0.95, 1.0, 1.0)
pad_title.SetMargin(0.05, 0.05, 0.05, 0.05), pad_title.Draw()

pad_main = ROOT.TPad("pad_main", "pad_main", 0.0, 0.35, 1.0, 0.95)
pad_main.SetMargin(0.06, 0.025, 0.0, 0.05), pad_main.Draw()

pad_sf = ROOT.TPad("pad_sf", "pad_sf", 0.0, 0.0, 1.0, 0.35)
pad_sf.SetMargin(0.06, 0.025, 0.25, 0.0), pad_sf.Draw()


pad_title.cd()
titlebox = ROOT.TPaveText(0, 0.1, 1, 0.9, "NDC NB")
titlebox.SetFillColor(0)
# titlebox.SetTextFont(42)
titlebox.SetTextSize(0.7)
titlebox.AddText(0.5, 0.5, f"Unrolled efficiency")
titlebox.Draw()
c.Update()


eff_list, d_eff_list, bin_list, d_bin_list = [], [], [], []

eff_mc_list, d_eff_mc_list = [], []
sf_list, d_sf_list = [], []


bin_text = []


idx_plot = 0
for bin_key, [n_bin, idx_pt, idx_eta] in bin_dict.items():


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
    d_bin_list.append(1)

    # data_eff_graph.SetPoint(n_bin, n_bin, eff)
    # data_eff_graph.SetPointError(n_bin, d_eff, d_eff)

    mc_pass, mc_fail = ws_bmark.data(f"Minv_mc_pass_{bin_key}"), ws_bmark.data(f"Minv_mc_fail_{bin_key}")

    eff_mc, d_eff_mc = eval_efficiency(mc_pass.sumEntries(), mc_fail.sumEntries(), 
                                       sumw2_error(mc_pass), sumw2_error(mc_fail))

    eff_mc_list.append(eff_mc)
    d_eff_mc_list.append(d_eff_mc)

    sf_list.append(eff/eff_mc)
    d_sf_list.append(((d_eff/eff)**2 + (d_eff_mc/eff_mc)**2)**0.5)



pad_main.cd()
data_eff_graph = ROOT.TGraphErrors(720, 
                                   array("d", bin_list), array("d", eff_list), 
                                   array("d", d_bin_list), array("d", d_eff_list))

mc_eff_graph = ROOT.TGraphErrors(720,
                                 array("d", bin_list), array("d", eff_mc_list), 
                                 array("d", d_bin_list), array("d", d_eff_mc_list))


data_eff_graph.SetMarkerStyle(20)
data_eff_graph.SetMarkerSize(0.7)
data_eff_graph.SetMarkerColor(ROOT.kBlack)
data_eff_graph.SetLineWidth(1)
data_eff_graph.SetTitle("")
data_eff_graph.GetXaxis().SetTitle("")
data_eff_graph.GetXaxis().SetLabelSize(0)
data_eff_graph.GetYaxis().SetTitle("Efficiency")
data_eff_graph.GetYaxis().SetTickSize(0.00)
data_eff_graph.GetYaxis().SetLabelSize(0.04)
data_eff_graph.GetYaxis().SetTitleOffset(0.4)
data_eff_graph.GetYaxis().SetTitleSize(0.06)

data_eff_graph.GetYaxis().SetRangeUser(0.885, 1.045)
data_eff_graph.GetXaxis().SetRangeUser(0, 721)

data_eff_graph.Draw("ZAP")

mc_eff_graph.SetLineWidth(1)
mc_eff_graph.SetMarkerStyle(20)
mc_eff_graph.SetMarkerSize(0.7)
mc_eff_graph.SetMarkerColor(ROOT.kRed)
mc_eff_graph.SetTitle("")
mc_eff_graph.SetLineColor(ROOT.kRed)

mc_eff_graph.GetYaxis().SetRangeUser(0.885, 1.045)
mc_eff_graph.GetXaxis().SetRangeUser(0, 721)

mc_eff_graph.Draw("PZ same")
c.Update()



offsetXaxis = 0.5

vertline = ROOT.TLine(36, 0, 36, 1.2)
vertline.SetLineColor(11)
vertline.SetLineStyle(4)

for i in range(1, 15): # do not need line at canvas borders 
        vertline.DrawLine(48*i-offsetXaxis, 0.885, 48*i-offsetXaxis, 1.045)

bintext = ROOT.TLatex()
#bintext.SetNDC()
bintext.SetTextSize(0.035)  # 0.03
bintext.SetTextFont(42)
bintext.SetTextAngle(45)

nptBins = 15
etarange = 48

sliceLabelOffset = 12.0
    
for i in range(0, 15):
    ytext = 1.0
    bintext.DrawLatex(etarange*i + etarange/sliceLabelOffset, ytext, bin_text[i])


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


legend = ROOT.TLegend(0.8, 0.25, 0.95, 0.5)
# legend.SetBorderSize(0)
legend.SetFillStyle(1001)
legend.SetFillColor(0)
legend.SetTextFont(42)
legend.SetTextSize(0.065)
legend.AddEntry(data_eff_legend, "#varepsilon Data", "lp")
legend.AddEntry(mc_eff_legend, "#varepsilon MC", "lp")
legend.Draw()


CMS_lumi(pad_main, 5, 0, simulation=False)
pad_main.Update()



pad_sf.cd()
sf_graph = ROOT.TGraphErrors(720,
                            array("d", bin_list), array("d", sf_list), 
                            array("d", d_bin_list), array("d", d_sf_list))
sf_graph.SetMarkerStyle(20)
sf_graph.SetMarkerSize(0.7)
sf_graph.SetMarkerColor(ROOT.kBlue)
sf_graph.SetLineColor(ROOT.kBlue)
sf_graph.SetLineWidth(1)
sf_graph.SetTitle("")
sf_graph.GetYaxis().SetRangeUser(0.985, 1.035)
sf_graph.GetXaxis().SetRangeUser(0, 721)
sf_graph.GetXaxis().SetTitle("Bin number")
sf_graph.GetXaxis().SetLabelSize(0.09)
sf_graph.GetXaxis().SetTitleSize(0.11)
sf_graph.GetXaxis().SetTickSize(0.05)
sf_graph.GetXaxis().SetTitleOffset(0.9)
sf_graph.GetYaxis().SetTitle("Data / MC")
sf_graph.GetYaxis().SetTitleSize(0.08)
sf_graph.GetYaxis().SetTickSize(0.00)
sf_graph.GetYaxis().SetLabelSize(0.07)
sf_graph.GetYaxis().SetTitleOffset(0.3)
sf_graph.GetYaxis().SetNdivisions(505)

sf_graph.Draw("ZAP")


hline_gr = ROOT.TGraph(2, array("d", [0, 721]), array("d", [1, 1]))
hline_gr.SetLineColor(ROOT.kBlack)
hline_gr.SetLineStyle(2)
hline_gr.Draw("L same")

for i in range(1, 15): # do not need line at canvas borders 
        vertline.DrawLine(48*i-offsetXaxis, 0.985, 48*i-offsetXaxis, 1.035)


c.Update()








c.SaveAs("unrolled_eff_bmark.pdf")

