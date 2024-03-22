import ROOT
from utilities.base_library import bin_dictionary
from utilities.results_utils import results_manager
from utilities.CMS_lumi import CMS_lumi
from utilities.dataset_utils import import_pdf_library
from utilities.plot_utils import style_settings
from array import array
from copy import copy

ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

gen_folder = "/scratch/rforti/tnp_efficiencies_results/tracking"

check_biasMC_agreement = True

import_pdf_library("RooCMSShape")

file_bmark = ROOT.TFile(gen_folder + "/benchmark/ws_tracking.root", "READ")
ws_bmark = file_bmark.Get("w")

file_bb = ROOT.TFile(gen_folder + "/BBlight_legacySettings/ws_tracking_BBlight.root", "READ")
ws_bb = file_bb.Get("w")

#file_altSig = ROOT.TFile(gen_folder + "/../egm_tnp_results/tracking/current_results.root", "READ")
#hist_rm1_altSig = file_altSig.Get("h_bias_altSig_2d")

res_bmark = results_manager("indep", "pt_tracking", "eta", import_ws=ws_bmark)
res_bb = results_manager("indep", "pt_tracking", "eta", import_ws=ws_bb)


with open(gen_folder+"/../egm_tnp_results/tracking/allEfficiencies.txt", "r") as file_bmark:
    row_list = file_bmark.readlines()

res_nomi = results_manager("indep", "pt_tracking", "eta", import_txt=row_list)
res_altSig = results_manager("indep", "pt_tracking", "eta", import_txt=row_list, altSig_check=True)


if check_biasMC_agreement is True:
    file_pseudo_bmark = ROOT.TFile(gen_folder + "/pseudodata/ws_tracking_pseudodata.root", "READ")
    file_pseudo_bb = ROOT.TFile(gen_folder + "/pseudodata_BBlight_legacySettings/ws_tracking_pseudodata.root", "READ")
    res_pseudo_bmark = results_manager("indep", "pt_tracking", "eta", import_ws=file_pseudo_bmark.Get("w"))
    res_pseudo_bb = results_manager("indep", "pt_tracking", "eta", import_ws=file_pseudo_bb.Get("w"))


c = ROOT.TCanvas("c", "c", 3600, 1800)
c.cd()

pad_title = ROOT.TPad("pad_title", "pad_title", 0.0, 0.96, 1.0, 1.0)
pad_title.SetMargin(0.05, 0.05, 0.05, 0.05), pad_title.Draw()

pad_plot = ROOT.TPad("pad_plot", "pad_plot", 0.0, 0.0, 1.0, 0.96)
pad_plot.SetMargin(0.11, 0.05, 0.1, 0.05), pad_plot.Draw()

# pad_eta = ROOT.TPad("pad_plot", "pad_plot", 0.0, 0.0, 1.0, 0.48)
# pad_eta.SetMargin(0.11, 0.05, 0.1, 0.05), pad_eta.Draw()


pad_title.cd()
titlebox = ROOT.TPaveText(0, 0.1, 1, 0.9, "NDC NB")
titlebox.SetFillColor(0)
# titlebox.SetTextFont(42)
titlebox.SetTextSize(0.7)
plot_title = "Systematic effects compared" if check_biasMC_agreement is False else "Agreement of systematic measured on data w.r.t. the expected bias from MC"
# plot_title = "Systematic effects compared" if check_biasMC_agreement is False else "Comparison between Barlow-Beeston and Benchmark results"

titlebox.AddText(0.5, 0.5, plot_title)
titlebox.Draw()
c.Update()


bin_list_pt, d_bin_list_pt = [], []
bin_list_eta, d_bin_list_eta = [], []

eff_old_list, d_eff_old_list = [], []
eff_list, d_eff_list = [], []


bin_text = []

bin_dict = bin_dictionary("pt_tracking", "eta")

diff_list = [0]*48

idx_plot = 0
for bin_key, [n_bin, idx_pt, idx_eta] in bin_dict.items():
    
    if idx_eta==40: print(bin_key)

    if idx_eta == 1:
        pt_text = bin_key.split("][")[0]
        pt_min, pt_max = pt_text.split("to")
        # print(pt_min, pt_max)
        pt_min = pt_min.replace("[", "")
        pt_min = pt_min.replace(".0", "")
        pt_max = pt_max.replace(".0", "")    
        bin_text.append(f"p_{{T}} #in  [{pt_min},{pt_max}] GeV/c")


    eff_bb, d_eff_bb = res_bb.getEff(bin_key)
    eff_bmark, d_eff_bmark = res_bmark.getEff(bin_key)
    eff_altSig, d_eff_altSig = res_altSig.getEff(bin_key)
    eff_nomi, d_eff_nomi = res_nomi.getEff(bin_key)

    ratio_new = eff_bb/eff_bmark
    ratio_old = eff_altSig/eff_nomi

    eff_list.append(ratio_new - ratio_old)
    d_eff_list.append(ratio_new*((d_eff_bb/eff_bb)**2 + (d_eff_bmark/eff_bmark)**2
                                 #- (2*d_eff_bb*d_eff_bmark/(eff_bb*eff_bmark))
                                )**0.5)

    eff_old_list.append(0)
    d_eff_old_list.append(ratio_old*((d_eff_altSig/eff_altSig)**2 + (d_eff_nomi/eff_nomi)**2)**0.5)

    bin_list_pt.append(n_bin)
    d_bin_list_pt.append(0)

    idx_eta_new = 4*(idx_eta-1) + idx_pt
    bin_list_eta.append(idx_eta_new)
    d_bin_list_eta.append(0)
    
    '''
    if idx_pt==1 or idx_pt==4:
        stat_diff = ratio_old-ratio_new if idx_pt==1 else ratio_new-ratio_old
    else:
        stat_diff = -abs(ratio_new-ratio_old)
    diff_list[idx_eta-1] += stat_diff
    '''
    
    if check_biasMC_agreement is True:
        effMC_bmark, d_effMC_bmark = res_pseudo_bmark.getEff(bin_key)
        effMC_bb, d_effMC_bb = res_pseudo_bb.getEff(bin_key)
        eff_list[-1] = eff_bb*(effMC_bmark/effMC_bb) - eff_bmark
        # eff_list[-1] = eff_bb - eff_bmark
        d_eff_list[-1] = (eff_bb*effMC_bmark)*(((d_eff_bb/eff_bb)**2 + (d_effMC_bmark/effMC_bmark)**2)**0.5)

        d_eff_old_list[-1] = d_eff_bmark

print("Mean : ", sum(eff_list)/len(eff_list))
print("Std dev : ", (sum([(x - sum(eff_list)/len(eff_list))**2 for x in eff_list])/len(eff_list))**0.5)

    
'''
eta_chosen_idx = -1
max_diff = 0
for i in range(len(diff_list)):
    if diff_list[i] > max_diff:
        max_diff = diff_list[i]
        eta_chosen = i+1
        
print(eta_chosen, max_diff)
'''

pad_plot.cd()

new_eff_graph = ROOT.TGraphErrors(len(bin_list_pt), 
                                  array("d", bin_list_pt), array("d", eff_list), 
                                  array("d", d_bin_list_pt), array("d", d_eff_list))

old_eff_graph = ROOT.TGraphErrors(len(bin_list_pt),
                                  array("d", bin_list_pt), array("d", eff_old_list),
                                  array("d", d_bin_list_pt), array("d", d_eff_old_list))


lowlim, uplim = -0.028, 0.04

new_eff_graph.SetMarkerStyle(20)
new_eff_graph.SetMarkerSize(1.7)
new_eff_graph.SetMarkerColor(ROOT.kBlack)
new_eff_graph.SetLineWidth(2)
new_eff_graph.SetTitle("")

new_eff_graph.GetXaxis().SetTitle("Bin number (p_{T})")
new_eff_graph.GetXaxis().SetRangeUser(-1, len(bin_list_pt)+2)
new_eff_graph.GetXaxis().SetTitleSize(0.04)
new_eff_graph.GetXaxis().SetTitleOffset(1.1)
new_eff_graph.GetXaxis().SetLabelSize(0.025)
new_eff_graph.GetXaxis().SetTickSize(0.0)

# new_eff_graph.GetYaxis().SetTitle("#varepsilon_{B.B.}/#varepsilon_{Bmark} - #varepsilon_{AltSig}/#varepsilon_{Nomi}")

new_eff_graph.GetYaxis().SetTitle("Syst. effect")
new_eff_graph.GetYaxis().SetTitleSize(0.04)
new_eff_graph.GetYaxis().SetRangeUser(lowlim, uplim)
new_eff_graph.GetYaxis().SetTickSize(0.0)
new_eff_graph.GetYaxis().SetLabelSize(0.025)
new_eff_graph.GetYaxis().SetTitleOffset(1)

new_eff_graph.Draw("ZAP")

old_eff_graph.SetTitle("")
old_eff_graph.SetLineWidth(1)
# old_eff_graph.SetFillColor(4)
old_eff_graph.SetFillColor(ROOT.kOrange-2)
# old_eff_graph.SetFillStyle(3015)
#old_eff_graph.SetFillStyle(3001)

# old_eff_graph.GetYaxis().SetRangeUser(0.885, 1.045)
# old_eff_graph.GetXaxis().SetRangeUser(0, 721)
old_eff_graph.Draw("3 same Z")

new_eff_graph.Draw("same ZP")



c.Update()


nptBins = 4
etarange = 48

offsetXaxis = 0.5

vertline = ROOT.TLine(36, 0, 36, 1.2)
vertline.SetLineColor(11)
vertline.SetLineStyle(4)
vertline.SetLineWidth(4)

for i in range(1, nptBins): # do not need line at canvas borders 
    vertline.DrawLine(48*i-offsetXaxis, lowlim, 48*i-offsetXaxis, uplim)


bintext = ROOT.TLatex()
#bintext.SetNDC()
bintext.SetTextSize(0.025)  # 0.03
bintext.SetTextFont(42)
bintext.SetTextAngle(45)



sliceLabelOffset = 6.0
    
for i in range(0, nptBins):
    ytext = lowlim + 0.003
    bintext.DrawLatex(etarange*i + etarange/sliceLabelOffset, ytext, bin_text[i])


new_eff_legend = copy(new_eff_graph)
new_eff_legend.SetMarkerColor(ROOT.kBlack)
new_eff_legend.SetLineColor(ROOT.kBlack)
new_eff_legend.SetMarkerSize(2)
new_eff_legend.SetLineWidth(2)

old_eff_legend = copy(old_eff_graph)
# old_eff_legend.SetMarkerColor(ROOT.kRed)
#old_eff_legend.SetLineColor(ROOT.kRed)
old_eff_legend.SetMarkerSize(0)
old_eff_legend.SetLineWidth(0)


legend = ROOT.TLegend(0.12, 0.78, 0.42, 0.9)
# legend.SetBorderSize(0)
legend.SetFillStyle(1001)
legend.SetFillColor(0)
legend.SetTextFont(42)
legend.SetTextSize(0.03)
if not check_biasMC_agreement:
    legend.AddEntry(new_eff_legend, "#varepsilon_{B.B.}/#varepsilon_{Benchmark} - #varepsilon_{AltSig}/#varepsilon_{Nominal}", "lp")
    legend.AddEntry(old_eff_legend, "Old syst. effect", "F")
else:
    legend.AddEntry(new_eff_legend, "#varepsilon_{B.B.}(#varepsilon_{MC, Benchmark} / #varepsilon_{MC, B.B}) - #varepsilon_{Benchmark}", "lp")
    # legend.AddEntry(new_eff_legend, "#varepsilon_{B.B.} - #varepsilon_{Benchmark}", "lp")
    legend.AddEntry(old_eff_legend, "#sigma_{#varepsilon_{Benchmark}}", "F")

legend.Draw()

CMS_lumi(pad_plot, 5, 0, simulation=False)
pad_plot.Update()


c.SaveAs("check_mc_agreement_unrolled.png")

'''

c_hist = ROOT.TCanvas("c_hist", "c_hist", 1800, 1800)
c_hist.cd()

style_settings()

pad_title.Draw()
c_hist.Update()

pad_hist = ROOT.TPad("pad_hist", "pad_hist", 0.0, 0.0, 1.0, 0.96)
pad_hist.SetMargin(0.11, 0.05, 0.1, 0.05), pad_hist.Draw()

lowlim_hist, uplim_hist = -0.02, 0.022

hist_bin_1 = ROOT.TH1D("hist_bin_1", "hist_bin_1", 25, lowlim_hist, uplim_hist)
hist_bin_2 = ROOT.TH1D("hist_bin_2", "hist_bin_2", 25, lowlim_hist, uplim_hist)
hist_bin_3 = ROOT.TH1D("hist_bin_3", "hist_bin_3", 25, lowlim_hist, uplim_hist)
hist_bin_4 = ROOT.TH1D("hist_bin_4", "hist_bin_4", 25, lowlim_hist, uplim_hist)
hist_total = ROOT.TH1D("hist_total", "hist_total", 25, lowlim_hist, uplim_hist)

for i in range(len(eff_list)):
    hist_total.Fill(eff_list[i])
    if bin_list_eta[i] < 48:
        hist_bin_1.Fill(eff_list[i])
    elif bin_list_eta[i] < 96:
        hist_bin_2.Fill(eff_list[i])
    elif bin_list_eta[i] < 144:
        hist_bin_3.Fill(eff_list[i])
    else:
        hist_bin_4.Fill(eff_list[i])

hist_total.SetLineColor(ROOT.kBlack)
hist_bin_1.SetLineColor(ROOT.kOrange)
hist_bin_2.SetLineColor(ROOT.kRed)
hist_bin_3.SetLineColor(ROOT.kBlue)
hist_bin_4.SetLineColor(ROOT.kGreen)

hist_total.SetLineWidth(0)
hist_bin_1.SetLineWidth(2)
hist_bin_2.SetLineWidth(2)
hist_bin_3.SetLineWidth(2)
hist_bin_4.SetLineWidth(2)

hist_total.SetContour(51)
hist_total.Draw("")
hist_total.GetZaxis().SetTitle("")
hist_total.SetTitleSize(0)
hist_total.SetTitle("")

hist_total.GetXaxis().SetTitle("Syst. effect")
hist_total.GetXaxis().SetTitleOffset(1.2)
hist_total.GetXaxis().SetTitleSize(0.04)
hist_total.GetXaxis().SetLabelSize(0.02)


hist_bin_1.Draw("same")
hist_bin_2.Draw("same")
hist_bin_3.Draw("same")
hist_bin_4.Draw("same")

legend_hist = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legend_hist.AddEntry(hist_total, "All bins", "l")
legend_hist.AddEntry(hist_bin_1, "pT #in 24-35", "l")
legend_hist.AddEntry(hist_bin_2, "pT #in 35-45", "l")
legend_hist.AddEntry(hist_bin_3, "pT #in 45-55", "l")
legend_hist.AddEntry(hist_bin_4, "pT #in 55-65", "l")
legend_hist.Draw()





c_hist.Update()
c_hist.SaveAs("hist_eff.png")


'''