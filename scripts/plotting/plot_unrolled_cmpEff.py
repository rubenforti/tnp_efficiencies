import ROOT
import sys
import argparse
from utilities.base_lib import binnings, import_pdf_library
from utilities.results_manager import results_manager
from utilities.CMS_lumi import CMS_lumi
from utilities.binning_utils import bin_dictionary
from utilities.plot_utils import style_settings
from array import array
from copy import copy


ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

gen_folder = "/scratch/rforti/tnp_efficiencies_results/reco"

binning_pt_name = "pt_reco"
binning_eta_name = "eta"

import_pdf_library("RooCMSShape")

parser = argparse.ArgumentParser(description="Plot the unrolled comparison of efficiencies or scale factors")

parser.add_argument("-sf", "--scale_factors", action="store_true", help="Analyze scale factors")
parser.add_argument("-cmp_syst", "--compare_systematic_effects", action="store_true", help="Compare systematic effects")
parser.add_argument("-biasMC", "--check_biasMC_agreement", action="store_true", help="Check agreement of syst effect measured on data w.r.t. the expected bias from MC")
parser.add_argument("-yLims", "--yLims", nargs=2, type=float, default=[-0.028, 0.04], help="Set the y-axis limits")
parser.add_argument(           "--output_subscript", default="", help="Subscript to be added to the output plot name") 


args = parser.parse_args()

analyze_scale_factors = args.scale_factors
if analyze_scale_factors:
    qnt = "Scale factors"
else:
    qnt = "Efficiencies"

cmp_syst_effects = args.compare_systematic_effects
check_biasMC_agreement = args.check_biasMC_agreement

if check_biasMC_agreement and(analyze_scale_factors or cmp_syst_effects): 
    sys.exit("Cannot check biasMC agreement and analyze scale factors or compare systematic effects at the same time. Exiting...")


file_nomi = ROOT.TFile(gen_folder + "/legacy_fit/ws_reco.root", "READ")
ws_nomi = file_nomi.Get("w")
res_nomi = results_manager("indep", binning_pt_name, binning_eta_name, import_ws=ws_nomi)

'''
with open(gen_folder+"/../egm_tnp_results/trackingplus/allEfficiencies.txt", "r") as file_nomi:
    row_list = file_nomi.readlines()

res_nomi = results_manager("indep", "pt_tracking", "eta", import_txt=row_list)
'''

file_alt = ROOT.TFile(gen_folder + "/BBlight_legacySettings/ws_reco_BBlight.root", "READ")
ws_alt = file_alt.Get("w")
res_alt = results_manager("indep", binning_pt_name, binning_eta_name, import_ws=ws_alt)

plot_name = qnt + " comparison"


nptBins = len(binnings[binning_pt_name])-1
netaBins = len(binnings[binning_eta_name])-1
ntotBins = nptBins*netaBins


if cmp_syst_effects is True:
    plot_name = "Syst effects comparison on " + qnt.lower()
    '''
    res_nomi_old = res_nomi
    file_alt_old = ROOT.TFile(gen_folder + "/BBlight_legacySettings_dRsig01/ws_tracking_BBlight.root", "READ")
    ws_alt_old = file_alt_old.Get("w")
    res_alt_old = results_manager("indep", "pt_tracking", "eta", import_ws=ws_alt_old)
    '''
    with open(gen_folder+"/../egm_tnp_results/tracking/allEfficiencies.txt", "r") as file_nomi:
        row_list = file_nomi.readlines()

    res_nomi_old = results_manager("indep", binning_pt_name, binning_eta_name, import_txt=row_list)
    res_alt_old = results_manager("indep",  binning_pt_name, binning_eta_name, import_txt=row_list, altSig_check=True)
    

if check_biasMC_agreement is True:
    plot_name = "Agreement of syst effect measured on data wrt expected bias from MC"
    file_pseudo_nomi = ROOT.TFile(gen_folder + "/pseudodata/ws_tracking_pseudodata.root", "READ")
    file_pseudo_alt = ROOT.TFile(gen_folder + "/pseudodata_BBlight_legacySettings/ws_tracking_pseudodata.root", "READ")
    res_pseudo_nomi = results_manager("indep", binning_pt_name, binning_eta_name, import_ws=file_pseudo_nomi.Get("w"))
    res_pseudo_alt = results_manager("indep", binning_pt_name, binning_eta_name, import_ws=file_pseudo_alt.Get("w"))


c = ROOT.TCanvas("c", "c", 3600, 1800)
c.cd()

#style_settings()

pad_title = ROOT.TPad("pad_title", "pad_title", 0.0, 0.95, 1.0, 1.0)
pad_title.SetMargin(0.05, 0.05, 0.05, 0.05), pad_title.Draw()

pad_plot = ROOT.TPad("pad_plot", "pad_plot", 0.0, 0.0, 1.0, 0.95)
pad_plot.SetMargin(0.08, 0.035, 0.1, 0.05), pad_plot.Draw()

# pad_eta = ROOT.TPad("pad_plot", "pad_plot", 0.0, 0.0, 1.0, 0.netaBins)
# pad_eta.SetMargin(0.11, 0.05, 0.1, 0.05), pad_eta.Draw()


pad_title.cd()
titlebox = ROOT.TPaveText(0, 0.1, 1, 0.9, "NDC NB")
titlebox.SetFillColor(0)
titlebox.SetTextSize(0.7)
titlebox.AddText(0.5, 0.5, plot_name + f" ({args.output_subscript})" if args.output_subscript else plot_name)
titlebox.Draw()
c.Update()


bin_list_pt, d_bin_list_pt = [], []
bin_list_eta, d_bin_list_eta = [], []

plot_list, d_plot_list, err_band_list = [], [], []


bin_text = []

diff_list = [0]*netaBins

x_values = [0]*netaBins

idx_plot = 0
for bin_key, [n_bin, idx_pt, idx_eta] in bin_dictionary(binning_pt_name, binning_eta_name).items():
    
    #if idx_eta==40: print(bin_key)

    if idx_eta == 1:
        pt_text = bin_key.split("][")[0]
        pt_min, pt_max = pt_text.split("to")
        # print(pt_min, pt_max)
        pt_min = pt_min.replace("[", "")
        pt_min = pt_min.replace(".0", "")
        pt_max = pt_max.replace(".0", "")    
        bin_text.append(f"p_{{T}} #in  [{pt_min},{pt_max}] GeV/c")


    eff_alt, d_eff_alt = res_alt.getEff(bin_key)
    eff_nomi, d_eff_nomi = res_nomi.getEff(bin_key)

    if n_bin == 58:
        print(bin_key, eff_alt, d_eff_alt, eff_nomi, d_eff_nomi)

    band_amplitude = d_eff_nomi

    if not check_biasMC_agreement:

        store_qnt, d_store_qnt = eff_alt-eff_nomi, d_eff_alt
        if cmp_syst_effects is True:
            eff_alt_old, d_eff_alt_old = res_alt_old.getEff(bin_key)
            eff_nomi_old, d_eff_nomi_old = res_nomi_old.getEff(bin_key)
            store_qnt = (eff_alt-eff_nomi)/eff_nomi
            d_store_qnt = (eff_alt/(eff_nomi**2))*((eff_nomi*d_eff_alt)**2 + (eff_alt*d_eff_nomi)**2)**0.5
            # x_values[idx_eta-1] = (eff_alt_old-eff_nomi_old)/eff_nomi_old
            band_amplitude = (eff_alt_old-eff_nomi_old)/eff_nomi_old
                #(eff_alt_old/(eff_nomi_old**2))*((eff_nomi_old*d_eff_alt_old)**2 + (eff_alt_old*d_eff_nomi_old)**2)**0.5

        if analyze_scale_factors:
            effMC_alt, d_effMC_alt = res_alt.getEffMC(bin_key)
            effMC_nomi, d_effMC_nomi = res_nomi.getEffMC(bin_key)
            sf_alt, sf_nomi = eff_alt/effMC_alt, eff_nomi/effMC_nomi
            if cmp_syst_effects is False:
                store_qnt, d_store_qnt = sf_alt - sf_nomi, d_eff_alt/effMC_alt
                band_amplitude = d_eff_nomi/effMC_nomi
            else:
                effMC_alt_old, d_effMC_alt_old = res_alt_old.getEffMC(bin_key)
                effMC_nomi_old, d_effMC_nomi_old = res_nomi_old.getEffMC(bin_key)
                sf_alt_old, sf_nomi_old = eff_alt_old/effMC_alt_old, eff_nomi_old/effMC_nomi_old
                store_qnt = ((sf_alt-sf_nomi)/sf_nomi)
                x_values[idx_eta-1] = (sf_alt_old-sf_nomi_old)/sf_nomi_old
                # errors on MC eff are negligible
    
    else:
        eff_pseudo_nomi, d_eff_pseudo_nomi = res_pseudo_nomi.getEff(bin_key)
        eff_pseudo_alt, d_eff_pseudo_alt = res_pseudo_alt.getEff(bin_key)
        store_qnt = eff_alt*(eff_pseudo_nomi/eff_pseudo_alt) - eff_nomi
        d_store_qnt = d_eff_alt*(eff_pseudo_nomi/eff_pseudo_alt)



    plot_list.append(store_qnt), d_plot_list.append(d_store_qnt), err_band_list.append(band_amplitude)

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

print("Mean : ", sum(plot_list)/len(plot_list))
print("Std dev : ", (sum([(x - sum(plot_list)/len(plot_list))**2 for x in plot_list])/len(plot_list))**0.5)

    
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

float_points_graph = ROOT.TGraphErrors(len(bin_list_pt), 
                                  array("d", bin_list_pt), array("d", plot_list), 
                                  array("d", d_bin_list_pt), array("d", d_plot_list))

baseline_graph = ROOT.TGraphErrors(len(bin_list_pt),
                                  array("d", bin_list_pt), array("d", x_values),
                                  array("d", d_bin_list_pt), array("d", err_band_list))


lowlim, uplim = args.yLims

float_points_graph.SetMarkerStyle(20)
float_points_graph.SetMarkerSize(1.7)
float_points_graph.SetMarkerColor(ROOT.kBlack)
float_points_graph.SetLineWidth(2)
float_points_graph.SetTitle("")

float_points_graph.GetXaxis().SetTitle("Bin number")
float_points_graph.GetXaxis().SetRangeUser(-1, ntotBins+2)
float_points_graph.GetXaxis().SetTitleSize(0.04)
float_points_graph.GetXaxis().SetTitleOffset(1.1)
float_points_graph.GetXaxis().SetLabelSize(0.025)
float_points_graph.GetXaxis().SetTickSize(0.0)

# float_points_graph.GetYaxis().SetTitle("#varepsilon_{B.B.}/#varepsilon_{Bmark} - #varepsilon_{AltSig}/#varepsilon_{Nomi}")

float_points_graph.GetYaxis().SetTitle("Rel bias on " + qnt.lower() if cmp_syst_effects else qnt)
float_points_graph.GetYaxis().SetTitleSize(0.04)
float_points_graph.GetYaxis().SetRangeUser(lowlim, uplim)
float_points_graph.GetYaxis().SetTickSize(0.0)
float_points_graph.GetYaxis().SetLabelSize(0.025)
float_points_graph.GetYaxis().SetTitleOffset(1)

float_points_graph.Draw("ZAP")

baseline_graph.SetTitle("")
baseline_graph.SetLineWidth(1)
# baseline_graph.SetFillColor(4)
baseline_graph.SetFillColor(ROOT.kOrange-2)
# baseline_graph.SetFillStyle(3015)
#baseline_graph.SetFillStyle(3001)

# baseline_graph.GetYaxis().SetRangeUser(0.885, 1.045)
# baseline_graph.GetXaxis().SetRangeUser(0, 721)
baseline_graph.Draw("3 same Z")

float_points_graph.Draw("same ZP")



c.Update()


offsetXaxis = 0.5

vertline = ROOT.TLine(36, 0, 36, 1.2)
vertline.SetLineColor(11)
vertline.SetLineStyle(4)
vertline.SetLineWidth(6)

for i in range(1, nptBins): # do not need line at canvas borders 
    vertline.DrawLine(netaBins*i-offsetXaxis, lowlim, netaBins*i-offsetXaxis, uplim)


bintext = ROOT.TLatex()
#bintext.SetNDC()
bintext.SetTextSize(0.025)  # 0.03
bintext.SetTextFont(42)
bintext.SetTextAngle(45)


hline = ROOT.TLine(-1, 0, ntotBins+2, 0)
hline.SetLineColor(12)
hline.SetLineStyle(10)
hline.SetLineWidth(2)
hline.DrawLine(-1, 0, ntotBins+2, 0)


sliceLabelOffset = 6.0
    
for i in range(0, nptBins):
    ytext = lowlim + 0.003
    bintext.DrawLatex(netaBins*i + netaBins/sliceLabelOffset, ytext, bin_text[i])


float_pts_legend = copy(float_points_graph)
float_pts_legend.SetMarkerColor(ROOT.kBlack)
float_pts_legend.SetLineColor(ROOT.kBlack)
float_pts_legend.SetMarkerSize(2)
float_pts_legend.SetLineWidth(2)

baseline_legend = copy(baseline_graph)
# baseline_legend.SetMarkerColor(ROOT.kRed)
#baseline_legend.SetLineColor(ROOT.kRed)
baseline_legend.SetMarkerSize(0)
baseline_legend.SetLineWidth(0)


legend = ROOT.TLegend(0.09, 0.78, 0.25, 0.92)
legend.SetTextSize(0.03)


if analyze_scale_factors:
    if not cmp_syst_effects:
        leg_str_floating = "s_{Alt} - s_{Nomi}"
        leg_str_baseline = "#sigma_{s_{Nomi}}"
    else:
        leg_str_floating = "s_{Alt}/s_{Nomi} - s_{Alt, REF.}/s_{Nomi, REF.}"
        leg_str_baseline = "#sigma_{s_{Alt, REF.}/s_{Nomi, REF.}}"

elif cmp_syst_effects:
    leg_str_floating = "(#varepsilon_{B.B.}-#varepsilon_{Benchmark})/#varepsilon_{Benchmark}"
    leg_str_baseline = "(#varepsilon_{AltSig}-#varepsilon_{Nomi})/#varepsilon_{Nomi}"

elif check_biasMC_agreement:
    leg_str_floating = "#varepsilon_{Alt}(#varepsilon_{MC, Benchmark} / #varepsilon_{MC, Alt}) - #varepsilon_{Benchmark}"
    leg_str_baseline = "#sigma_{#varepsilon_{Benchmark}}"

else:
    leg_str_floating = "#varepsilon_{Alt} - #varepsilon_{Nomi}"
    leg_str_baseline = "#sigma_{#varepsilon_{Nomi}}"


legend.AddEntry(float_pts_legend, leg_str_floating, "lp")
legend.AddEntry(baseline_legend, leg_str_baseline, "F")
    

legend.Draw()

CMS_lumi(pad_plot, 5, 0, simulation=False)
pad_plot.Update()

c.SaveAs(f"{plot_name.replace(' ', '_').lower()}_{args.output_subscript}.png" if args.output_subscript else f"{plot_name.replace(' ', '_').lower()}.png")

