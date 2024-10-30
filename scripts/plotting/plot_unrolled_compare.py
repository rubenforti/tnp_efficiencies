import ROOT
import sys
import argparse
from utilities.base_lib import binnings, import_pdf_library, default_binnings_by_eff
from utilities.results_manager import results_manager
from utilities.CMS_lumi import CMS_lumi
from utilities.binning_utils import bin_dictionary
from utilities.plot_utils import style_settings
from array import array
from copy import copy


ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

parser = argparse.ArgumentParser(description="Plot the unrolled comparison of efficiencies or scale factors")
parser.add_argument('-i', '--input_file', type=str, help='Input file', required=True)
parser.add_argument('-b', '--benchmark_input_file', type=str, help='Benchmark input file')
parser.add_argument('-e', '--eff', type=str, help='Efficiency type', required=True)
parser.add_argument('-o', '--output', type=str, help='Output file', required=True)
parser.add_argument('-sf', '--scale_factors', action="store_true", help="Analyze scale factors")
parser.add_argument(       "--altModel_txt", type=str, default="", choices=["", "altSig", "altBkg"],
                    help='Benchmark efficiencies (in a .txt file) are selected with the specified method')
parser.add_argument('-pseudo', '--pseudodata', action="store_true", help="Analyzing pseudodata")
parser.add_argument("-ptBins", "--ptBins", nargs=2, type=int, default=[0,-1], help="Select only a specific range of pt bins to plot")
parser.add_argument("-yLims",  "--yLims", nargs=2, type=float, default=[-0.028,0.04], help="Set the y-axis limits")
parser.add_argument("-nl", "--names_legend", nargs=2, type=str, default=["Nomi", "Alt"], help="Names for the legend (in order: nominal, alternative)")
parser.add_argument(       "--compare_systematics", type=str, 
                    help="The systematic effect in the input file with the specified one in the .txt file here provided")
parser.add_argument(       "--check_biasMC_agreement", nargs=2, type=str, 
                    help="Check agreement of syst effect measured on data w.r.t. the expected bias from MC. Two files must be provided, containing the measured bias (on pseudodata) for NOMINAL and ALTERNATIVE strategies, in this order")
parser.add_argument(       "--title_subscript", type=str, default="", help="Subscript to add to the title")
args = parser.parse_args()


binning_pt, binning_eta, _ = default_binnings_by_eff[args.eff]

res_test = results_manager(args.input_file, "indep", binning_pt, binning_eta)
res_bmark = results_manager(args.benchmark_input_file, "indep", binning_pt, binning_eta, altModel_check=args.altModel_txt)

qnt = "Scale factors" if args.scale_factors else "Efficiencies"

plot_name = qnt + " comparison"


study_biasMC_agreement = True if args.check_biasMC_agreement is not None else False
study_systematics_cmp = True if args.compare_systematics is not None else False
study_diff = not (study_biasMC_agreement or study_systematics_cmp)

print(args.check_biasMC_agreement)

if args.pseudodata and not study_diff:
    sys.exit("Pseudodata can be used only for the comparison of the test and benchmark results. Exiting...")

if study_biasMC_agreement and study_systematics_cmp: 
    sys.exit("Cannot check biasMC agreement and compare systematic effects at the same time. Exiting...")

if study_systematics_cmp:
    plot_name = "Syst effects comparison on " + qnt.lower()
    oldSyst_file = args.compare_systematics
    if oldSyst_file.endswith(".txt") and hasattr(args, "altModel_txt"):
        res_old_nomi = results_manager(oldSyst_file, "indep", binning_pt, binning_eta)
        res_old_syst = results_manager(oldSyst_file, "indep",  binning_pt, binning_eta, altModel_check=args.altModel_txt)
    else:
        sys.exit("The file with the old systematic effects must be a .txt file, and the alternative method to analyize must be specified. Exiting...")

if study_biasMC_agreement:
    if args.pseudodata:
        sys.exit("Cannot check biasMC agreement when the tested effect is studied on pseudodata. Exiting...")
    plot_name = "Agreement of syst effect measured on data wrt expected bias from MC"
    if args.check_biasMC_agreement[0].endswith(".root") and args.check_biasMC_agreement[1].endswith(".root"):
        res_pseudo_nomi = results_manager(args.check_biasMC_agreement[0], "indep", binning_pt, binning_eta)
        res_pseudo_alt =  results_manager(args.check_biasMC_agreement[1], "indep", binning_pt, binning_eta)
    else:
        sys.exit("The files with the bias on pseudodata must be .root files. Exiting...")


# nptBins = len(binnings[binning_pt])-1


#Counting of pt bins starts from 0 and ends at nptBins-1
if args.ptBins[1]==-1: args.ptBins[1] = (len(binnings[binning_pt])-1) - 1

nptBins = args.ptBins[1] - args.ptBins[0] +1
netaBins = len(binnings[binning_eta])-1
ntotBins = nptBins*netaBins


c = ROOT.TCanvas("c", "c", 3600, 1800)
c.cd()

#style_settings()

pad_title = ROOT.TPad("pad_title", "pad_title", 0.0, 0.95, 1.0, 1.0)
pad_title.SetMargin(0.05, 0.05, 0.05, 0.1), pad_title.Draw()

pad_plot = ROOT.TPad("pad_plot", "pad_plot", 0.0, 0.0, 1.0, 0.95)
pad_plot.SetMargin(0.08, 0.035, 0.1, 0.05), pad_plot.Draw()

# pad_eta = ROOT.TPad("pad_plot", "pad_plot", 0.0, 0.0, 1.0, 0.netaBins)
# pad_eta.SetMargin(0.11, 0.05, 0.1, 0.05), pad_eta.Draw()


pad_title.cd()
titlebox = ROOT.TPaveText(0, 0.1, 1, 0.9, "NDC NB")
titlebox.SetFillColor(0)
titlebox.SetTextSize(0.7)
titlebox.AddText(0.5, 0.5, f"{plot_name} ({args.title_subscript})" if args.title_subscript!="" else plot_name)
titlebox.Draw()
c.Update()


bin_list_pt, d_bin_list_pt = [], []
bin_list_eta, d_bin_list_eta = [], []

plot_list, d_plot_list, err_band_list = [], [], []


bin_text = []

diff_list = [0]*netaBins

x_values = [0]*netaBins

for bin_key, [n_bin, idx_pt, idx_eta] in bin_dictionary(binning_pt, binning_eta).items():
    

    if idx_pt-1 < args.ptBins[0] or idx_pt-1 > args.ptBins[1]: continue

    #print(idx_pt, idx_eta)

    if idx_eta == 1:
        pt_text = bin_key.split("][")[0]
        pt_min, pt_max = pt_text.split("to")
        # print(pt_min, pt_max)
        pt_min = pt_min.replace("[", "")
        pt_min = pt_min.replace(".0", "")
        pt_max = pt_max.replace(".0", "")    
        bin_text.append(f"p_{{T}} #in  [{pt_min},{pt_max}] GeV/c")


    eff_test,  d_eff_test  = res_test.getEff(bin_key)
    eff_bmark, d_eff_bmark = res_bmark.getEff(bin_key) if not args.pseudodata else res_bmark.getEffMC(bin_key)

    # if d_eff_bmark/eff_bmark > 0.015: print(bin_key)

    band_amplitude = d_eff_bmark

    if study_diff:  # It is shown the difference between the two results
        if args.scale_factors:
            effMC_test, d_effMC_test = res_test.getEffMC(bin_key)
            effMC_bmark, d_effMC_bmark = res_bmark.getEffMC(bin_key)
            store_qnt, d_store_qnt = (eff_test/effMC_test) - (eff_bmark/effMC_bmark), d_eff_test/effMC_test
            band_amplitude = d_eff_bmark/effMC_bmark
        else:
            store_qnt, d_store_qnt = eff_test - eff_bmark, d_eff_test

    elif study_systematics_cmp:  # The relative difference between the two results is shown, and the uncertainty band represents the rel. difference of the old results
        '''
        if args.scale_factors:
            effMC_test, d_effMC_test = res_test.getEffMC(bin_key)
            effMC_bmark, d_effMC_bmark = res_bmark.getEffMC(bin_key)

            store_qnt, d_store_qnt = store_qnt/effMC_test, d_eff_test/effMC_test

            sf_alt, sf_nomi = eff_test/effMC_test, eff_bmark/effMC_bmark

            effMC_test_old, d_effMC_test_old = res_old_syst.getEffMC(bin_key)
            effMC_bmark_old, d_effMC_bmark_old = res_old_nomi.getEffMC(bin_key)
            sf_alt_old, sf_nomi_old = eff_test_old/effMC_test_old, eff_bmark_old/effMC_bmark_old
            store_qnt = ((sf_alt-sf_nomi)/sf_nomi)
            # x_values[idx_eta-1] = (sf_alt_old-sf_nomi_old)/sf_nomi_old
        '''
        eff_test_old, d_eff_test_old = res_old_syst.getEff(bin_key)
        eff_bmark_old, d_eff_bmark_old = res_old_nomi.getEff(bin_key)
        store_qnt = (eff_test-eff_bmark)/eff_bmark
        d_store_qnt = (eff_test/eff_bmark) * ((d_eff_test/eff_test)**2 + (d_eff_bmark/eff_bmark)**2)**0.5
        band_amplitude_syst = (eff_test_old-eff_bmark_old)/eff_bmark_old
    
    elif study_biasMC_agreement:  # The "test" results are corrected by the bias measured on pseudodata
        eff_pseudo_nomi, d_eff_pseudo_nomi = res_pseudo_nomi.getEff(bin_key)
        eff_pseudo_alt, d_eff_pseudo_alt = res_pseudo_alt.getEff(bin_key)
        biasMC_measured = eff_pseudo_alt/eff_pseudo_nomi -1
        store_qnt = (eff_test-eff_bmark) / biasMC_measured
        d_store_qnt = d_eff_test/ biasMC_measured
        band_amplitude = store_qnt * ( (d_eff_bmark/(eff_test-eff_bmark))**2 +
                                       (d_eff_pseudo_alt/eff_pseudo_alt)**2 + (d_eff_pseudo_nomi/eff_pseudo_nomi)**2 )**0.5

    
    plot_list.append(store_qnt), d_plot_list.append(d_store_qnt), err_band_list.append(band_amplitude)

    bin_list_pt.append(n_bin-1)
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
#float_points_graph.GetXaxis().SetRangeUser(-1, ntotBins+2)
float_points_graph.GetXaxis().SetRangeUser(-1 + args.ptBins[0]*48, 48*(args.ptBins[1]+1) +1)
float_points_graph.GetXaxis().SetTitleSize(0.04)
float_points_graph.GetXaxis().SetTitleOffset(1.1)
float_points_graph.GetXaxis().SetLabelSize(0.025)
float_points_graph.GetXaxis().SetTickSize(0.0)

# float_points_graph.GetYaxis().SetTitle("#varepsilon_{B.B.}/#varepsilon_{Bmark} - #varepsilon_{AltSig}/#varepsilon_{Nomi}")

float_points_graph.GetYaxis().SetTitle("Rel bias on " + qnt.lower() if study_systematics_cmp else qnt)
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


hline = ROOT.TLine(netaBins*args.ptBins[0]-1, 0, netaBins*(args.ptBins[1]+1)+1, 0)
hline.SetLineColor(12)
hline.SetLineStyle(10)
hline.SetLineWidth(2)
hline.DrawLine(netaBins*args.ptBins[0]-1, 0, netaBins*(args.ptBins[1]+1)+1, 0)


vertline = ROOT.TLine(36, 0, 36, 1.2)
vertline.SetLineColor(11)
vertline.SetLineStyle(4)
vertline.SetLineWidth(6) 
offsetXaxis = 0.5

bintext = ROOT.TLatex()
#bintext.SetNDC()
bintext.SetTextSize(0.025)  # 0.03
bintext.SetTextFont(42)
bintext.SetTextAngle(45)
sliceLabelOffset = 10       # 12.0



j=0
for i in range(args.ptBins[0], args.ptBins[1]+1):
    ytext = lowlim + 0.002
    # ytext = uplim - 0.01
    bintext.DrawLatex(netaBins*(i) + netaBins/sliceLabelOffset, ytext, bin_text[j])
    if i!=args.ptBins[1]:
        vertline.DrawLine(netaBins*(i+1)-offsetXaxis, lowlim, netaBins*(i+1)-offsetXaxis, uplim)
    # print(i, netaBins*(i+1)-offsetXaxis)
    j+=1


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


# legend = ROOT.TLegend(0.8, 0.8, 0.96, 0.94)  # Upper-right corner
legend = ROOT.TLegend(0.09, 0.8, 0.25, 0.94)  # Upper-left corner
# legend = ROOT.TLegend(0.8, 0.12, 0.96, 0.28)  # Lower-right corner
legend.SetTextSize(0.03)


in_var = "#varepsilon" if not args.scale_factors else "s"
nl = args.names_legend
str_nomi, str_alt = in_var+"_{"+nl[0]+"}", in_var+"_{"+nl[1]+"}"

if study_diff and not args.pseudodata:
    leg_str_floating = f"{str_alt} - {str_nomi}"
    leg_str_baseline = f"#sigma({str_nomi})"
elif args.pseudodata:
    leg_str_floating = in_var+"_{Fit} - "+in_var+"_{True}"
    leg_str_baseline = "#sigma("+in_var+"_{True})"
elif study_systematics_cmp:
    leg_str_floating = f"({str_alt}/{str_nomi}) - 1"
    leg_str_baseline = in_var+"_{Alt, REF.}/"+in_var+"_{Nomi, REF.} - 1"
elif study_biasMC_agreement:
    ### HAS TO BE CORRECTED!!!!
    leg_str_floating = "#varepsilon_{Alt}(#varepsilon_{MC, Benchmark} / #varepsilon_{MC, Alt}) - #varepsilon_{Benchmark}"
    leg_str_baseline = "#sigma_{#varepsilon_{Benchmark}}"



legend.AddEntry(float_pts_legend, leg_str_floating, "lp")
legend.AddEntry(baseline_legend, leg_str_baseline, "F")

legend.Draw()

CMS_lumi(pad_plot, 5, 0, simulation=False)
pad_plot.Update()

writename = "prova"

if (args.ptBins[0]!=0 or args.ptBins[1]!=(len(binnings[binning_pt])-2)):
    writename+=f"_ptBins{args.ptBins[0]}to{args.ptBins[1]}"

c.SaveAs(writename+".png")

