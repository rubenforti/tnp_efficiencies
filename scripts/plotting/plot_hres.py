"""
"""
import ROOT
from copy import copy
from utilities.CMS_lumi import CMS_lumi
from utilities.plot_utils import style_settings
from utilities.res_tools import results_manager
from utilities.base_lib import binning, bin_dictionary
import statistics


ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True


folder = "/scratch/rforti/tnp_efficiencies_results/tracking"
analysis = "pseudodata/ws_tracking_pseudodata.root"
cmp_flag = "legacy-onlyFailSA"

cmpBmark = False

bmark_flag = "Bmark" if cmpBmark else ""

build_results_from_txt = False


if build_results_from_txt:

    with open(folder+"/allEfficiencies.txt", "r") as file_bmark:
        row_list = file_bmark.readlines()

    res_nomi = results_manager("indep", "pt_tracking", "eta", import_txt=row_list)
    res_altSig = results_manager("indep", "pt_tracking", "eta", import_txt=row_list, altSig_check=True)

    bins_pt, bins_eta = binning("pt_tracking"), binning("eta")


    hist_dict = {
        "Efficiency": ROOT.TH2D("h_efficiency_2d", "Efficiency", len(bins_pt)-1, bins_pt, len(bins_eta)-1, bins_eta),
        "Rel. error on efficiency": ROOT.TH2D("h_efficiency_err_2d", "Rel. error on efficiency", len(bins_pt)-1, bins_pt, len(bins_eta)-1, bins_eta),
        "Rel. bias AltSig 2d": ROOT.TH2D("h_bias_altSig_2d", "Rel. bias AltSig", len(bins_pt)-1, bins_pt, len(bins_eta)-1, bins_eta),
        "Rel. bias AltSig" : ROOT.TH1D("h_bias_altSig", "Rel. bias AltSig", 25, -0.01, 0.01),
        "Pull AltSig 2d": ROOT.TH2D("h_pull_altSig_2d", "Pull AltSig", len(bins_pt)-1, bins_pt, len(bins_eta)-1, bins_eta),
        "Pull AltSig" : ROOT.TH1D("h_pull_altSig", "Pull AltSig", 25, -5, 5),
        "Ratio error AltSig/Nomi" : ROOT.TH2D("h_ratio_err_altSigNomi_2d", "Ratio error AltSig/Nomi", len(bins_pt)-1, bins_pt, len(bins_eta)-1, bins_eta)
    }


    pull_list = []
    bias_list = []
    rel_error_list = []



    for b_key, [_, idx_pt, idx_eta] in bin_dictionary("pt_tracking", "eta").items():

        eff_nomi, d_eff_nomi = res_nomi.getEff(b_key)
        eff_altSig, d_eff_altSig = res_altSig.getEff(b_key)

        hist_dict["Efficiency"].SetBinContent(idx_pt, idx_eta, eff_nomi)
        hist_dict["Rel. error on efficiency"].SetBinContent(idx_pt, idx_eta, d_eff_nomi/eff_nomi)
        hist_dict["Rel. bias AltSig 2d"].SetBinContent(idx_pt, idx_eta, (eff_altSig-eff_nomi)/eff_nomi)
        hist_dict["Pull AltSig 2d"].SetBinContent(idx_pt, idx_eta, (eff_altSig-eff_nomi)/d_eff_nomi)
        hist_dict["Ratio error AltSig/Nomi"].SetBinContent(idx_pt, idx_eta, d_eff_altSig/d_eff_nomi)

        hist_dict["Rel. bias AltSig"].Fill((eff_altSig-eff_nomi)/eff_nomi)
        hist_dict["Pull AltSig"].Fill((eff_altSig-eff_nomi)/d_eff_nomi)

        pull_list.append((eff_altSig-eff_nomi)/d_eff_nomi)
        bias_list.append((eff_altSig-eff_nomi)/eff_nomi)
        rel_error_list.append(d_eff_nomi/eff_nomi)


    file_out = ROOT.TFile(f"{folder}/current_results.root", "RECREATE")
    file_out.cd()
    [hist.Write() for hist in hist_dict.values()]
    file_out.Close()


    print("Pull:", f"Mean = {statistics.mean(pull_list)}", " ", f"Std. = {statistics.pvariance(pull_list)**0.5}")
    print("Bias:", f"Mean = {statistics.mean(bias_list)}", " ", f"Std. = {statistics.pvariance(bias_list)**0.5}")
    print("Rel. error:", f"Mean = {statistics.mean(rel_error_list)}", " ", f"Std. = {statistics.pvariance(rel_error_list)**0.5}")

else:
    # file_hres_cmp = ROOT.TFile(f"{folder}/{analysis.replace('ws', f'res_cmp-{cmp_flag}')}", "READ")
    file_hres_cmp = ROOT.TFile(f"{folder}/{analysis.replace('ws', f'res_cmpMC')}", "READ")

    hist_dict = {}
    # [histos.update({key.GetName() : file_hres.Get(key.GetName())}) for key in file_hres.GetListOfKeys()]
    [hist_dict.update({key.GetName() : file_hres_cmp.Get(key.GetName())}) for key in file_hres_cmp.GetListOfKeys()]
    # [histos.update({key.GetName() : file_other.Get(key.GetName())}) for key in file_other.GetListOfKeys()]



for hist_key, hist in hist_dict.items():

    style_settings() 

    if "2d" in hist.GetName():
        x_dim = 1600
        y_dim = 1200
        c = ROOT.TCanvas(hist_key, hist_key, x_dim, y_dim)
        c.cd()
        pad_plot = ROOT.TPad("pad_plot", "pad_plot", 0, 0, 1, 0.92)
        pad_title = ROOT.TPad("pad_title", "pad_title", 0, 0.92, 1, 1)
        pad_plot.SetMargin(0.08, 0.12, 0.12, 0.05), pad_plot.Draw()
        pad_title.SetMargin(0.15, 0.15, 0.15, 0.15), pad_title.Draw()
        plot_name = hist.GetTitle()
    else:
        ROOT.gStyle.SetOptStat("emruo")
        ROOT.gStyle.SetStatBorderSize(2)
        if hist_key == "h_efficiency": 
            ROOT.gStyle.SetStatX(0.6)
        else:
            ROOT.gStyle.SetStatX(0.9)
        ROOT.gStyle.SetStatY(0.95)
        x_dim = 1600
        y_dim = 1200
        c = ROOT.TCanvas(hist_key, hist_key, x_dim, y_dim)
        c.cd()
        pad_plot = ROOT.TPad("pad_plot", "pad_plot", 0, 0, 1, 0.95)
        pad_title = ROOT.TPad("pad_title", "pad_title", 0, 0.95, 1, 1)
        pad_plot.SetMargin(0.12, 0.08, 0.12, 0.04), pad_plot.Draw()
        pad_title.SetMargin(0., 0., 0., 0.), pad_title.Draw()
        x_title = copy(hist.GetTitle())
        plot_name = ""

    pad_plot.cd()
    hist.SetContour(51)
    hist.Draw("COLZ")
    hist.GetZaxis().SetTitle("")
    hist.SetTitleSize(0)
    hist.SetTitle("")

    x_axis, y_axis, z_axis = hist.GetXaxis(), hist.GetYaxis(), hist.GetZaxis()
    
    if "2d" in hist.GetName():
        
        h1d = hist_dict[hist_key.replace("_2d", "")]
        print(type(h1d))
        lowValZ = h1d.GetXaxis().GetXmin()
        highValZ = h1d.GetXaxis().GetXmax()
        z_axis.SetRangeUser(lowValZ, highValZ)


        x_axis.SetTitle("p_{T}^{#mu} [GeV]")
        x_axis.SetTitleOffset(1.2)
        x_axis.SetTitleSize(0.04)
        y_axis.SetTitle("#eta^{#mu}")
        y_axis.SetTitleOffset(0.8)
        y_axis.SetTitleSize(0.04)
    else:
        hist.GetXaxis().SetTitle(x_title)
        hist.GetXaxis().SetTitleOffset(1.2)
        hist.GetXaxis().SetTitleSize(0.04)
        hist.GetXaxis().SetLabelSize(0.025)
        y_axis.SetTitle("Events")
        y_axis.SetTitleOffset(1.25)
        y_axis.SetTitleSize(0.04)
        y_axis.SetLabelSize(0.025)
        #pad_plot.SetLogy()
    
    
    
    CMS_lumi(pad_plot, 5, 0)
    pad_plot.Update()
    
    pad_title.cd()
    titlebox = ROOT.TPaveText(0, 0.05, 1., 0.95, "NDC NB")
    titlebox.SetFillColor(0)
    titlebox.SetTextSize(0.6)
    titlebox.AddText(plot_name)
    titlebox.Draw()
    pad_title.Update()

    c.Update()

    if not build_results_from_txt:
        save_folder = folder+"/"+analysis.split("/")[0]
    else:
        save_folder = folder

    c.SaveAs(f"{save_folder}/{hist.GetName()}.png")

    # c.SaveAs(f"{save_folder}/{hist.GetName().replace('h', 'legacy')}.png")


