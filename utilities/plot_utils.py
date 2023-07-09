"""
"""
import ROOT
import sys
import os
import pickle
from array import array
from utilities.results_utils import results_manager
from utilities.base_library import lumi_factors, binning, sumw2_error
from utilities.CMS_lumi import CMS_lumi

bkg_categories= ["WW", "WZ", "ZZ", "TTSemileptonic", "Ztautau"]


def style_settings():
    """
    General style settings for plots
    """
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetStatBorderSize(0)
    ROOT.gStyle.SetStatX(1.)
    ROOT.gStyle.SetStatY(1.)

###############################################################################

def plot_minv_distrib_with_fit(pad, flag, axis, data, pdf_fit, pull=False):
    """
    """

    pad.cd()

    if not pull: 
        pad.Divide(1, 2, 0, 0)
        subpad_title = pad.GetPad(1)
        subpad_title.SetPad("subpad_title", "subpad_title", 0, 0.9, 1, 1, 0, 0, 0)
        subpad_title.SetMargin(0.15, 0.05, 0.05, 0), subpad_title.Draw()
        subpad_plot = pad.GetPad(2)
        subpad_plot.SetPad("subpad_plot", "subpad_plot", 0, 0, 1, 0.9, 0, 0, 0)
        subpad_plot.SetMargin(0.15, 0.05, 0.15, 0.05), subpad_plot.Draw()
    else:
        pad.Divide(1, 3, 0, 0)
        subpad_title = pad.GetPad(1)
        subpad_title.SetPad("subpad_title", "subpad_title", 0, 0.9, 1, 1, 0, 0, 0)
        subpad_title.SetMargin(0.15, 0.05, 0.05, 0), subpad_title.Draw()
        subpad_plot = pad.GetPad(2)
        subpad_plot.SetPad("subpad_plot", "subpad_plot", 0, 0.3, 1, 0.9, 0, 0, 0)
        subpad_plot.SetMargin(0.15, 0.05, 0, 0.05), subpad_plot.Draw()
        subpad_pull = pad.GetPad(3)
        subpad_pull.SetPad("subpad_pull", "subpad_pull", 0, 0, 1, 0.3, 0, 0, 0)
        subpad_pull.SetMargin(0.15, 0.05, 0.15, 0), subpad_pull.Draw()

    
    # In questa maniera il titolo non viene disegnatom in realtà manca proprio la box che dovrebbe contenere
    # il testo. Forse si può bypassare il problema con la legenda
    subpad_title.cd()
    titlebox = ROOT.TPaveText(0, 0, 1, 1, "NDC NB")
    titlebox.SetFillColor(5)
    titlebox.SetTextSize(0.5)
    titlebox.SetTextAlign(12)
    titlebox.AddText(pad.GetTitle())
    titlebox.Draw()
    subpad_title.Update()
    subpad_title.Draw()
    pad.Update()


    subpad_plot.cd()
    plot_frame = axis.frame(ROOT.RooFit.Bins(axis.getBins("plot_binning")))
    
    # COSÌ COMPARE ALL'INTERNO DEL PLOT, CHE PALLE
    # BISOGNERÀ CREARE UN TPAVETEXT DA PIAZZARE BENE NEL PAD COME TITOLO
    plot_frame.SetTitle("")
    plot_frame.SetTitleSize(0)
    frame_yaxis = plot_frame.GetYaxis()
    frame_yaxis.SetTitleOffset(1.2)

    if pull:
        frame_axis = plot_frame.GetXaxis()
        frame_axis.SetLabelSize(0)
        frame_axis.SetTitleSize(0)

    data.plotOn(plot_frame,
                ROOT.RooFit.Binning("plot_binning"),
                ROOT.RooFit.Name("Data"))

    pdf_fit.plotOn(plot_frame, 
                   ROOT.RooFit.LineColor(ROOT.kRed),
                   ROOT.RooFit.Name("Fit"))
        
    for comp in pdf_fit.getComponents():
        if "bkg" in comp.GetName():
            bkg_set = ROOT.RooArgSet(comp)
            pdf_fit.plotOn(plot_frame,
                ROOT.RooFit.Components(bkg_set),
                ROOT.RooFit.LineColor(ROOT.kBlue),
                ROOT.RooFit.LineStyle(ROOT.kDashed))
    plot_frame.Draw()

    CMS_lumi(subpad_plot, 5, 0, simulation=True)
    pad.Update() 

    if pull:
        subpad_pull.cd()
        pull_frame = plot_frame.pullHist("Data", "Fit", False)
        pull_frame.SetTitle("")
        pull_axis = pull_frame.GetXaxis()
        pull_axis.SetTitle(axis.GetTitle())
        pull_axis.SetTitleSize(0.07)
        pull_axis.SetLabelSize(0.07)
        pull_axis.SetTitleOffset(1)
        pull_frame.Draw()
        pad.Update()

###############################################################################

def plot_pass_fail(type_analysis, plot_objects, bin_key, pull=False, figpath=""): 
    """
    """
    c = ROOT.TCanvas(f"pf_plot_{bin_key}", f"pf_plot_{bin_key}", 900, 1200)
    c.cd()

    style_settings()

    pad_title_eff = ROOT.TPad("pad_title_eff", "pad_title_eff", 0, 0.9, 1, 1)
    pad_plot_pass = ROOT.TPad("pad_plot_pass", "Passing probes", 0, 0.3, 0.5, 0.9)
    pad_plot_fail = ROOT.TPad("pad_plot_fail", "Failing probes", 0.5, 0.3, 1, 0.9)
    pad_info_pass = ROOT.TPad("pad_info_pass", "pad_info_pass", 0, 0, 0.5, 0.3)
    pad_info_fail = ROOT.TPad("pad_info_fail", "pad_info_fail", 0.5, 0, 1, 0.3)


    for pad_core in [pad_plot_pass, pad_plot_fail, pad_info_pass, pad_info_fail]:
        pad_core.SetMargin(0.15, 0.05, 0.15, 0.05)
        pad_core.Draw()


    pad_title_eff.cd()
    main_info_box = ROOT.TPaveText(0, 0, 1, 1, "NDC NB")
    main_info_box.SetFillColor(0)
    # titlebox.SetTextFont(42)
    main_info_box.SetTextSize(0.35)
    main_info_box.AddText(f"Bin {bin_key}")
    eff, deff = plot_objects["efficiency"]
    main_info_box.AddText(f"Efficiency = {round(eff*100,2)} #pm {round(deff*100,2)} %")
    main_info_box.Draw()
    c.Update()

    pass_obj = plot_objects["pass"]
    fail_obj = plot_objects["fail"]

    pass_obj["axis"].setBins(60, "plot_binning")
    fail_obj["axis"].setBins(60, "plot_binning")


    plot_minv_distrib_with_fit(pad_plot_pass, "pass", pass_obj["axis"], pass_obj["data"], pass_obj["model"], pull=pull)
    plot_minv_distrib_with_fit(pad_plot_fail, "fail", fail_obj["axis"], fail_obj["data"], fail_obj["model"], pull=pull)
    c.SaveAs(f"{figpath}/fit_pf_{type_analysis}_{bin_key}.pdf")

###############################################################################

def plot_eff_results(filename, results, binning_pt=(), binning_eta=()):
    """
    """

    file = ROOT.TFile.Open(filename, "UPDATE")

    if len(binning_pt) == 0:
        binning_pt = binning("pt")
    if len(binning_eta) == 0:
        binning_eta = binning("eta")

    ROOT.gStyle.SetOptStat("n")

    res_dict = results.dictionary()

    idx_pt, idx_eta = 0, 0

    h_eff = ROOT.TH2D("efficiency_th2", "Efficiency (pt,eta)", len(binning_pt)-1, binning_pt,
                      len(binning_eta)-1, binning_eta)
    h_deff = ROOT.TH2D(
        "eff_rel_error_th2", "Relative error on efficiency (pt,eta)",
        len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta)
    
    '''
    if results._analysis == 'indep':
        corr_pass = ROOT.TH2
    '''

    for key in res_dict.keys():
        idx_pt = int(key.split(",")[0])
        idx_eta = int(key.split(",")[1])
        eff = res_dict[key]["efficiency"]
        # print(eff[0], eff[1])
        h_eff.SetBinContent(idx_pt, idx_eta, eff[0])
        h_eff.SetBinError(idx_pt, idx_eta, eff[1])
        h_deff.SetBinContent(idx_pt, idx_eta, eff[1]/eff[0])

    c1 = ROOT.TCanvas("efficiency", "efficiency", 1600, 900)
    c1.cd()
    ROOT.gPad.SetRightMargin(0.15)
    h_eff.Draw("COLZ")
    c1.SaveAs("figs/efficiency_sim.pdf")

    c2 = ROOT.TCanvas("error efficiency", "error efficiency", 1600, 900)
    c2.cd()
    ROOT.gPad.SetRightMargin(0.15)
    h_deff.Draw("COLZ")
    c2.SaveAs("figs/efficiency_sim_rel_errors.pdf")

    #file.Write("efficiency_th2")
    # file.Write("eff_rel_error_th2")
    file.Close()

###############################################################################

def plot_bkg_on_histo(plot_objects, flag, bin_key, figpath=''):
    """
    """

    imported_data, imported_mc_sig, imported_pdf_bkg = False, False, False

    c = ROOT.TCanvas(f"c_{bin_key}_{flag}", "c", 900, 1200)
    c.cd()

    style_settings()

    pad_title = ROOT.TPad("pad_title", "pad_title", 0, 0.95, 1, 1)
    pad_plot = ROOT.TPad("pad_plot", "pad_plot", 0, 0.25, 0.7, 0.95)
    pad_pull = ROOT.TPad("pad_pull", "pad_pull", 0, 0, 0.7, 0.25)
    pad_info = ROOT.TPad("pad_info", "pad_info", 0.7, 0, 1, 0.9)

    pad_title.SetMargin(0.1, 0.1, 0.1, 0.1), pad_title.Draw()
    pad_plot.SetMargin(0.15, 0.05, 0, 0.05), pad_plot.Draw()
    pad_pull.SetMargin(0.15, 0.05, 0.2, 0), pad_pull.Draw()
    pad_info.SetMargin(0, 0, 0.05, 0.05), pad_info.Draw()

    pad_title.cd()
    titlebox = ROOT.TPaveText(0, 0, 1, 1, "NDC NB")
    titlebox.SetFillColor(0)
    # titlebox.SetTextFont(42)
    titlebox.SetTextSize(0.35)
    titlebox.AddText(f"Bin {bin_key} - {flag}ing probes")
    titlebox.Draw()
    c.Update()

    axis = plot_objects.pop("axis")


    pad_plot.cd()
    # pad_plot.SetLogy()
    frame = axis.frame(ROOT.RooFit.Bins(60))
    frame.SetTitle("")
    frame.SetTitleSize(0)
    frame_yaxis = frame.GetYaxis()
    frame_yaxis.SetTitle("Events / (1 GeV)")
    frame_axis = frame.GetXaxis()
    frame_axis.SetLabelSize(0)
    xtitle = frame_axis.GetTitle()
    frame_axis.SetTitleSize(0)

    if "fit_pars" in plot_objects.keys():
        fit_pars = plot_objects.pop("fit_pars")

    if "data" in plot_objects.keys():
        imported_data = True
        datahist = plot_objects.pop("data")
        datahist.plotOn(frame, 
                        ROOT.RooFit.Binning("plot_binning"),
                        ROOT.RooFit.Name("Data"))
    
    if "MC_signal" in plot_objects.keys():
        imported_mc_sig = True
        mc_sig_dict = plot_objects.pop("MC_signal")

        mc_sig_dict["histo_pdf"].plotOn(frame,
                                        ROOT.RooFit.Name("MC signal"),
                                        # ROOT.RooFit.Binning("plot_binning"),
                                        ROOT.RooFit.LineColor(mc_sig_dict["color"]),
                                        ROOT.RooFit.Normalization(mc_sig_dict["integral"], ROOT.RooAbsReal.NumEvent))

    if "pdf_bkg_fit" in plot_objects.keys():
        imported_pdf_bkg = True
        pdf_bkg_obj = plot_objects.pop("pdf_bkg_fit")
        pdf_bkg_obj["pdf"].plotOn(frame, 
                                ROOT.RooFit.Name("Fitted bkg"),
                                ROOT.RooFit.LineColor(pdf_bkg_obj["color"]),
                                ROOT.RooFit.Normalization(pdf_bkg_obj["norm"].getVal(), ROOT.RooAbsReal.NumEvent))

    total_bkg = plot_objects.pop("total_bkg")
    total_bkg["roohisto"].plotOn(frame,
                                 ROOT.RooFit.Name("Total bkg"),
                                 ROOT.RooFit.Binning("plot_binning"),
                                 # ROOT.RooFit.Invisible(),
                                 ROOT.RooFit.LineColor(total_bkg["color"]),
                                 ROOT.RooFit.MarkerColor(total_bkg["color"]),
                                 # ROOT.RooFit.Normalization(total_bkg["integral"], ROOT.RooAbsReal.NumEvent)
                                 )
    
    pad_info.cd()
    textbox = ROOT.TPaveText(0, 0.4, 1, 0.8, "NDC NB")
    textbox.SetFillColor(0)
    textbox.SetTextSize(0.05)
    textbox.SetTextAlign(12)
    
    sigma_histo_bkg = sumw2_error(total_bkg['roohisto'])
    if imported_data:
        textbox.AddText(f"Data entries: {datahist.sumEntries()}")
    if imported_mc_sig:
        sigma_histo_signal = sumw2_error(mc_sig_dict['roohisto'])
        textbox.AddText(f"MC signal entries: {round(mc_sig_dict['integral'],2)} #pm {round(sigma_histo_signal,2)}") 
    if imported_pdf_bkg:
        textbox.AddText(f"Nbkg from fit: {round(pdf_bkg_obj['norm'].getVal(),2)} #pm {round(pdf_bkg_obj['norm'].getError(),2)}")

    textbox.AddText(f"Nbkg from MC: {round(total_bkg['integral'],2)} #pm {round(sigma_histo_bkg,2)}")
        
    if imported_pdf_bkg:
        delta_bkg = total_bkg['integral'] - pdf_bkg_obj['norm'].getVal()
        sigma_on_delta_bkg = (sigma_histo_bkg**2 + pdf_bkg_obj['norm'].getError()**2)**0.5
        nsigma = delta_bkg/sigma_on_delta_bkg
        textbox.AddText(f"Distance in #sigma = {round(nsigma, 2)}")
    
    textbox.AddText("--------------------")
    c.Update()

    print(plot_objects.keys())
    for bkg_key in plot_objects:
        pad_plot.cd()
        bkg_obj = plot_objects[bkg_key]
        bkg_obj["histo_pdf"].plotOn(frame, 
                         ROOT.RooFit.Name(bkg_key),
                         ROOT.RooFit.LineColor(bkg_obj["color"]),
                         ROOT.RooFit.LineStyle(1),
                         ROOT.RooFit.Normalization(bkg_obj["integral"], ROOT.RooAbsReal.NumEvent))
        pad_info.cd()
        bkg_error = sumw2_error(bkg_obj["roohisto"])
        textbox.AddText(f"Num {bkg_key} = {round(bkg_obj['integral'],2)} #pm {round(bkg_error, 2)}")

    pad_plot.cd()
    frame.Draw()
    CMS_lumi(pad_plot, 5, 0, simulation=True)

    if imported_pdf_bkg:
        pad_pull.cd()
        pull = frame.pullHist("Total bkg", "Fitted bkg", True)
        pull.SetTitle("")
        pull_axis = pull.GetXaxis()
        pull_axis.SetTitle(xtitle)
        pull_axis.SetTitleSize(0.07)
        pull_axis.SetLabelSize(0.07)
        pull_axis.SetTitleOffset(1)
        pull.Draw()

    pad_info.cd()
    textbox.Draw()  
    c.Update()

    #print(total_bkg["integral"], pdf_bkg_obj["norm"].getVal())

    c.SaveAs(f"{figpath}/bkg_{bin_key}_{flag}.pdf")

        





if __name__ == '__main__':

    
    '''
    file = ROOT.TFile('results/benchmark_iso/ws_iso_indep.root', "READ")

    ws = file.Get("w")
    
    res = results_manager("indep")
    res.Open("results/benchmark_iso_sim/results_iso_sim.pkl")

    for bin_pt in range(1, 16):
        for bin_eta in range(1, 49):
            res.add_result(ws, bin_pt, bin_eta, check=False)
    '''
    # res.Write('results/benchmark_iso/new_results_2.pkl')

    # plot_results("root_files/ws/ws_iso_indep.root", res)

    style_settings()


    axis = ROOT.RooRealVar("mass", "mass", 60, 120)

    mu = ROOT.RooRealVar("mu", "mu", 90, 60, 120)
    sigma = ROOT.RooRealVar("sigma", "sigma", 5, 0.5, 10)
    gaus = ROOT.RooGaussian("gaus", "gaus", axis, mu, sigma)

    data = gaus.generateBinned(ROOT.RooArgSet(axis), 10000)

    c = ROOT.TCanvas("c", "c", 800, 600)
    c.cd()

    pad_plot = ROOT.TPad("pad_plot", "pad_plot", 0, 0.25, 0.7, 0.95)
    pad_plot.Draw()
    plot_minv_distrib_with_fit(pad_plot, axis, data, gaus, pull=True)
    
    c.SaveAs("test.pdf")

