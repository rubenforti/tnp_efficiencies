"""
"""
import ROOT
import sys
import os
import pickle
from copy import copy
from array import array
from utilities.base_library import lumi_factors, binning, sumw2_error
from utilities.CMS_lumi import CMS_lumi

bkg_categories= ["WW", "WZ", "ZZ", "TTSemileptonic", "Ztautau", "SameCharge"]

colors = { 
    "WW_bkg" : ROOT.kOrange+1,
    "WZ_bkg" : ROOT.kYellow+3,
    "ZZ_bkg" : ROOT.kGreen+1,
    "TTSemileptonic_bkg" : ROOT.kCyan+1,
    "Ztautau_bkg" : ROOT.kMagenta+1,
    "SameCharge_bkg" : ROOT.kOrange+10,
    "total_bkg" : ROOT.kRed,
    "pdf_bkg_fit" : ROOT.kRed,
    "signal" : ROOT.kBlue
}

def style_settings():
    """
    General style settings for plots
    """
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetStatBorderSize(0)
    ROOT.gStyle.SetStatX(1.)
    ROOT.gStyle.SetStatY(1.)

###############################################################################

def plot_minv_distrib_with_fit(pad, axis, data, pdf_fit):
    """
    Function that creates and saves the plot containing the invariant mass 
    distribution, together with the post-fit pdf. The pull is not implemented
    yet, since there has been some memory problems with the pads.
    """

    pad.cd()

    plot_frame = axis.frame(ROOT.RooFit.Bins(axis.getBins("plot_binning")))
    plot_frame.SetTitle("")
    # plot_frame.SetTitleSize(0)
    frame_yaxis = plot_frame.GetYaxis()
    frame_yaxis.SetTitleOffset(1.2)
    frame_yaxis.SetTitle("Events / (1 GeV)")
    # frame_yaxis.SetTitleSize(0.08)
    frame_axis = plot_frame.GetXaxis()
    frame_axis.SetTitle("M_{inv} TP [GeV]")
    # frame_axis.SetTitleSize(0.05)
    frame_axis.SetTitleOffset(1.2)


    data.plotOn(plot_frame,
                ROOT.RooFit.Binning("plot_binning"),
                ROOT.RooFit.Name("Data"))
    data.statOn(plot_frame,
                ROOT.RooFit.Label(pad.GetTitle()),
                ROOT.RooFit.What("H"),
                ROOT.RooFit.Layout(0.68, 0.95, 0.95))

    pdf_fit.plotOn(plot_frame,
                   ROOT.RooFit.Name("Fit"),
                   ROOT.RooFit.NormRange("fitRange"),
                   ROOT.RooFit.Range("fitRange"),
                   ROOT.RooFit.LineColor(ROOT.kRed))
    
    for comp in pdf_fit.getComponents():
        if "bkg" in comp.GetName():
            bkg_set = ROOT.RooArgSet(comp)
            pdf_fit.plotOn(plot_frame,
                ROOT.RooFit.NormRange("fitRange"),
                ROOT.RooFit.Name("Bkg pdf"),
                ROOT.RooFit.Components(bkg_set),
                ROOT.RooFit.LineColor(ROOT.kBlue),
                ROOT.RooFit.LineStyle(ROOT.kDashed),
                ROOT.RooFit.LineWidth(4))

    plot_frame.Draw()

    CMS_lumi(pad, 5, 0, simulation=False)
    pad.Update()
    

###############################################################################

def plot_fitted_pass_fail(type_analysis, plot_objects, bin_key, pull=False, figpath=""): 
    """
    Function that creates and saves the plot containing the passing and failing
    probes histograms, together with the post-fit pdfs and the efficiency 
    inferred. The objects to be plotted are contained in the "plot_objects" 
    dictionary, which must be structured as follows:
        plot_objects = { "pass(fail)" : { "axis" : RooRealVar,  "data" : RooDataHist, 
                                          "model" : RooAbsPdf,  "res" : RooFitResult },
                         "efficiency" : [eff, deff],
                         "efficiency_mc" : [eff_mc, deff_mc],
                         "scale_factor" : [scale_factor, dscale_factor]
                        }
    """
    c = ROOT.TCanvas(f"pf_plot_{bin_key}", f"pf_plot_{bin_key}", 1600, 1200)
    c.cd()

    style_settings()

    main_info_edge = 0.865
    plot_info_edge = 0.285

    pad_title_eff = ROOT.TPad("pad_title_eff", "pad_title_eff", 0.35, 0.93, 0.65, 0.995)
    pad_main_info = ROOT.TPad("pad_main_info", "main_info", 0.1, main_info_edge, 0.9, 0.93)
    pad_plot_pass = ROOT.TPad("pad_plot_pass", "Passing probes", 0, plot_info_edge, 0.5, main_info_edge)
    pad_plot_fail = ROOT.TPad("pad_plot_fail", "Failing probes", 0.5, plot_info_edge, 1, main_info_edge)
    pad_info = ROOT.TPad("pad_info", "pad_info", 0, 0, 1, plot_info_edge)

    
    pad_title_eff.SetMargin(0.15, 0.15, 0.15, 0.15)
    pad_title_eff.Draw()

    pad_main_info.SetMargin(0.15, 0.15, 0.15, 0.15)
    pad_main_info.Draw()

    pad_plot_pass.SetMargin(0.15, 0.05, 0.1, 0.05)
    pad_plot_pass.Draw()

    pad_plot_fail.SetMargin(0.15, 0.05, 0.1, 0.05)
    pad_plot_fail.Draw()

    pad_info.SetMargin(0.15, 0.05, 0.05, 0.05), 
    pad_info.Draw()
    
    pad_title_eff.cd()
    titlebox = ROOT.TPaveText(0, 0.1, 1, 0.9, "NDC NB")
    titlebox.SetFillColor(33)
    # titlebox.SetTextFont(42)
    titlebox.SetTextSize(0.3)
    titlebox.AddText(0.5, 0.5, f"Bin {bin_key}")
    titlebox.Draw()
    c.Update()

    pad_main_info.cd()
    main_info_box = ROOT.TPaveText(0, 0.1, 1., 0.9, "NDC NB")
    main_info_box.SetFillColor(33)
    # titlebox.SetTextFont(42)
    main_info_box.SetTextSize(0.35)
    eff, deff = plot_objects["efficiency"]
    eff_mc, deff_mc = plot_objects["efficiency_mc"]
    scale_factor, dscale_factor = plot_objects["scale_factor"]
    main_info_box.AddText(0.2, 0.5, f"#varepsilon = {(eff*100):.2f} #pm {(deff*100):.2f} %")
    main_info_box.AddText(0.5, 0.5, f"#varepsilon MC = {(eff_mc*100):.2f} #pm {(deff_mc*100):.2f} %")
    main_info_box.AddText(0.8, 0.5, f"SF = {scale_factor:.4f} #pm {dscale_factor:.4f}")
    main_info_box.Draw()
    c.Update()

    pass_obj = plot_objects["pass"]
    fail_obj = plot_objects["fail"]

    binning_pass = pass_obj["axis"].getBinning()
    binning_fail = fail_obj["axis"].getBinning()
    
    pass_obj["axis"].setBins(binning_pass.numBins(), "plot_binning")
    fail_obj["axis"].setBins(binning_fail.numBins(), "plot_binning")

    plot_minv_distrib_with_fit(pad_plot_pass, pass_obj["axis"], pass_obj["data"], pass_obj["model"])
    plot_minv_distrib_with_fit(pad_plot_fail, fail_obj["axis"], fail_obj["data"], fail_obj["model"])


    pad_info.cd()
    uplim_stats_gen = 0.05 if type_analysis == "indep" else 0.6
    stats_pass = ROOT.TPaveText(0, uplim_stats_gen, 0.5, 0.95, "NDC NB")
    stats_fail = ROOT.TPaveText(0.5, uplim_stats_gen, 1., 0.95, "NDC NB")
    stats_gen = ROOT.TPaveText(0, 0.05, 1., 0.6, "NDC NB")
    stats_pass.SetFillColor(0), stats_pass.SetTextSize(0.07)
    stats_fail.SetFillColor(0), stats_fail.SetTextSize(0.07)
    stats_gen.SetFillColor(0), stats_gen.SetTextSize(0.07)

    if type_analysis == "indep":
        [stats_pass.AddText(f"{par.GetTitle()} = {par.getVal():.3f} #pm {par.getError():.3f}") 
         for par in pass_obj["res"].floatParsFinal()]
        stats_pass.AddText("------------------------------")
        stats_pass.AddText(f"Status = {pass_obj['res'].status()}")
        stats_pass.AddText(f"Cov quality = {pass_obj['res'].covQual()}")
        stats_pass.AddText(f"Edm = {pass_obj['res'].edm():.5f}")
        ndof = pass_obj["data"].numEntries() - pass_obj["res"].floatParsFinal().getSize()
        chi2 = float(pass_obj["res"].GetTitle())
        stats_pass.AddText(f"Chi2/ndof = {chi2:.1f} / {ndof}")
        stats_pass.Draw()
        c.Update()

        [stats_fail.AddText(f"{par.GetTitle()} = {par.getVal():.3f} #pm {par.getError():.3f}") 
         for par in fail_obj["res"].floatParsFinal()]
        stats_fail.AddText("------------------------------")
        # stats_box.AddText(f"Chi2/ndof = {res.chi2(), 2) / pars.getSize()}")
        stats_fail.AddText(f"Status = {fail_obj['res'].status()}")
        stats_fail.AddText(f"Cov quality = {fail_obj['res'].covQual()}")
        stats_fail.AddText(f"Edm = {fail_obj['res'].edm():.5f}")
        ndof = fail_obj["data"].numEntries() - fail_obj["res"].floatParsFinal().getSize()
        chi2 = float(fail_obj["res"].GetTitle())
        stats_fail.AddText(f"Chi2/ndof = {chi2:.1f} / {ndof}")
        stats_fail.Draw()
        c.Update()

    else:
        for par in pass_obj["res"].floatParsFinal():
            if "pass" in par.GetName():
                stats_pass.AddText(f"{par.GetTitle()} = {par.getVal():.2f} #pm {par.getError():.2f}")
            elif "fail" in par.GetName():
                stats_fail.AddText(f"{par.GetTitle()} = {par.getVal():.2f} #pm {par.getError():.2f}")
            else:
                if "efficiency" in par.GetName():
                    stats_gen.AddText(f"{par.GetTitle()} = {par.getVal():.4f} #pm {par.getError():.4f}")
                else:
                    stats_gen.AddText(f"{par.GetTitle()} = {par.getVal():.2f} #pm {par.getError():.2f}")

        stats_gen.AddText("------------------------------------------------------------")
        stats_gen.AddText(f"Status = {pass_obj['res'].status()}")
        stats_gen.AddText(f"Cov quality = {pass_obj['res'].covQual()}")
        stats_gen.AddText(f"Edm = {pass_obj['res'].edm():.5f}")
        ndof = pass_obj["data"].numEntries() + fail_obj["data"].numEntries() - \
            pass_obj["res"].floatParsFinal().getSize()
        chi2 = float(pass_obj["res"].GetTitle())
        stats_gen.AddText(f"Chi2/ndof = {chi2:.1f} / {ndof}")
        stats_pass.Draw(), stats_fail.Draw(), stats_gen.Draw()
        c.Update()

    c.SaveAs(f"{figpath}/fit_pf_{type_analysis}_{bin_key}.pdf")

###############################################################################

def plot_bkg(plot_dictionary, flag, bin_key, logscale=True, figpath=''):
    """
    Function that creates and saves the plot containing the mass distributions 
    of the various bkg processes and the total bkg; a reference histogram 
    (data or signal MC) can be plotted in the same figure. For graphical 
    choice, the bkg distributions of single processes (and the signal MC 
    distribution) are plotted as RooHistPdf.
    The objects to be plotted are contained in the dictionary "plot_objects", 
    which must be structured as follows:
        plot_objects = { 
            "axis" : RooRealVar,
            "total_bkg" : { "roohisto":RooDataHist, "lumi_scale":float, "integral":float},
            "{bkg_process}_bkg" : { "roohisto":RooDataHist, "histo_pdf":RooHistPdf,
                                    "lumi_scale":float, "integral":float},
            "MC_signal (optional)" : as above,
            "data (optional)" : RooDataHist,
            "fit_pars (optional)" : { "nsig":RooRealVar, "nbkg":RooRealVar },
            "pdf_bkg_fit (optional)" : {"pdf":RooAbsPdf, "norm":RooRealVar}
            }
        }
    """

    imported_data, imported_mc_sig, imported_pdf_bkg = False, False, False

    c = ROOT.TCanvas(f"c_{bin_key}_{flag}", "c", 1600, 900)
    c.cd()

    style_settings()

    pad_title = ROOT.TPad("pad_title", "pad_title", 0, 0.9, 1, 1)
    pad_plot = ROOT.TPad("pad_plot", "pad_plot", 0, 0, 0.6, 0.9)
    # pad_pull = ROOT.TPad("pad_pull", "pad_pull", 0, 0, 0.65, 0.25)
    pad_info = ROOT.TPad("pad_info", "pad_info", 0.6, 0.1, 1, 0.9)

    pad_title.SetMargin(0.1, 0.1, 0.1, 0.1), pad_title.Draw()
    pad_plot.SetMargin(0.12, 0.05, 0.12, 0.05), pad_plot.Draw()
    pad_info.SetMargin(0.05, 0.05, 0.1, 0.9), pad_info.Draw()

    pad_title.cd()
    titlebox = ROOT.TPaveText(0, 0, 1, 1, "NDC NB")
    titlebox.SetFillColor(0)
    # titlebox.SetTextFont(42)
    titlebox.SetTextSize(0.35)
    titlebox.AddText(f"Bin {bin_key} - {flag}ing probes")
    titlebox.Draw()
    c.Update()

    plot_objects = copy(plot_dictionary)

    axis = plot_objects.pop("axis")

    pad_plot.cd()
    frame = axis.frame(ROOT.RooFit.Bins(axis.getBins("plot_binning")))
    frame.SetTitle("")
    frame.SetTitleSize(0)
    frame_yaxis = frame.GetYaxis()
    frame_yaxis.SetTitle("Events / (1 GeV)")
    frame_yaxis.SetTitleSize(0.035)
    frame_yaxis.SetTitleOffset(1.2)
    frame_axis = frame.GetXaxis()
    # xtitle = frame_axis.GetTitle(), frame_axis.SetTitleSize(0)
    frame_axis.SetTitle("M_{inv} TP [GeV]")
    frame_axis.SetTitleSize(0.035)
    frame_axis.SetTitleOffset(1.2)

    ctrl_plot_max = 0


    if "fit_pars" in plot_objects.keys():
        fit_pars = plot_objects.pop("fit_pars")

    if "data" in plot_objects.keys():
        imported_data = True
        datahist = plot_objects.pop("data")
        datahist.plotOn(frame, 
                        ROOT.RooFit.Binning("plot_binning"),
                        ROOT.RooFit.Name("Data"))
        for bin_idx in range(datahist.numEntries()):
            if datahist.weight(bin_idx) > ctrl_plot_max:
                ctrl_plot_max = datahist.weight(bin_idx)
    
    if "MC_signal" in plot_objects.keys():
        imported_mc_sig = True
        mc_sig_dict = plot_objects.pop("MC_signal")
        mc_sig_dict["histo_pdf"].plotOn(frame,
                                        ROOT.RooFit.Name("MC signal"),
                                        ROOT.RooFit.LineColor(colors["signal"]),
                                        ROOT.RooFit.Normalization(mc_sig_dict["integral"], ROOT.RooAbsReal.NumEvent))
        for bin_idx in range(mc_sig_dict["roohisto"].numEntries()):
            if mc_sig_dict["roohisto"].weight(bin_idx) > ctrl_plot_max:
                ctrl_plot_max = mc_sig_dict["roohisto"].weight(bin_idx)

    if "pdf_bkg_fit" in plot_objects.keys():
        imported_pdf_bkg = True
        pdf_bkg_obj = plot_objects.pop("pdf_bkg_fit")
        pdf_bkg_obj["pdf"].plotOn(frame, 
                                ROOT.RooFit.Name("Fitted bkg"),
                                ROOT.RooFit.LineColor(colors["pdf_bkg_fit"]),
                                ROOT.RooFit.Normalization(pdf_bkg_obj["norm"].getVal(), ROOT.RooAbsReal.NumEvent))


    total_bkg = plot_objects.pop("total_bkg")
    
    for nbins_total_bkg in [60, 30, 20]:
            
            ctrl_plot_binning = 0
            axis.setBins(nbins_total_bkg, f"plot_binning_total_bkg")
            tmp_histo = ROOT.RooDataHist(total_bkg["roohisto"].GetName(), total_bkg["roohisto"].GetTitle(), 
                                            ROOT.RooArgSet(axis), f"plot_binning_total_bkg")
            tmp_histo.add(total_bkg["roohisto"])

            for i in range(tmp_histo.numEntries()):
                if tmp_histo.weight(i)<=0:
                    ctrl_plot_binning += 1
                    
            if ctrl_plot_binning==0 or nbins_total_bkg==20:
                total_bkg.update({"roohisto" : tmp_histo})
                break
    
    for bin_idx in range(total_bkg["roohisto"].numEntries()):
        if total_bkg["roohisto"].weight(bin_idx) > ctrl_plot_max:
            ctrl_plot_max = total_bkg["roohisto"].weight(bin_idx)

    total_bkg["roohisto"].plotOn(frame,
                                 ROOT.RooFit.Name("Total bkg"),
                                 ROOT.RooFit.Binning("plot_binning_total_bkg"),
                                 # ROOT.RooFit.Invisible(),
                                 ROOT.RooFit.LineColor(colors["total_bkg"]),
                                 ROOT.RooFit.MarkerColor(colors["total_bkg"])
                                 )
    
    pad_info.cd()
    textbox = ROOT.TPaveText(0, 0., 0.9, 0.55, "NDC NB")
    textbox.SetFillColor(0)
    textbox.SetTextSize(0.04)
    textbox.SetTextAlign(12)
    
    legend = ROOT.TLegend(0.1, 0.65, 0.9, 0.95)
    legend.SetNColumns(2)
    legend.SetFillColor(0)
    legend.SetTextSize(0.04)
    legend.SetTextAlign(12)
    legend.SetBorderSize(0)
    if imported_data:
        legend.AddEntry("Data", "Data", "lep")
        legend_obj = legend.GetListOfPrimitives().Last()
        legend_obj.SetLineColor(ROOT.kBlack)
        legend_obj.SetLineWidth(3)
    if imported_mc_sig:
        legend.AddEntry("MC signal", "MC signal", "l")
        legend_obj = legend.GetListOfPrimitives().Last()
        legend_obj.SetLineColor(colors["signal"])
        legend_obj.SetLineWidth(3)
    if imported_pdf_bkg:
        legend.AddEntry("Fitted bkg", "Fitted bkg", "l")
        legend_obj = legend.GetListOfPrimitives().Last()
        legend_obj.SetLineColor(colors["pdf_bkg_fit"])
        legend_obj.SetLineWidth(3)

    legend.AddEntry("Total bkg", "Total bkg", "lp")
    legend_obj = legend.GetListOfPrimitives().Last()
    legend_obj.SetMarkerStyle(ROOT.kFullCircle)
    legend_obj.SetMarkerColor(colors["total_bkg"])
    legend_obj.SetMarkerSize(2)
    legend_obj.SetLineColor(colors["total_bkg"])
    legend_obj.SetLineWidth(2)

    sigma_histo_bkg = sumw2_error(total_bkg['roohisto'])
    if imported_data:
        textbox.AddText(f"Data entries: {datahist.sumEntries()}")
    if imported_mc_sig:
        sigma_histo_signal = sumw2_error(mc_sig_dict['roohisto'])
        textbox.AddText(
            f"MC signal entries: {mc_sig_dict['integral']:.2f} #pm {sigma_histo_signal:.2f}") 
    if imported_pdf_bkg:
        textbox.AddText(
            f"Nbkg from fit: {pdf_bkg_obj['norm'].getVal():.2f} #pm {pdf_bkg_obj['norm'].getError():.2f}")

    textbox.AddText(f"Nbkg from MC: {total_bkg['integral']:.2f} #pm {sigma_histo_bkg:.2f}")
        
    if imported_pdf_bkg:
        delta_bkg = total_bkg['integral'] - pdf_bkg_obj['norm'].getVal()
        sigma_on_delta_bkg = (sigma_histo_bkg**2 + pdf_bkg_obj['norm'].getError()**2)**0.5
        nsigma = delta_bkg/sigma_on_delta_bkg
        textbox.AddText(f"Distance in #sigma = {nsigma:.2f}")
    
    textbox.AddText("--------------------")
    textbox.AddText("Num bkg:")
    c.Update()

    
    for bkg_key in plot_objects:
        
        bkg_obj = plot_objects[bkg_key]
        
        for nbins_bkg in [60, 30, 20, 15, 10]:
            ctrl_plot_binning = 0
            axis.setBins(nbins_bkg, f"plot_binning_{bkg_key}")
            tmp_histo = ROOT.RooDataHist(bkg_obj["roohisto"].GetName(), bkg_obj["roohisto"].GetTitle(), 
                                            ROOT.RooArgSet(axis), f"plot_binning_{bkg_key}")

            tmp_histo.add(bkg_obj["roohisto"])
            for i in range(tmp_histo.numEntries()):
                tmp_histo.get(i)
                if tmp_histo.weight(i) <= 0.01:
                    ctrl_plot_binning += 1

            if ctrl_plot_binning==0 or nbins_bkg==10:
                axis.setBins(nbins_bkg)
                tmp_histpdf = ROOT.RooHistPdf(f"{bkg_key}_pdf", f"{bkg_key}_pdf", 
                                              ROOT.RooArgSet(axis), tmp_histo)
                bkg_obj.update({"histo_pdf" : tmp_histpdf})
                break

        pad_plot.cd()
        bkg_obj["histo_pdf"].plotOn(frame, 
                         ROOT.RooFit.Name(bkg_key),
                         ROOT.RooFit.LineColor(colors[bkg_key]),
                         ROOT.RooFit.LineStyle(1),
                         ROOT.RooFit.Normalization(bkg_obj["integral"], ROOT.RooAbsReal.NumEvent))

        pad_info.cd()
        bkg_error = sumw2_error(bkg_obj["roohisto"])
        bkg_key_plot = copy(bkg_key).replace("_bkg", "")
        textbox.AddText(
            f"  {bkg_key_plot} = {bkg_obj['integral']:.2f} #pm {bkg_error:.2f}")
        legend.AddEntry(bkg_key_plot, bkg_key_plot, "l")
        legend_obj = legend.GetListOfPrimitives().Last()
        legend_obj.SetLineColor(colors[bkg_key])
        legend_obj.SetLineWidth(3)

    
    legend.Draw()

    pad_plot.cd()
    pad_plot.SetLogy() if logscale is True else pad_plot.SetLogy(False)
    
    frame.SetMaximum(1.5*ctrl_plot_max)      
    frame.SetMinimum(5e-2)
    frame.Draw()
    pad_plot.Update()
 
 
    CMS_lumi(pad_plot, 5, 0, simulation=True)

    '''
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
    '''
    
    pad_info.cd()
    textbox.Draw()  
    c.Update()

    c.SaveAs(f"{figpath}/bkg_{bin_key}_{flag}.pdf")

###############################################################################

def plot_projected_bkg(plot_dictionary, binning_pt, binning_eta, flag, logscale=True, figpath=''):
    """
    """
    bins_pt, bins_eta = binning(binning_pt), binning(binning_eta)

    if not (len(bins_pt) ==2 or len(bins_eta) == 2):
        print("ERROR: the projection has to be made by collapsing one of the two dimensions")
        sys.exit()

    if len(bins_pt) > len(bins_eta):
        proj_axis = "pt"
        proj_binning = bins_pt
    else:
        proj_axis = "eta"
        proj_binning = bins_eta

    c = ROOT.TCanvas(f"bkg_{proj_axis}_{flag}", "c", 900, 900)
    c.cd()

    style_settings()

    pad_title = ROOT.TPad("pad_title", "pad_title", 0, 0.9, 1, 1)
    pad_plot = ROOT.TPad("pad_plot", "pad_plot", 0, 0, 1, 0.9)

    pad_title.SetMargin(0.1, 0.1, 0.1, 0.1), pad_title.Draw()
    pad_plot.SetMargin(0.12, 0.05, 0.12, 0.05), pad_plot.Draw()

    pad_title.cd()
    titlebox = ROOT.TPaveText(0, 0, 1, 1, "NDC NB")
    titlebox.SetFillColor(0)
    # titlebox.SetTextFont(42)
    titlebox.SetTextSize(0.35)
    titlebox.AddText(f"Bkg {proj_axis} distrib. - {flag}ing probes")
    titlebox.Draw()
    c.Update()

    plot_objects = copy(plot_dictionary)

    # ref_histo = plot_objects.pop("ref_histo")

    pad_plot.cd()
    axis = ROOT.RooRealVar(proj_axis, proj_axis, proj_binning[0], proj_binning[-1])
    axis.setBinning(ROOT.RooBinning(len(proj_binning)-1, proj_binning))

    tot_bkg_roohist = ROOT.RooDataHist("tot_bkg_roohist", "tot_bkg_roohist", ROOT.RooArgSet(axis), "")

    frame = axis.frame(ROOT.RooFit.Bins(axis.getBins()))
    frame.SetTitle("")
    frame.SetTitleSize(0)
    frame_yaxis = frame.GetYaxis()
    if proj_axis == "pt":
        frame_yaxis.SetTitle("Events / (2 GeV)")
    elif proj_axis == "eta":
        frame_yaxis.SetTitle("Events / (0.1)")
    frame_yaxis.SetTitleSize(0.035)
    frame_yaxis.SetTitleOffset(1.2)
    frame_axis = frame.GetXaxis()
    # xtitle = frame_axis.GetTitle(), frame_axis.SetTitleSize(0)
    if proj_axis == "pt":
        frame_axis.SetTitle("p_{T}^{#mu} [GeV]")
    elif proj_axis == "eta":
        frame_axis.SetTitle("#eta^{#mu}")
    frame_axis.SetTitleSize(0.035)
    frame_axis.SetTitleOffset(1.2)


    for bkg_key in plot_objects.keys():

        bkg_hist = plot_objects[bkg_key]

        if proj_axis == "pt":
            nbins_eta = bkg_hist.GetNbinsY()
            bkg_hist = bkg_hist.ProjectionX(f"{bkg_hist.GetName()}_{proj_axis}", 1, nbins_eta, "e")
        elif proj_axis == "eta":
            nbins_pt = bkg_hist.GetNbinsX()
            bkg_hist = bkg_hist.ProjectionY(f"{bkg_hist.GetName()}_{proj_axis}", 1, nbins_pt, "e")
        
        roohisto_bkg = ROOT.RooDataHist(f"{bkg_key}_roohist", "", ROOT.RooArgSet(axis), bkg_hist)

        roohistpdf = ROOT.RooHistPdf(f"{bkg_key}_pdf", "", ROOT.RooArgSet(axis), roohisto_bkg)

        roohistpdf.plotOn(frame,
                          ROOT.RooFit.Name(bkg_key),
                          ROOT.RooFit.LineColor(colors[bkg_key]),
                          ROOT.RooFit.LineStyle(1),
                          ROOT.RooFit.Normalization(roohisto_bkg.sumEntries(), ROOT.RooAbsReal.NumEvent))
        
        tot_bkg_roohist.add(roohisto_bkg)

    tot_bkg_roohist.plotOn(frame,
                           ROOT.RooFit.Name("Total bkg"),
                           ROOT.RooFit.Binning("plot_binning_total_bkg"),
                           # ROOT.RooFit.Invisible(),
                           ROOT.RooFit.LineColor(colors["total_bkg"]),
                           ROOT.RooFit.MarkerColor(colors["total_bkg"]))
    
    pad_plot.SetLogy() if logscale is True else pad_plot.SetLogy(False)   

    ctrl_plot_max=0
    for bin_idx in range(tot_bkg_roohist.numEntries()):
            if tot_bkg_roohist.weight(bin_idx) > ctrl_plot_max:
                ctrl_plot_max = tot_bkg_roohist.weight(bin_idx)
    frame.SetMaximum(1.5*ctrl_plot_max)
    frame.SetMinimum(1)
    frame.Draw()
    pad_plot.Update()
 
    CMS_lumi(pad_plot, 5, 0, simulation=True) 

    c.Update()

    c.SaveAs(f"{figpath}/bkg_{proj_axis}_{flag}.pdf")


    
        





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
    plot_minv_distrib_with_fit(pad_plot, "pass",axis, data, gaus, pull=False)
    
    c.SaveAs("../prova_plot.pdf")

