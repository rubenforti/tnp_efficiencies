"""
"""
import ROOT
import sys
import os
import pickle
from copy import deepcopy
from array import array
from utilities.base_lib import binnings, sumw2_error
from utilities.CMS_lumi import CMS_lumi



colors = { 
    "bkg_QCD" : ROOT.kOrange+1,
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
    "bkg_total" : ROOT.kRed,
    "bkg_total_SS" : ROOT.kRed-2, #ROOT.kBlack, 
    "mc" : ROOT.kBlue, #ROOT.kRed,
    "mc_SS" : ROOT.kBlue-2, #ROOT.kOrange+4, 
    "pdf_bkg_fit" : ROOT.kRed-2,

    "QCD" : ROOT.kOrange+1,
    "Diboson" : ROOT.kGreen,
    "Top" : ROOT.kCyan+3,
    "Ztautau" : ROOT.kMagenta+1,
    "Wjets" : ROOT.kOrange+4,
    "Zjets" : ROOT.kBlue+4,
    "Diboson_SS" : ROOT.kGreen,
    "Top_SS" : ROOT.kCyan+3,
    "Ztautau_SS" : ROOT.kMagenta+1,
    "Wjets_SS" : ROOT.kOrange+4,
    "Zjets_SS" : ROOT.kBlue+4,
    "SameCharge" : ROOT.kYellow+2,
}


bkg_grouping = {
    "QCD" : ["bkg_QCD"],

    "SameCharge" : ["bkg_SameCharge"],

    "Diboson" : ["bkg_WW", "bkg_WZ", "bkg_ZZ"],
    "Top" : ["bkg_TTSemileptonic", "bkg_TTFullyleptonic"],
    "Ztautau" : ["bkg_Ztautau"],
    "Wjets" : ["bkg_WplusJets", "bkg_WminusJets"],
    "Zjets" : ["bkg_Zjets"],

    "Diboson_SS" : ["bkg_WW_SS", "bkg_WZ_SS", "bkg_ZZ_SS"],
    "Top_SS" : ["bkg_TTSemileptonic_SS", "bkg_TTFullyleptonic_SS"],
    "Ztautau_SS" : ["bkg_Ztautau_SS"],
    "Wjets_SS" : ["bkg_WplusJets_SS", "bkg_WminusJets_SS"],
    "Zjets_SS" : ["bkg_Zjets_SS"], 
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
                   #ROOT.RooFit.NormRange("x_binning"),
                   #ROOT.RooFit.Range("x_binning"),
                   ROOT.RooFit.LineColor(ROOT.kRed))
    
    isBB_bkg = False
    for serv in pdf_fit.servers():
        if "BB" in serv.GetName():
            isBB_bkg = True
            break
       
    for comp in pdf_fit.getComponents():
        if (isBB_bkg is False and "bkg" in comp.GetName()) or \
           (isBB_bkg is True  and "BB_bkg" in comp.GetName() and "_paramPdf" not in comp.GetName()):
            bkg_set = ROOT.RooArgSet(comp)
            pdf_fit.plotOn(plot_frame,
                # ROOT.RooFit.NormRange("x_binning"),
                ROOT.RooFit.Name("Bkg pdf"),
                ROOT.RooFit.Components(bkg_set),
                ROOT.RooFit.LineColor(ROOT.kBlue),
                ROOT.RooFit.LineStyle(ROOT.kDashed),
                ROOT.RooFit.LineWidth(4))   

    plot_frame.Draw()

    CMS_lumi(pad, 5, 0, simulation=False)
    pad.Update()
    
###############################################################################


def plot_fitted_pass_fail(type_analysis, plot_objects, bin_key, figpath=""): 
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

    pad_info.SetMargin(0.15, 0.05, 0.05, 0.05)
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

    binning_pass = pass_obj["axis"].getBinning("x_binning")
    binning_fail = fail_obj["axis"].getBinning("x_binning")
    
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
        try:
            chi2 = float(pass_obj["res"].GetTitle().split("_")[0])
            ndof = int(pass_obj["res"].GetTitle().split("_")[1])
            stats_pass.AddText(f"Chi2/ndof = {chi2:.1f} / {ndof}")
        except:
            print("Chi2/ndof not available")
        stats_pass.Draw()
        c.Update()

        [stats_fail.AddText(f"{par.GetTitle()} = {par.getVal():.3f} #pm {par.getError():.3f}") 
         for par in fail_obj["res"].floatParsFinal()]
        stats_fail.AddText("------------------------------")
        # stats_box.AddText(f"Chi2/ndof = {res.chi2(), 2) / pars.getSize()}")
        stats_fail.AddText(f"Status = {fail_obj['res'].status()}")
        stats_fail.AddText(f"Cov quality = {fail_obj['res'].covQual()}")
        stats_fail.AddText(f"Edm = {fail_obj['res'].edm():.5f}")
        try:
            chi2 = float(fail_obj["res"].GetTitle().split("_")[0])
            ndof = int(fail_obj["res"].GetTitle().split("_")[1])
            stats_fail.AddText(f"Chi2/ndof = {chi2:.1f} / {ndof}")
        except:
            print("Chi2/ndof not available")
        stats_fail.Draw()
        c.Update()

    elif "sim" in type_analysis:
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
        try:
            chi2 = float(pass_obj["res"].GetTitle().split("_")[0])
            ndof = int(pass_obj["res"].GetTitle().split("_")[1])
            stats_gen.AddText(f"Chi2/ndof = {chi2:.1f} / {ndof}")
        except:
            print("Chi2/ndof not available")
        stats_pass.Draw(), stats_fail.Draw(), stats_gen.Draw()
        c.Update()

    else:
        sys.exit("ERROR: type_analysis not recognized")

    c.SaveAs(f"{figpath}/fit_pf_{type_analysis}_{bin_key}.pdf")

###############################################################################


def plot_bkg_object(frame, axis, hist_list, label, list_nbins_plot):
    """
    """
    isHistPlot = True if "SameCharge" in label else False
    label = deepcopy(label)#.replace("_SS", "")
    minNBins = 20 if isHistPlot else 15
    
    for nbins_total_bkg in list_nbins_plot:

        ctrl_plot_binning = 0
        axis.setBins(nbins_total_bkg, f"plot_binning_{label}")
        tmp_histo = ROOT.RooDataHist(f"{label}_histo", f"{label}_histo", 
                                     ROOT.RooArgSet(axis), f"plot_binning_{label}")
        
        for hist in hist_list: tmp_histo.add(hist)

        for i in range(tmp_histo.numEntries()):
            if tmp_histo.weight(i) < 0: ctrl_plot_binning += 1
                
        if ctrl_plot_binning==0 or nbins_total_bkg==minNBins: break

    plot_commands = ROOT.RooLinkedList()
    plot_commands.Add(ROOT.RooFit.Name(label))
    plot_commands.Add(ROOT.RooFit.LineColor(colors[label.replace("_SS", "")]))


    if isHistPlot:
        bkg_obj_plot = tmp_histo
        plot_commands.Add(ROOT.RooFit.Binning(f"plot_binning_{label}"))
        plot_commands.Add(ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))
        plot_commands.Add(ROOT.RooFit.MarkerColor(colors[label.replace("_SS", "")]))
    else:
        bkg_obj_plot = ROOT.RooHistPdf(f"{label}_pdf", f"{label}_pdf", ROOT.RooArgSet(axis), tmp_histo)
        plot_commands.Add(ROOT.RooFit.LineColor(colors[label.replace("_SS", "")]))
        plot_commands.Add(ROOT.RooFit.LineStyle(1))
        plot_commands.Add(ROOT.RooFit.Normalization(tmp_histo.sumEntries(), ROOT.RooAbsReal.NumEvent))

    bkg_obj_plot.plotOn(frame, plot_commands)

    # print(f"{label} plotted with {nbins_total_bkg} bins")

###############################################################################


def plot_bkg(plot_dictionary, flag, bin_key,
             charge_separation="",
             group_backgrounds=True, logscale=False, figpath=''):
    """
    Function that creates and saves the plot containing the mass distributions 
    of the various bkg processes and the total bkg; a reference histogram 
    (data or signal MC) can be plotted in the same figure. For graphical 
    choice, the bkg distributions of single processes (and the signal MC 
    distribution) are plotted as RooHistPdf.
    The objects to be plotted are contained in the dictionary "plot_objects", 
    which must be structured as follows:
        plot_objects = { "axis" : RooRealVar,
                         "bkg_total" : RooDataHist,
                         "bkg_{process}" : RooDataHist, 
                         "mc (optional)" : RooDataHist,
                         "data (optional)" : RooDataHist,
                         "fit_pars (optional)" : {"nsig":RooRealVar, "nbkg":RooRealVar},
                         "pdf_bkg_fit (optional)" : {"pdf":RooAbsPdf, "norm":RooRealVar}
                        }
    """

    imported_data, imported_mc_sig, imported_pdf_bkg = False, False, False

    c = ROOT.TCanvas(f"c_{bin_key}_{flag}", "c", 1600, 900)
    c.cd()

    style_settings()

    pad_title = ROOT.TPad("pad_title", "pad_title", 0, 0.9, 1, 1)
    pad_plot = ROOT.TPad("pad_plot", "pad_plot", 0, 0, 0.6, 0.9)
    # pad_pull = ROOT.TPad("pad_pull", "pad_pull", 0, 0, 0.65, 0.25)
    pad_info = ROOT.TPad("pad_info", "pad_info", 0.6, 0.05, 1, 0.9)

    pad_title.SetMargin(0.1, 0.1, 0.1, 0.1), pad_title.Draw()
    pad_plot.SetMargin(0.12, 0.05, 0.12, 0.05), pad_plot.Draw()
    pad_info.SetMargin(0.05, 0.05, 0.1, 0.9), pad_info.Draw()

    pad_title.cd()
    titlebox = ROOT.TPaveText(0, 0, 1, 1, "NDC NB")
    titlebox.SetFillColor(0)
    # titlebox.SetTextFont(42)
    titlebox.SetTextSize(0.35)

    title_str = f"Background {flag}ing probes"
    if charge_separation != "": 
        title_str += f" {charge_separation}"
    if bin_key != "[24.0to65.0][-2.4to2.4]": 
        title_str += f" - bin {bin_key}"

    titlebox.AddText(title_str)
    titlebox.Draw()
    c.Update()

    plot_objects = deepcopy(plot_dictionary)

    axis = plot_objects["axis"]

    binning_plot = axis.getBinning("plot_binning")
    nbins_plot_max = binning_plot.numBins()

    if nbins_plot_max == 80: 
        list_nbins_plot = [80, 40, 20, 16, 10]
    elif nbins_plot_max == 60:
        list_nbins_plot = [60, 30, 20, 15, 10]
    else:
        list_nbins_plot = [nbins_plot_max]

    plotted_histos = []


    # Initializing the plotting pad
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

    # Initializing the info pad
    pad_info.cd()
    # Create the legend in the upper part of the pad
    legend = ROOT.TLegend(0.1, 0.65, 0.9, 0.95)
    legend.SetNColumns(2)
    legend.SetFillColor(0)
    legend.SetTextSize(0.04)
    legend.SetTextAlign(12)
    legend.SetBorderSize(0)
    # Create the textbox in the lower part of the pad
    textbox = ROOT.TPaveText(0, 0., 0.9, 0.6, "NDC NB")
    textbox.SetFillColor(0)
    textbox.SetTextSize(0.04)
    textbox.SetTextAlign(12)
        
    # Plotting the total bkg

    postfix = "" if charge_separation == "" else f"_{charge_separation}"
    ### Temporary, just for consistency
    if postfix == "_OS": postfix = ""
    ###
    print(charge_separation, postfix)
    total_bkg = plot_objects[f"bkg_total{postfix}"]
    pad_plot.cd()
    plot_bkg_object(frame, axis, [total_bkg], "bkg_total", list_nbins_plot)
    plotted_histos.append(total_bkg)
    
    pad_info.cd()
    legend.AddEntry(f"Total bkg {charge_separation}", f"Total bkg {charge_separation}", "l")
    legend_obj = legend.GetListOfPrimitives().Last()
    legend_obj.SetLineColor(colors["bkg_total"])
    legend_obj.SetLineWidth(3)

    totbkg_label = "Nbkg total from MC" if charge_separation == "" else f"Nbkg {charge_separation} total from MC"
    textbox.AddText(f"{totbkg_label}: {total_bkg.sumEntries():.0f} #pm {sumw2_error(total_bkg):.0f}")


    # Plotting the single bkg processes
    pad_info.cd()
    textbox.AddText("--------------------")
    textbox.AddText("Background events for category:")
        
    
    bkg_plot_it_dict = {k : v for k , v in plot_objects.items() if "bkg_" in k}

    if group_backgrounds:
        bkg_obj_tmp_container = deepcopy(bkg_plot_it_dict)
        bkg_grouping_tmp = {}
        for bkg_group, bkg_cat_list in bkg_grouping.items():
            bkg_grouping_tmp[bkg_group] = []
            for bkg_cat in bkg_cat_list:
                if bkg_cat in bkg_obj_tmp_container.keys(): bkg_grouping_tmp[bkg_group].append(bkg_cat)
        bkg_plot_it_dict = bkg_grouping_tmp

    for bkg_cat, bkg_obj in bkg_plot_it_dict.items():

        if bkg_cat == "bkg_total": continue
    
        bkg_key_print = deepcopy(bkg_cat).replace("bkg_", "")

        nEntries = 0
        if group_backgrounds:
            hist_plotting = [plot_objects[key] for key in bkg_obj if type(plot_objects[key]) is ROOT.RooDataHist]
            not_found_bkgs = [key for key in bkg_obj if type(plot_objects[key]) is not ROOT.RooDataHist]
            for nf_bkg in not_found_bkgs:
                print(f"WARNING: {nf_bkg} not found as RooDataHist object in the workspace. Skipping...")
            if len(hist_plotting) == 0: 
                continue
            for hist in hist_plotting:
                nEntries += hist.sumEntries()
            
        else:
            hist_plotting = [bkg_obj]
            nEntries = bkg_obj.sumEntries()       
        
        if nEntries < 0: continue

        '''
        if "SameCharge" not in bkg_cat:
            #if (bkg_cat == list(bkg_plot_it_dict.keys())[-1]) or (bkg_cat=="Zjets"): continue
            #if bkg_cat == list(bkg_plot_it_dict.keys())[-2]: continue
            textbox.AddText(f"\n")
            legend.AddEntry(bkg_cat, " ", "l")
            legend_obj = legend.GetListOfPrimitives().Last()
            legend_obj.SetLineColor(ROOT.kWhite)
            legend_obj.SetLineWidth(3)
            continue
        '''

        
        pad_plot.cd()
        plot_bkg_object(frame, axis, hist_plotting, bkg_cat, list_nbins_plot)
        bkg_obj_new = hist_plotting[0]
        for hist in hist_plotting[1:]: bkg_obj_new.add(hist)
        plotted_histos.append(bkg_obj_new)
        pad_info.cd()
        bkg_error = sumw2_error(bkg_obj_new)
        textbox.AddText(f"  {bkg_key_print} = {bkg_obj_new.sumEntries():.2f} #pm {bkg_error:.2f}")
        leg_marker_opt = "l" if "SameCharge" not in bkg_cat else "lp"
        legend.AddEntry(bkg_key_print, bkg_key_print, leg_marker_opt)
        legend_obj = legend.GetListOfPrimitives().Last()
        legend_obj.SetLineColor(colors[bkg_cat.replace("_SS", "")])
        legend_obj.SetLineWidth(3)
        if "SameCharge" in bkg_cat: 
            legend_obj.SetMarkerStyle(ROOT.kFullCircle)
            legend_obj.SetMarkerColor(colors["bkg_SameCharge"])
            legend_obj.SetMarkerSize(2)

    # Plotting other objects
    pad_info.cd()
    textbox.AddText("--------------------")

    if "data" in plot_objects.keys():
        datahist = plot_objects["data"]
        pad_plot.cd()
        datahist.plotOn(frame, 
                        ROOT.RooFit.Binning("plot_binning"),
                        ROOT.RooFit.Name("Data"),
                        ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
        plotted_histos.append(datahist)
        pad_info.cd()
        legend.AddEntry("Data", "Data", "lep")
        legend_obj = legend.GetListOfPrimitives().Last()
        legend_obj.SetLineColor(ROOT.kBlack)
        legend_obj.SetLineWidth(3)
        textbox.AddText(f"Data events: {datahist.sumEntries():.0f} #pm {datahist.sumEntries()**0.5:.0f}")
    
    if "mc" in plot_objects.keys() and not (flag == "pass" and logscale is False):
        mc_sig_hist = plot_objects["mc"]
        pad_plot.cd()
        mc_sig_roohistpdf = ROOT.RooHistPdf(f"MC signal", f"MC signal", ROOT.RooArgSet(axis), mc_sig_hist)
        mc_sig_roohistpdf.plotOn(frame,
                                 ROOT.RooFit.Name("MC signal"),
                                 ROOT.RooFit.LineColor(colors["mc"]),
                                 ROOT.RooFit.Normalization(mc_sig_hist.sumEntries(), ROOT.RooAbsReal.NumEvent))
        plotted_histos.append(mc_sig_hist)
        pad_info.cd()
        legend.AddEntry("MC signal", "MC signal", "l")
        legend_obj = legend.GetListOfPrimitives().Last()
        legend_obj.SetLineColor(colors["mc"])
        legend_obj.SetLineWidth(3)
        textbox.AddText(
            f"MC signal entries: {mc_sig_hist.sumEntries():.2f} #pm {sumw2_error(mc_sig_hist):.2f}")
    
    if "mc_SS" in plot_objects.keys():
        mc_SS_hist = plot_objects["mc_SS"]
        pad_plot.cd()
        mc_SS_roohistpdf = ROOT.RooHistPdf(f"MC_SS", f"MC_SS", ROOT.RooArgSet(axis), mc_SS_hist)
        mc_SS_roohistpdf.plotOn(frame,
                                 ROOT.RooFit.Name("MC_SS"),
                                 ROOT.RooFit.LineColor(colors["mc_SS"]),
                                 ROOT.RooFit.Normalization(mc_SS_hist.sumEntries(), ROOT.RooAbsReal.NumEvent))
        plotted_histos.append(mc_SS_hist)
        pad_info.cd()
        legend.AddEntry("MC_SS", "MC_SS", "l")
        legend_obj = legend.GetListOfPrimitives().Last()
        legend_obj.SetLineColor(colors["mc_SS"])
        legend_obj.SetLineWidth(3)
        textbox.AddText(
            f"MC SS entries: {mc_SS_hist.sumEntries():.2f} #pm {sumw2_error(mc_SS_hist):.2f}") 

    if "pdf_bkg_fit" in plot_objects.keys():
        pad_plot.cd()
        pdf_bkg_obj = plot_objects.pop("pdf_bkg_fit")
        pdf_bkg_obj["pdf"].plotOn(frame, 
                                  ROOT.RooFit.Name("Fitted bkg"),
                                  ROOT.RooFit.LineColor(colors["pdf_bkg_fit"]),
                                  ROOT.RooFit.Normalization(pdf_bkg_obj["norm"].getVal(), ROOT.RooAbsReal.NumEvent))
        pad_info.cd()
        legend.AddEntry("Fitted bkg", "Fitted bkg", "l")
        legend_obj = legend.GetListOfPrimitives().Last()
        legend_obj.SetLineColor(colors["pdf_bkg_fit"])
        legend_obj.SetLineWidth(3)
        delta_bkg = total_bkg.sumEntries() - pdf_bkg_obj['norm'].getVal()
        sigma_on_delta_bkg = (sumw2_error(total_bkg)**2 + pdf_bkg_obj['norm'].getError()**2)**0.5
        nsigma = delta_bkg/sigma_on_delta_bkg
        textbox.AddText(f"Distance in #sigma = {nsigma:.2f}")
        textbox.AddText(
            f"Nbkg from fit: {pdf_bkg_obj['norm'].getVal():.2f} #pm {pdf_bkg_obj['norm'].getError():.2f}")

    
    # Setting the control variable for the maximum value of the y-axis
    ctrl_plot_max = 0  
    for plot_hist in plotted_histos:
        #print(plot_hist.GetName())
        for bin_idx in range(plot_hist.numEntries()):
            #print(plot_hist.weight(bin_idx), ctrl_plot_max)
            if plot_hist.weight(bin_idx) > ctrl_plot_max:
                ctrl_plot_max = plot_hist.weight(bin_idx)


    # c.Update()

    pad_plot.cd()
    if logscale is True:
        pad_plot.SetLogy()
        frame.SetMaximum(2.5*ctrl_plot_max)
    
    frame.SetMinimum(0.1)
    
    # print(ctrl_plot_max, frame.GetMaximum())
    frame.Draw()
    #pad_plot.Update()

 
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
    legend.Draw()
    textbox.Draw()  
    c.Update()

    filename = f"bkg_{bin_key}_{flag}" if bin_key!="[24.0to65.0][-2.4to2.4]" else f"bkg_total_{flag}"
    #if charge_separation != "": filename += f"_{charge_separation}"

    c.SaveAs(f"{figpath}/{flag}/{filename}.png")

###############################################################################


def plot_bkg_comparison(plot_dictionary, flag, bin_key, figpath=''):
    """
    Function that creates and saves the plot containing the mass distributions 
    of two different backgrounds. The objects to be plotted are contained in 
    the dictionary "plot_objects", which must be structured as follows:
        plot_objects = { 
            "axis":RooRealVar, "bkg_1":RooDataHist, "bkg_2":RooDataHist 
            }
    The order of the dictionary has to be the one shown here; the bkg keys are
    the labels used in the plots to refer to the two histograms.
    """

    plot_objects = deepcopy(plot_dictionary)

    hist_list = list(plot_objects.values())
    hist_names = list(plot_objects.keys())

    c = ROOT.TCanvas(f"c_{bin_key}_{flag}", "c", 900, 900)
    c.cd()

    style_settings()

    pad_title = ROOT.TPad("pad_title", "pad_title", 0, 0.9, 1, 1)
    pad_plot = ROOT.TPad("pad_plot", "pad_plot", 0, 0.3, 1, 0.9)
    pad_ratio = ROOT.TPad("pad_ratio", "pad_ratio", 0, 0, 1, 0.3)

    pad_title.SetMargin(0.1, 0.1, 0.1, 0.1), pad_title.Draw()
    pad_plot.SetMargin(0.1, 0.05, 0, 0.05), pad_plot.Draw()
    pad_ratio.SetMargin(0.1, 0.05, 0.2, 0.0), pad_ratio.Draw()

    pad_title.cd()
    titlebox = ROOT.TPaveText(0, 0, 1, 1, "NDC NB")
    titlebox.SetFillColor(0)
    # titlebox.SetTextFont(42)
    titlebox.SetTextSize(0.25)

    title_str = f"{hist_names[1].replace('bkg_','')}-{hist_names[2].replace('bkg_','')} comparison ({flag}ing probes)"
    if bin_key != "[24.0to65.0][-2.4to2.4]": 
        title_str += f" - bin {bin_key}"

    titlebox.AddText(title_str)
    titlebox.Draw()
    c.Update()


    axis = plot_objects["axis"]

    binning_plot = axis.getBinning("plot_binning")
    nbins = binning_plot.numBins()


    # Initializing the plotting pad
    pad_plot.cd()
    frame = axis.frame(ROOT.RooFit.Bins(axis.getBins("plot_binning")))
    frame.SetTitle("")
    frame.SetTitleSize(0)
    frame_yaxis = frame.GetYaxis()
    frame_yaxis.SetTitle("Events / (1 GeV)")
    frame_yaxis.SetTitleSize(0.032)
    frame_yaxis.SetTitleOffset(1.2)
    frame_axis = frame.GetXaxis()
    frame_axis.SetTitle("")
    frame_axis.SetTitleSize(0)

    # Create the legend in the upper part of the plot pad
    pad_plot.cd()
    legend = ROOT.TLegend(0.65, 0.8, 0.95, 0.95)
    legend.SetNColumns(1)
    legend.SetFillColor(0)
    legend.SetTextSize(0.03)
    legend.SetTextAlign(12)
    legend.SetBorderSize(1)


    #Initializing the ratio pad
    pad_ratio.cd()    
    frame_ratio = axis.frame(ROOT.RooFit.Bins(nbins))
    frame_ratio.SetTitle("")
    frame_ratio.SetTitleSize(0)
    frame_ratio_yaxis = frame_ratio.GetYaxis()

    hist_names_label = [h.replace("bkg_","")+"_MC" if "total" in h else h.replace("bkg_","") for h in hist_names]

    frame_ratio_yaxis.SetTitle(f"Ratio {hist_names_label[1]} / {hist_names_label[2]}")
    frame_ratio_yaxis.SetTitleSize(0.064)
    frame_ratio_yaxis.SetTitleOffset(0.6)
    frame_ratio_yaxis.SetLabelSize(0.07)
    frame_ratio_axis = frame_ratio.GetXaxis()
    frame_ratio_axis.SetTitle("M_{inv} TP [GeV]")
    frame_ratio_axis.SetTitleSize(0.07)
    frame_ratio_axis.SetTitleOffset(1.2)
    frame_ratio_axis.SetLabelSize(0.07)

    for bkg_key, bkg_hist in plot_objects.items():
        if bkg_key=="axis": continue

        plot_commands = ROOT.RooLinkedList()
        plot_commands.Add(ROOT.RooFit.Name(bkg_key))
        plot_commands.Add(ROOT.RooFit.LineColor(colors[bkg_key.replace("_SS", "")]))
        
        if "total" in bkg_key:
            bkg_obj_plot = ROOT.RooHistPdf(f"{bkg_key}_pdf", f"{bkg_key}_pdf", ROOT.RooArgSet(axis), bkg_hist)
            plot_commands.Add(ROOT.RooFit.LineColor(colors[bkg_key]))
            plot_commands.Add(ROOT.RooFit.LineStyle(1))
            plot_commands.Add(ROOT.RooFit.Normalization(bkg_hist.sumEntries(), ROOT.RooAbsReal.NumEvent))
        else:
            bkg_obj_plot = bkg_hist
            plot_commands.Add(ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson))
            plot_commands.Add(ROOT.RooFit.MarkerColor(colors[bkg_key.replace("_SS", "")]))

        pad_plot.cd()
        bkg_obj_plot.plotOn(frame, plot_commands)


        legend.AddEntry(bkg_key, f"{bkg_key.replace('bkg_','')}: {bkg_hist.sumEntries():.0f} #pm {sumw2_error(bkg_hist):.0f}", "l" if "MC" in bkg_key else "lp")
        legend_obj = legend.GetListOfPrimitives().Last()
        legend_obj.SetLineColor(colors[bkg_key.replace("_SS", "")])
        legend_obj.SetMarkerColor(colors[bkg_key.replace("_SS", "")])
        legend_obj.SetLineWidth(3)


    print(hist_list[1].sumEntries(), hist_list[2].sumEntries())

    total_ratio = hist_list[1].sumEntries()/hist_list[2].sumEntries()

    legend.AddEntry(" ", f"Ratio = {total_ratio:.2f}", "")

    pad_plot.cd()
    frame.Draw()
    legend.Draw()
    CMS_lumi(pad_plot, 5, 0, simulation=(True if ("MC" in hist_names[1] and "MC" in hist_names[2]) else False))


    array_ratio, array_ratio_err = [1]*nbins, [0]*nbins
    max_ratio, min_ratio = 1, 1
    for i in range(nbins):
        if hist_list[2].weight(i) == 0: continue
        array_ratio[i] = hist_list[1].weight(i)/hist_list[2].weight(i)
        if array_ratio[i] == 0: continue
        array_ratio_err[i] = array_ratio[i] * (
            (hist_list[1].weightError(ROOT.RooAbsData.SumW2)/hist_list[1].weight(i))**2 + 
            (hist_list[2].weightError(ROOT.RooAbsData.SumW2)/hist_list[2].weight(i))**2 )**0.5
        if array_ratio[i] > max_ratio: max_ratio = array_ratio[i]
        if array_ratio[i] < min_ratio: min_ratio = array_ratio[i]
        
    pad_ratio.cd()
    histo_ratio = ROOT.RooDataHist("ratio", "ratio", ROOT.RooArgSet(axis), "plot_binning")
    for i in range(nbins):
        axis.setVal(binning_plot.binCenter(i))
        histo_ratio.add(ROOT.RooArgSet(axis), array_ratio[i], array_ratio_err[i]**2)
    histo_ratio.plotOn(frame_ratio, ROOT.RooFit.Binning("plot_binning"), ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))
    
    frame_ratio.SetMaximum(max(max_ratio+0.2, 1.2))
    frame_ratio.SetMinimum(min(min_ratio-0.2, 0.8))

    frame_ratio.Draw()

    hline = ROOT.TLine(axis.getBinning("plot_binning").binLow(0), 1, axis.getBinning("plot_binning").binHigh(nbins-1), 1)
    hline.SetLineStyle(2)
    hline.Draw()

    hline_avg = ROOT.TLine(axis.getBinning("plot_binning").binLow(0), total_ratio,
                           axis.getBinning("plot_binning").binHigh(nbins-1), total_ratio)
    hline_avg.SetLineColor(ROOT.kRed)
    hline_avg.Draw()

    filename = f"bkgComp_{bin_key}_{flag}" if bin_key!="[24.0to65.0][-2.4to2.4]" else f"bkgComp_total_{flag}"
    #if charge_separation != "": filename += f"_{charge_separation}"

    c.SaveAs(f"{figpath}/{flag}/{filename}.png")



###############################################################################

def plot_2d_bkg_distrib(histos_dict, bkg_cat, figpath=""):
    """
    """
    c = ROOT.TCanvas(f"c_{bkg_cat}", "c", 1600, 900)
    c.cd()

    if bkg_cat != "TTFullyleptonic_bkg":
        bkg_key_print = deepcopy(bkg_cat).replace("_bkg", "") 
    else:
        bkg_key_print = "TTLeptonic"

    hist_pass, hist_fail = histos_dict["pass"], histos_dict["fail"]

    style_settings()

    print(hist_pass.GetEntries(), hist_fail.GetEntries())

    pad_title = ROOT.TPad("pad_title", "pad_title", 0, 0.94, 1, 1)
    pad_subtitle_pass = ROOT.TPad("pad_subtitle_pass", "pad_subtitle_pass", 0, 0.9, 0.5, 0.94)
    pad_subtitle_fail = ROOT.TPad("pad_subtitle_fail", "pad_subtitle_fail", 0.5, 0.9, 1, 0.94)
    pad_pass = ROOT.TPad("pad_pass", "pad_plot", 0, 0, 0.5, 0.9)
    pad_fail = ROOT.TPad("pad_fail", "pad_plot", 0.5, 0, 1, 0.9)

    pad_title.SetMargin(0.1, 0.1, 0.1, 0.1), pad_title.Draw()
    pad_subtitle_pass.SetMargin(0.1, 0.1, 0.1, 0.1), pad_subtitle_pass.Draw()
    pad_subtitle_fail.SetMargin(0.1, 0.1, 0.1, 0.1), pad_subtitle_fail.Draw()
    pad_pass.SetMargin(0.1, 0.135, 0.12, 0.05), pad_pass.Draw()
    pad_fail.SetMargin(0.1, 0.135, 0.12, 0.05), pad_fail.Draw()


    pad_title.cd()
    titlebox = ROOT.TPaveText(0, 0, 1, 1, "NDC NB")
    titlebox.SetFillColor(0)
    # titlebox.SetTextFont(42)
    titlebox.SetTextSize(0.6)
    titlebox.AddText(f"{bkg_key_print}  2D distribution")
    titlebox.Draw()
    c.Update()

    pad_subtitle_pass.cd()
    title_pass = ROOT.TPaveText(0, 0, 1, 1, "NDC NB")
    title_pass.SetFillColor(0)
    title_pass.SetTextSize(0.6)
    title_pass.AddText(f"Passing probes")
    title_pass.Draw()
    c.Update()

    pad_subtitle_fail.cd()
    title_fail = ROOT.TPaveText(0, 0, 1, 1, "NDC NB")
    title_fail.SetFillColor(0)
    title_fail.SetTextSize(0.6)
    title_fail.AddText(f"Failing probes")
    title_fail.Draw()
    c.Update()

    pad_pass.cd()
    hist_pass.Draw("colz")
    hist_pass.SetTitle("")
    hist_pass.SetTitleSize(0)
    Yaxis_pass = hist_pass.GetYaxis()
    Yaxis_pass.SetTitle("#eta^{#mu}")
    Yaxis_pass.SetTitleSize(0.035)
    Yaxis_pass.SetTitleOffset(1.2)
    Xaxis_pass = hist_pass.GetXaxis()
    Xaxis_pass.SetTitle("p_{T}^{#mu} [GeV]")
    Xaxis_pass.SetTitleSize(0.035)
    Xaxis_pass.SetTitleOffset(1.2)
    CMS_lumi(pad_pass, 5, 0, simulation=True)
    pad_pass.Update()


    pad_fail.cd()
    hist_fail.Draw("colz")
    hist_fail.SetTitle("")
    hist_fail.SetTitleSize(0)
    Yaxis_fail = hist_fail.GetYaxis()
    Yaxis_fail.SetTitle("#eta^{#mu}")
    Yaxis_fail.SetTitleSize(0.035)
    Yaxis_fail.SetTitleOffset(1.2)
    Xaxis_fail = hist_fail.GetXaxis()
    Xaxis_fail.SetTitle("p_{T}^{#mu} [GeV]")
    Xaxis_fail.SetTitleSize(0.035)
    Xaxis_fail.SetTitleOffset(1.2)
    CMS_lumi(pad_fail, 5, 0, simulation=True)
    pad_fail.Update()

    c.SaveAs(f"{figpath}/{bkg_cat}_distrib_2d.pdf")

###############################################################################


def plot_projected_bkg(plot_dictionary, binning_pt, binning_eta, flag, logscale=True, 
                       figpath=''):
    """
    """
    bins_pt, bins_eta = binnings[binning_pt], binnings[binning_eta]

    if not (len(bins_pt)==2 or len(bins_eta)==2):
        sys.exit("ERROR: the projection has to be made by collapsing one of the two dimensions")

    if len(bins_pt) > len(bins_eta):
        proj_axis = "pt"
        proj_binning = bins_pt
    else:
        proj_axis = "eta"
        proj_binning = bins_eta

    c = ROOT.TCanvas(f"bkg_{proj_axis}_{flag}", "c", 1200, 1200)
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

    plot_objects = deepcopy(plot_dictionary)

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

        if bkg_key != "TTFullyleptonic_bkg":
            bkg_key_print = deepcopy(bkg_key).replace("_bkg", "") 
        else:
            bkg_key_print = "TTLeptonic"

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
                          ROOT.RooFit.Name(bkg_key_print),
                          ROOT.RooFit.LineColor(colors[bkg_key]),
                          ROOT.RooFit.LineStyle(1),
                          ROOT.RooFit.Normalization(roohisto_bkg.sumEntries(), ROOT.RooAbsReal.NumEvent))
        
        tot_bkg_roohist.add(roohisto_bkg)

    tot_bkg_roohist.plotOn(frame,
                           ROOT.RooFit.Name("Total bkg"),
                           ROOT.RooFit.Binning("plot_binning_total_bkg"),
                           # ROOT.RooFit.Invisible(),
                           ROOT.RooFit.LineColor(colors["bkg_total"]),
                           ROOT.RooFit.MarkerColor(colors["bkg_total"]))
    
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

    
    file = ROOT.TFile("/scratch/rforti/tnp_efficiencies_results/recoplus/bkg_studies/ws_recoplus.root", "READ")
    ws = file.Get("w")

    
    bin_key = "[24.0to26.0][-2.4to-2.3]"

    plot_dict = {
        "axis" : ws.var(f"x_fail_{bin_key}"),
        "WminusJets" : ws.data(f"Minv_bkg_fail_{bin_key}_WplusJets_SS"),
        "Zjets" : ws.data(f"Minv_bkg_fail_{bin_key}_Zjets_SS"),
    }

    plot_bkg_comparison(plot_dict, "fail", bin_key, logscale=True,
                        figpath="./")
