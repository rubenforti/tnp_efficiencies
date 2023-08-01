"""
"""
import ROOT
import sys
import os
import pickle
from array import array
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

def plot_minv_distrib_with_fit(pad, axis, data, pdf_fit):
    """
    Function that creates and saves the plot containing the invariant mass 
    distribution, together with the post-fit pdf. The pull is not implemented
    yet, since there has been some memory problems with the pads.
    """

    pad.cd()

    plot_frame = axis.frame(ROOT.RooFit.Bins(axis.getBins("plot_binning")))
    
    plot_frame.SetTitle("")
    plot_frame.SetTitleSize(0)
    frame_yaxis = plot_frame.GetYaxis()
    frame_yaxis.SetTitleOffset(1.75)

    data.plotOn(plot_frame,
                ROOT.RooFit.Binning("plot_binning"),
                ROOT.RooFit.Name("Data"))
    data.statOn(plot_frame,
                ROOT.RooFit.Label(pad.GetTitle()),
                ROOT.RooFit.What("H"),
                ROOT.RooFit.Layout(0.65, 0.95, 0.95))

    pdf_fit.plotOn(plot_frame,
                   ROOT.RooFit.Name("Fit"),
                   ROOT.RooFit.NormRange("fitRange"),
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
    c = ROOT.TCanvas(f"pf_plot_{bin_key}", f"pf_plot_{bin_key}", 900, 1200)
    c.cd()

    style_settings()

    title_plot_edge = 0.8
    plot_info_edge = 0.25

    pad_title_eff = ROOT.TPad("pad_title_eff", "pad_title_eff", 0, title_plot_edge, 1, 1)
    pad_plot_pass = ROOT.TPad("pad_plot_pass", "Passing probes", 0, plot_info_edge, 0.5, title_plot_edge)
    pad_plot_fail = ROOT.TPad("pad_plot_fail", "Failing probes", 0.5, plot_info_edge, 1, title_plot_edge)
    pad_info_pass = ROOT.TPad("pad_info_pass", "pad_info_pass", 0, 0, 0.5, plot_info_edge)
    pad_info_fail = ROOT.TPad("pad_info_fail", "pad_info_fail", 0.5, 0, 1, plot_info_edge)

    
    pad_title_eff.SetMargin(0.15, 0.15, 0.1, 0.1)
    pad_title_eff.Draw()

    pad_plot_pass.SetMargin(0.15, 0.05, 0.05, 0.05)
    pad_plot_pass.Draw()

    pad_plot_fail.SetMargin(0.15, 0.05, 0.05, 0.05)
    pad_plot_fail.Draw()

    pad_info_pass.SetMargin(0.15, 0.05, 0.05, 0.05)
    pad_info_pass.Draw()

    pad_info_fail.SetMargin(0.15, 0.05, 0.05, 0.05)
    pad_info_fail.Draw()

    pad_title_eff.cd()
    main_info_box = ROOT.TPaveText(0, 0.1, 1, 0.9, "NDC NB")
    main_info_box.SetFillColor(0)
    # titlebox.SetTextFont(42)
    main_info_box.SetTextSize(0.12)
    main_info_box.AddText(f"Bin {bin_key}")
    main_info_box.AddText( "----------------------------------------")
    eff, deff = plot_objects["efficiency"]
    main_info_box.AddText(f"Efficiency = {round(eff*100,2)} #pm {round(deff*100,2)} %")
    eff_mc, deff_mc = plot_objects["efficiency_mc"]
    main_info_box.AddText(f"Efficiency MC = {round(eff_mc*100,2)} #pm {round(deff_mc*100,2)} %")
    scale_factor, dscale_factor = plot_objects["scale_factor"]
    main_info_box.AddText(f"Scale factor = {round(scale_factor,4)} #pm {round(dscale_factor,4)}")
    # print(f"Scale factor = {round(scale_factor, 4)} #pm {round(dscale_factor, 4)}")   
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

   
    pad_info_pass.cd()
    stats_pass = ROOT.TPaveText(0, 0.05, 1., 0.95, "NDC NB")
    stats_pass.SetFillColor(0)
    stats_pass.SetTextSize(0.07)

    [stats_pass.AddText(f"{par.GetTitle()} = {round(par.getVal(),3)} #pm {round(par.getError(),3)}") 
     for par in pass_obj["res"].floatParsFinal()]

    stats_pass.AddText("------------------------------")
    # stats_box.AddText(f"Chi2/ndof = {round(res.chi2(), 2) / pars.getSize()}")
    stats_pass.AddText(f"Status = {pass_obj['res'].status()}")
    stats_pass.AddText(f"Cov quality = {pass_obj['res'].covQual()}")
    stats_pass.AddText(f"Edm = {pass_obj['res'].edm()}")
    stats_pass.Draw()
    c.Update()


    pad_info_fail.cd()
    stats_fail = ROOT.TPaveText(0, 0.05, 1., 0.95, "NDC NB")
    stats_fail.SetFillColor(0)
    stats_fail.SetTextSize(0.07)

    [stats_fail.AddText(f"{par.GetTitle()} = {round(par.getVal(),3)} #pm {round(par.getError(),3)}") 
     for par in fail_obj["res"].floatParsFinal()]

    stats_fail.AddText("------------------------------")
    # stats_box.AddText(f"Chi2/ndof = {round(res.chi2(), 2) / pars.getSize()}")
    stats_fail.AddText(f"Status = {fail_obj['res'].status()}")
    stats_fail.AddText(f"Cov quality = {fail_obj['res'].covQual()}")
    stats_fail.AddText(f"Edm = {fail_obj['res'].edm()}")
    stats_fail.Draw()
    c.Update()

    c.SaveAs(f"{figpath}/fit_pf_{type_analysis}_{bin_key}.pdf")

###############################################################################

def plot_bkg_on_histo(plot_objects, flag, bin_key, logscale=True, figpath=''):
    """
    Function that creates and saves the plot containing a reference histogram 
    (data or signal MC), the distributions of the various bkg processes and the
    total bkg; for graphical choice, the bkg distributions of single processes 
    and the signal MC distribution are plotted as RooHistPdf. The objects to be 
    plotted are contained in the dictionary "plot_objects", which must be 
    structured as follows:
        plot_objects = { "axis" : RooRealVar,
                         "total_bkg" : { "roohisto":RooDataHist, "lumi_scale":float, "integral":float, "color":int },
                         "{bkg_process}_bkg" : { "roohisto":RooDataHist, "histo_pdf":RooHistPdf,
                                                 "lumi_scale":float, "integral":float, "color":int },
                         "MC_signal (optional)" : as above,
                         "data (optional)" : RooDataHist,
                         "fit_pars (optional)" : { "nsig":RooRealVar, "nbkg":RooRealVar },
                         "pdf_bkg_fit (optional)" : {"pdf":RooAbsPdf, "norm":RooRealVar, "color":int }
                        }
            }
    
    """

    imported_data, imported_mc_sig, imported_pdf_bkg = False, False, False

    c = ROOT.TCanvas(f"c_{bin_key}_{flag}", "c", 1200, 900)
    c.cd()

    style_settings()

    pad_title = ROOT.TPad("pad_title", "pad_title", 0, 0.9, 1, 1)
    pad_plot = ROOT.TPad("pad_plot", "pad_plot", 0, 0, 0.7, 0.9)
    # pad_pull = ROOT.TPad("pad_pull", "pad_pull", 0, 0, 0.65, 0.25)
    pad_info = ROOT.TPad("pad_info", "pad_info", 0.7, 0.1, 1, 0.9)

    pad_title.SetMargin(0.1, 0.1, 0.1, 0.1), pad_title.Draw()
    pad_plot.SetMargin(0.15, 0.05, 0.05, 0.05), pad_plot.Draw()
    pad_info.SetMargin(0.05, 0.05, 0.1, 0.9), pad_info.Draw()

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
    frame = axis.frame(ROOT.RooFit.Bins(axis.getBins("plot_binning")))
    frame.SetTitle("")
    frame.SetTitleSize(0)
    frame_yaxis = frame.GetYaxis()
    frame_yaxis.SetTitle("Events / (1 GeV)")
    frame_axis = frame.GetXaxis()
    frame_axis.SetLabelSize(0)
    xtitle = frame_axis.GetTitle()
    frame_axis.SetTitleSize(0)

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
                                        ROOT.RooFit.LineColor(mc_sig_dict["color"]),
                                        ROOT.RooFit.Normalization(mc_sig_dict["integral"], ROOT.RooAbsReal.NumEvent))
        for bin_idx in range(mc_sig_dict["roohisto"].numEntries()):
            if mc_sig_dict["roohisto"].weight(bin_idx) > ctrl_plot_max:
                ctrl_plot_max = mc_sig_dict["roohisto"].weight(bin_idx)

    if "pdf_bkg_fit" in plot_objects.keys():
        imported_pdf_bkg = True
        pdf_bkg_obj = plot_objects.pop("pdf_bkg_fit")
        pdf_bkg_obj["pdf"].plotOn(frame, 
                                ROOT.RooFit.Name("Fitted bkg"),
                                ROOT.RooFit.LineColor(pdf_bkg_obj["color"]),
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
                                 ROOT.RooFit.LineColor(total_bkg["color"]),
                                 ROOT.RooFit.MarkerColor(total_bkg["color"])
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

    
    for bkg_key in plot_objects:
        
        bkg_obj = plot_objects[bkg_key]
        
        for nbins_bkg in [30, 20, 15, 10]:
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
                tmp_histpdf = ROOT.RooHistPdf(f"{bkg_key}_pdf", f"{bkg_key}_pdf", ROOT.RooArgSet(axis), tmp_histo)
                bkg_obj.update({"histo_pdf" : tmp_histpdf})
                break

        
        pad_plot.cd()
        bkg_obj["histo_pdf"].plotOn(frame, 
                         ROOT.RooFit.Name(bkg_key),
                         ROOT.RooFit.LineColor(bkg_obj["color"]),
                         ROOT.RooFit.LineStyle(1),
                         ROOT.RooFit.Normalization(bkg_obj["integral"], ROOT.RooAbsReal.NumEvent))
        pad_info.cd()
        bkg_error = sumw2_error(bkg_obj["roohisto"])
        textbox.AddText(f"Num {bkg_key} = {round(bkg_obj['integral'],2)} #pm {round(bkg_error, 2)}")

    pad_plot.cd()
    pad_plot.SetLogy() if logscale is True else pad_plot.SetLogy(False)
    
    frame.SetMaximum(1.5*ctrl_plot_max)      
    frame.SetMinimum(5e-2)
    frame.Draw()
    pad_plot.Update()
 

    
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
    plot_minv_distrib_with_fit(pad_plot, "pass",axis, data, gaus, pull=False)
    
    c.SaveAs("../prova_plot.pdf")

