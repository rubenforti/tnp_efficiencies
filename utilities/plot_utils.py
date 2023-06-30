"""
"""
import ROOT
import sys
import os
import pickle
from array import array
from utilities.results_utils import results_manager
from utilities.base_library import bkg_lumi_scales
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




def plot_pass_and_fail(axis, data, function, pdf_bkg, name='prova.png', pull=False):
    """
    """

    # path = os.path.dirname(__file__)
    # ROOT.gSystem.cd(path)

    c = ROOT.TCanvas()
    c.Divide(2)

    c.cd(1)
    ROOT.gPad.SetLeftMargin(0.15)
    frame1 = axis[1].frame(ROOT.RooFit.Title("PASS"))
    data[1].plotOn(frame1)
    
    function[1].plotOn(frame1, ROOT.RooFit.LineColor(4)) 
    
    for comp1 in function[1].getComponents():
        if comp1.GetName() == pdf_bkg[1]:
            bkg_set = ROOT.RooArgSet(comp1)
            function[1].plotOn(frame1, 
                               ROOT.RooFit.Components(bkg_set),
                               ROOT.RooFit.LineColor(2),
                               ROOT.RooFit.LineStyle(ROOT.kDashed))
    frame1.Draw()

    c.cd(2)
    ROOT.gPad.SetLeftMargin(0.15)
    frame0 = axis[0].frame(ROOT.RooFit.Title("FAIL"))
    data[0].plotOn(frame0)
    
    function[0].plotOn(frame0, ROOT.RooFit.LineColor(4))
    
    for comp0 in function[0].getComponents():
        if comp0.GetName() == pdf_bkg[0]:
            bkg_set = ROOT.RooArgSet(comp0)
            function[0].plotOn(frame0,
                ROOT.RooFit.Components(bkg_set),
                ROOT.RooFit.LineColor(2),
                ROOT.RooFit.LineStyle(ROOT.kDashed))
            
    frame0.Draw()

    c.SaveAs(name)



def plot_distr_with_fit(axis, data, function, pdf_bkg,
                        name='prova.png', pull=False):
    """
    """

    path = os.path.dirname(__file__)
    ROOT.gSystem.cd(path)

    c = ROOT.TCanvas()
    if pull is True:
        c.Divide(2)

    c.cd(1)
    ROOT.gPad.SetLeftMargin(0.15)
    frame = axis.frame(ROOT.RooFit.Title("TITLE"))
    data.plotOn(frame)
    
    function.plotOn(frame, ROOT.RooFit.LineColor(4))
        
    
    for comp in function.getComponents():
        if comp.GetName() == pdf_bkg:
            bkg_set = ROOT.RooArgSet(comp)
            function.plotOn(frame,
                ROOT.RooFit.Components(bkg_set),
                ROOT.RooFit.LineColor(2),
                ROOT.RooFit.LineStyle(ROOT.kDashed))
    
    #function.plotOn(frame, ROOT.RooCmdArg("LineColor", 4))
    # pdf_bkg.plotOn(frame, ROOT.RooCmdArg("LineColor", 2))
    frame.Draw()

    if pull:
        c.cd(2)
        ROOT.gPad.SetLeftMargin(0.15)
        hpull = frame.pullHist()
        frame2 = axis.frame(Title="Pull Distribution")
        frame2.addPlotable(hpull, "P")
        frame2.Draw()
    c.SaveAs(name)


def plot_results(filename, results, binning_pt=(), binning_eta=()):
    """
    """

    file = ROOT.TFile.Open(filename, "UPDATE")

    if len(binning_pt) == 0:
        binning_pt = array('d', [24., 26., 28., 30., 32., 34.,
                           36., 38., 40., 42., 44., 47., 50., 55., 60., 65.])
    if len(binning_eta) == 0:
        binning_eta = array('d', [round(-2.4 + i*0.1, 2) for i in range(49)])

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


def plot_bkg_on_data(plot_objects, flag, bin_key, figpath=''):
    """
    """

    c = ROOT.TCanvas("c", "c", 900, 900)
    c.cd()

    style_settings()

    pad_title = ROOT.TPad("pad_title", "pad_title", 0, 0.95, 1, 1)
    pad_plot = ROOT.TPad("pad_plot", "pad_plot", 0, 0.25, 0.7, 0.95)
    pad_pull = ROOT.TPad("pad_pull", "pad_pull", 0, 0, 0.7, 0.25)
    pad_info = ROOT.TPad("pad_info", "pad_info", 0.7, 0, 1, 0.9)

    pad_title.SetLeftMargin(0.1)
    pad_title.SetRightMargin(0.1)
    pad_title.SetTopMargin(0.1)
    pad_title.SetBottomMargin(0.1)
    pad_title.Draw()

    pad_plot.SetBottomMargin(0)
    pad_plot.SetLeftMargin(0.15)
    pad_plot.SetRightMargin(0.05)
    pad_plot.SetTopMargin(0.05)
    pad_plot.Draw()

    pad_pull.SetTopMargin(0)
    pad_pull.SetBottomMargin(0.2)
    pad_pull.SetLeftMargin(0.15)
    pad_pull.SetRightMargin(0.05)
    pad_pull.Draw()

    pad_info.SetLeftMargin(0)
    pad_info.SetRightMargin(0)
    pad_info.SetTopMargin(0.05)
    pad_info.SetBottomMargin(0)
    pad_info.Draw()


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
    frame = axis.frame(ROOT.RooFit.Bins(60))
    frame.SetTitle("")
    frame.SetTitleSize(0)
    frame_axis = frame.GetXaxis()
    frame_axis.SetLabelSize(0)
    xtitle = frame_axis.GetTitle()
    frame_axis.SetTitleSize(0)

    datahist = plot_objects.pop("data")
    datahist.plotOn(frame, ROOT.RooFit.Name("Data"))

    pdf_bkg_obj = plot_objects.pop("pdf_bkg_fit")
    pdf_bkg_obj["pdf"].plotOn(frame, 
                              ROOT.RooFit.Name("Fitted bkg"),
                              ROOT.RooFit.LineColor(pdf_bkg_obj["color"]),
                              ROOT.RooFit.Normalization(pdf_bkg_obj["norm"].getVal(), ROOT.RooAbsReal.NumEvent))

    total_bkg = plot_objects.pop("total_bkg")
    total_bkg["roohisto"].plotOn(frame,
                                 ROOT.RooFit.Name("Total bkg"),
                                 # ROOT.RooFit.Invisible(),
                                 ROOT.RooFit.LineColor(total_bkg["color"]),
                                 ROOT.RooFit.MarkerColor(total_bkg["color"]),
                                 ROOT.RooFit.Normalization(total_bkg["integral"]["val"], ROOT.RooAbsReal.NumEvent))
    
    pad_info.cd()
    textbox = ROOT.TPaveText(0, 0.4, 1, 0.8, "NDC NB")
    textbox.SetFillColor(0)
    textbox.SetTextSize(0.05)
    textbox.SetTextAlign(12)
    c.Update()
    textbox.AddText(f"Data entries: {datahist.sumEntries()}")
    textbox.AddText(f"Nbkg from fit: {round(pdf_bkg_obj['norm'].getVal(),2)} #pm {round(pdf_bkg_obj['norm'].getError(),2)}")
    textbox.AddText(f"Nbkg from MC: {round(total_bkg['integral']['val'],2)} #pm {round(total_bkg['integral']['err'],2)}")
    difference = total_bkg['integral']['val'] - pdf_bkg_obj['norm'].getVal()
    sigma_difference = ((total_bkg['integral']['err']**2)+(pdf_bkg_obj['norm'].getError()**2))**0.5
    nsigma = difference/sigma_difference
    textbox.AddText(f"Distance in #sigma = {round(nsigma, 2)}")
    textbox.AddText("--------------------")

   
    for bkg_key in plot_objects:
        pad_plot.cd()
        bkg_obj = plot_objects[bkg_key]
        bkg_obj["histo_pdf"].plotOn(frame, 
                         ROOT.RooFit.Name(bkg_key),
                         ROOT.RooFit.LineColor(bkg_obj["color"]),
                         ROOT.RooFit.LineStyle(1),
                         ROOT.RooFit.Normalization(bkg_obj["integral"], ROOT.RooAbsReal.NumEvent))
        pad_info.cd()
        bkg_error = (bkg_obj["lumi_scale"]*bkg_obj["integral"])**0.5
        textbox.AddText(f"Num {bkg_key} = {round(bkg_obj['integral'],2)} #pm {round(bkg_error, 2)}")


    pad_plot.cd()
    frame.Draw()
    CMS_lumi(pad_plot, 5, 0, simulation=True)

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

    print(total_bkg["integral"]["val"], pdf_bkg_obj["norm"].getVal())

    c.SaveAs(f"{figpath}/bkg_{bin_key}_{flag}.pdf")

    #c.Destructor()
        





if __name__ == '__main__':

    
    '''
    file = ROOT.TFile('results/benchmark_iso/ws_iso_indep.root', "READ")

    ws = file.Get("w")
    '''
    res = results_manager("indep")
    res.Open("results/benchmark_iso_sim/results_iso_sim.pkl")
    '''
    for bin_pt in range(1, 16):
        for bin_eta in range(1, 49):
            res.add_result(ws, bin_pt, bin_eta, check=False)
    '''
    # res.Write('results/benchmark_iso/new_results_2.pkl')

    plot_results("root_files/ws/ws_iso_indep.root", res)
