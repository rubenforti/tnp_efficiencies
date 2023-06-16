"""
"""
import ROOT
import sys
import os
import pickle
from array import array
from utilities.results_utils import results_manager



def plot_pass_and_fail(axis, data, function, bkg_pdf,
                       name='prova.png', pull=False):
    """
    """

    path = os.path.dirname(__file__)
    ROOT.gSystem.cd(path)

    c = ROOT.TCanvas()
    c.Divide(2)

    c.cd(1)
    ROOT.gPad.SetLeftMargin(0.15)
    frame1 = axis[1].frame(ROOT.RooFit.Title("PASS"))
    data[1].plotOn(frame1)
    
    function[1].plotOn(frame1, ROOT.RooFit.LineColor(4)) 
    
    for comp1 in function[1].getComponents():
        if comp1.GetName() == bkg_pdf[1]:
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
        if comp0.GetName() == bkg_pdf[0]:
            bkg_set = ROOT.RooArgSet(comp0)
            function[0].plotOn(frame0,
                ROOT.RooFit.Components(bkg_set),
                ROOT.RooFit.LineColor(2),
                ROOT.RooFit.LineStyle(ROOT.kDashed))
            
    frame0.Draw()

    c.SaveAs(name)



def plot_distr_with_fit(axis, data, function, bkg_pdf,
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
        if comp.GetName() == bkg_pdf:
            bkg_set = ROOT.RooArgSet(comp)
            function.plotOn(frame,
                ROOT.RooFit.Components(bkg_set),
                ROOT.RooFit.LineColor(2),
                ROOT.RooFit.LineStyle(ROOT.kDashed))
    
    #function.plotOn(frame, ROOT.RooCmdArg("LineColor", 4))
    # bkg_pdf.plotOn(frame, ROOT.RooCmdArg("LineColor", 2))
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
