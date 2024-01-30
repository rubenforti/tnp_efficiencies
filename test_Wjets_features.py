import ROOT
from utilities.base_library import bin_dictionary
import numpy as np
from scipy.stats import ks_2samp as KS
import matplotlib.pyplot as plt
from statistics import median
import sys, os
from array import array
from utilities.base_library import binning

path = os.path.dirname(__file__)
ROOT.gSystem.cd(path)

file = ROOT.TFile("/scratch/rforti/tnp_efficiencies_results/tracking/bkg_figs/OS/ws_tracking_bkg.root")
ws = file.Get("w")

file_SS = ROOT.TFile("/scratch/rforti/tnp_efficiencies_results/tracking/bkg_figs/SS/ws_tracking_bkg_SS.root")
ws_SS = file_SS.Get("w")

print(type(ws), type(ws_SS))

 


bkg_categories = [# "bkg_WW", "bkg_WZ", "bkg_ZZ", 
                  # "bkg_TTSemileptonic", "bkg_TTFullyleptonic", "bkg_Ztautau",
                  "bkg_WplusJets", "bkg_WminusJets"]
                  
bin_dict = bin_dictionary("pt_tracking", "eta")

ratio_mcbkg_samecharge, ratio_OS_SS_wjets = [], []
error_ratio_mcbkg_samecharge, error_ratio_OS_SS_wjets = [], []
pval_KS_samecharge, pval_KS_wjets_SS = [], []

ch = ["plus", "minus"]

NBINS = 40

axis = ROOT.RooRealVar("Minv", "Minv", 50, 130)
axis.setBins(NBINS, "")

h2d_ratio_OS_SS = ROOT.TH2D("h2d_ratio_OS_SS", "SS/OS ratio", 
                            4, binning("pt_tracking"), 48, binning("eta"))
h2d_ratio_OS_SameCharge = ROOT.TH2D("h2d_ratio_OS_SameCharge", "OS/SameCharge ratio",
                                    4, binning("pt_tracking"), 48, binning("eta"))
'''
h2d_KSpval_OS_SS = ROOT.TH2D("h2d_KSpval_OS_SS", "KS p-value OS/SS",
                                4, binning("pt_tracking"), 48, binning("eta"))
h2d_KSpval_OS_SameCharge = ROOT.TH2D("h2d_KSpval_OS_SameCharge", "KS p-value OS/SameCharge",
                                     4, binning("pt_tracking"), 48, binning("eta"))
'''


bin_pt = 1
bin_eta = 1

for b_key in bin_dict.keys():

    nbkg_mc_OS, nbkg_mc_SS = 0, 0
    array_mc_OS, array_mc_SS, array_samecharge,  = [0]*NBINS, [0]*NBINS, [0]*NBINS
    errors_mc_OS, errors_mc_SS = [0]*NBINS, [0]*NBINS

    array_ratio_OS_samecharge, array_ratio_OS_SS = [0]*NBINS, [0]*NBINS
    err_norm_ratio_OS_samecharge, err_norm_ratio_OS_SS = [0]*NBINS, [0]*NBINS

    hist_OS = ROOT.RooDataHist("hist_OS", "hist_OS", ROOT.RooArgList(axis), "")
    hist_SS = ROOT.RooDataHist("hist_SS", "hist_SS", ROOT.RooArgList(axis), "")
    hist_samecharge = ROOT.RooDataHist("hist_samecharge", "hist_samecharge", ROOT.RooArgList(axis), "")

    
    for cat in bkg_categories: 
        data = ws.data(f"Minv_bkg_fail_{b_key}_{cat.replace('bkg_', '')}")
        nbkg_mc_OS += data.sumEntries()
        data_SS = ws_SS.data(f"Minv_bkg_fail_{b_key}_{cat.replace('bkg_', '')}_SS")
        nbkg_mc_SS += data_SS.sumEntries()
        data_samecharge = ws.data(f"Minv_bkg_fail_{b_key}_SameCharge")
        nbkg_samecharge = data_samecharge.sumEntries()

        j=0
        for i in range(80): 
            data.get(i)
            array_mc_OS[j] += data.weight(i)
            errors_mc_OS[j] += data.weightError(ROOT.RooAbsData.SumW2)**2

            hist_OS.get(j)
            hist_OS.set(j, hist_OS.weight(j) + data.weight(i), 
                        (hist_OS.weightError(ROOT.RooAbsData.SumW2)**2 + data.weightError(ROOT.RooAbsData.SumW2)**2)**0.5)

            data_SS.get(i)
            array_mc_SS[j] += data_SS.weight(i)
            errors_mc_SS[j] += data_SS.weightError(ROOT.RooAbsData.SumW2)**2

            hist_SS.get(j)
            hist_SS.set(j, hist_SS.weight(j) + data_SS.weight(i), 
                        (hist_SS.weightError(ROOT.RooAbsData.SumW2)**2 + data_SS.weightError(ROOT.RooAbsData.SumW2)**2)**0.5)

            if i%int(80/NBINS)!=0: j+=1

    j=0
    for i in range(80):
        data_samecharge.get(i)
        array_samecharge[j] += data_samecharge.weight(i)
        hist_samecharge.get(j)
        hist_samecharge.set(j, hist_samecharge.weight(j) + data_samecharge.weight(i), 
                            (hist_samecharge.weightError(ROOT.RooAbsData.SumW2)**2 + data_samecharge.weightError(ROOT.RooAbsData.SumW2)**2)**0.5)
    
        if i%int(80/NBINS)!=0: j+=1

    ratio_mcbkg_samecharge.append(nbkg_mc_OS/nbkg_samecharge)
    ratio_OS_SS_wjets.append(nbkg_mc_OS/nbkg_mc_SS)

    for i in range(NBINS):
        hist_SS.get(i)
        hist_OS.get(i)
        hist_samecharge.get(i)
        # print(i, "", hist_SS.weight(i), hist_OS.weight(i), hist_samecharge.weight(i))

        # RATIO OS/SS
        try:
            array_ratio_OS_SS[i] = hist_OS.weight(i)/hist_SS.weight(i)
        except ZeroDivisionError:
            array_ratio_OS_SS[i] = 0
            err_norm_ratio_OS_SS[i] = 0
            continue
        try:
            err_norm_ratio_OS_SS[i] = array_ratio_OS_SS[i]*(
                (hist_SS.weightError(ROOT.RooAbsData.SumW2)/hist_SS.weight(i))**2 + 
                (hist_OS.weightError(ROOT.RooAbsData.SumW2)/hist_OS.weight(i))**2 )**0.5
            array_ratio_OS_SS[i] = (array_ratio_OS_SS[i]-1)/err_norm_ratio_OS_SS[i]
            err_norm_ratio_OS_SS[i] = 1.
        except ZeroDivisionError:
            err_norm_ratio_OS_SS[i] = 0

        # RATIO OS/SameCharge
        try:
            array_ratio_OS_samecharge[i] = hist_OS.weight(i)/hist_samecharge.weight(i)
        except ZeroDivisionError:
            array_ratio_OS_samecharge[i] = 0
            err_norm_ratio_OS_samecharge[i] = 0
            continue
        try:
            err_norm_ratio_OS_samecharge[i] = array_ratio_OS_samecharge[i]*(
                (hist_samecharge.weightError(ROOT.RooAbsData.Poisson)/hist_samecharge.weight(i))**2 + 
                (hist_OS.weightError(ROOT.RooAbsData.SumW2)/hist_OS.weight(i))**2 )**0.5
            array_ratio_OS_samecharge[i] = (array_ratio_OS_samecharge[i]-1)/err_norm_ratio_OS_samecharge[i]
            err_norm_ratio_OS_samecharge[i] = 1.
        except ZeroDivisionError:
            err_norm_ratio_OS_samecharge[i] = 0    
    

    roopdf_OS = ROOT.RooHistPdf("roopdf_OS", "roopdf_OS", ROOT.RooArgSet(axis), hist_OS)
    roopdf_SS = ROOT.RooHistPdf("roopdf_SS", "roopdf_SS", ROOT.RooArgSet(axis), hist_SS)
    roopdf_samecharge = ROOT.RooHistPdf("roopdf_samecharge", "roopdf_samecharge", ROOT.RooArgSet(axis), hist_samecharge)

    '''
    dataset_OS = roopdf_OS.generate(ROOT.RooArgSet(axis), NumEvents=int(nbkg_mc_OS))
    dataset_SS = roopdf_SS.generate(ROOT.RooArgSet(axis), NumEvents=int(nbkg_mc_SS))
    dataset_samecharge = roopdf_samecharge.generate(ROOT.RooArgSet(axis), NumEvents=int(nbkg_samecharge))

    array_OS_np, array_SS_np, array_samecharge_np = [], [], []

    for i in range(dataset_OS.numEntries()): 
        ev = dataset_OS.get(i)
        var = ev.find("Minv")        
        array_OS_np.append(var.getVal())
    for i in range(dataset_SS.numEntries()): 
        ev = dataset_SS.get(i)
        var = ev.find("Minv")
        array_SS_np.append(var.getVal())
    for i in range(dataset_samecharge.numEntries()): 
        ev = dataset_samecharge.get(i)
        var = ev.find("Minv")
        array_samecharge_np.append(var.getVal())

    array_OS_np = np.array(array_OS_np)
    array_SS_np = np.array(array_SS_np)
    array_samecharge_np = np.array(array_samecharge_np)
    '''

    '''
    array_OS_np = np.array(array_mc_OS).sort()
    array_SS_np = np.array(array_mc_SS).sort()
    array_samecharge_np = np.array(array_samecharge).sort()
    '''

    # pval_KS_samecharge.append(KS(array_OS_np, array_samecharge_np, alternative="twosided").pvalue)
    # pval_KS_wjets_SS.append(KS(array_OS_np, array_SS_np, alternative="twosided").pvalue)

    error_ratio_mcbkg_samecharge.append(ratio_mcbkg_samecharge[-1]*(
        (data_samecharge.sumEntriesW2()/(nbkg_samecharge**2)) + (data.sumEntriesW2()/(nbkg_mc_OS**2)))**0.5)
    
    error_ratio_OS_SS_wjets.append(ratio_OS_SS_wjets[-1]*(
        (data_SS.sumEntriesW2()/(nbkg_mc_SS**2)) + (data.sumEntriesW2()/(nbkg_mc_OS**2)))**0.5)

    h2d_ratio_OS_SS.SetBinContent(bin_pt, bin_eta, ratio_OS_SS_wjets[-1])
    h2d_ratio_OS_SS.SetBinError(bin_pt, bin_eta, error_ratio_OS_SS_wjets[-1])
    h2d_ratio_OS_SameCharge.SetBinContent(bin_pt, bin_eta, ratio_mcbkg_samecharge[-1])
    h2d_ratio_OS_SameCharge.SetBinError(bin_pt, bin_eta, error_ratio_mcbkg_samecharge[-1])
    # h2d_KSpval_OS_SS.SetBinContent(bin_pt, bin_eta, pval_KS_wjets_SS[-1])
    # h2d_KSpval_OS_SameCharge.SetBinContent(bin_pt, bin_eta, pval_KS_samecharge[-1])
    

    c = ROOT.TCanvas(f"c_{b_key}", f"c_{b_key}", 900, 1200)
    c.cd()

    pad_info = ROOT.TPad("pad_info", "pad_info", 0, 0.9, 1, 1)
    pad_info.SetMargin(0.25, 0.25, 0.2, 0.1)
    pad_info.Draw()


    pad_plot = ROOT.TPad("pad_plot", "pad_plot", 0, 0.4, 1, 0.9)
    pad_plot.SetMargin(0.1, 0.05, 0.1, 0.01)
    pad_plot.Draw()

    pad_ratio_samecharge = ROOT.TPad("pad_ratio_samecharge", "pad_ratio_samecharge", 0, 0.2, 1, 0.4)
    pad_ratio_samecharge.SetMargin(0.1, 0.05, 0.1, 0.01)
    pad_ratio_samecharge.Draw()

    pad_ratio_SS = ROOT.TPad("pad_ratio_SS", "pad_ratio_SS", 0, 0, 1, 0.2)
    pad_ratio_SS.SetMargin(0.1, 0.05, 0.1, 0.01)
    pad_ratio_SS.Draw()

    pad_info.cd()
    text = ROOT.TPaveText(0, 0.2, 1, 0.8, "NDC NB")
    text.AddText(f"Bin: {b_key}")
    text.AddText(f"Ratio: \tOS/SS = {nbkg_mc_OS/nbkg_mc_SS:.3f},  OS/SameCharge = {nbkg_mc_OS/nbkg_samecharge:.3f}")
    #text.AddText(0.5, 0.1, f"KS test p-value: \tOSvSS = {pval_KS_wjets_SS[-1]:.3f},  OSvSameCharge = {pval_KS_samecharge[-1]:.3f}")
    text.Draw()
    c.Update()  

    pad_plot.cd()
    frame = axis.frame()
    frame.SetTitle("")
    frame.GetYaxis().SetTitle("Events/(1 GeV)")
    frame.GetYaxis().SetTitleOffset(1.2)

    hist_samecharge.plotOn(frame, 
                           ROOT.RooFit.DataError(ROOT.RooAbsData.Poisson), 
                           ROOT.RooFit.MarkerColor(ROOT.kBlack))
    roopdf_OS.plotOn(frame, 
                     ROOT.RooFit.LineColor(ROOT.kRed), 
                     ROOT.RooFit.Normalization(nbkg_mc_OS, ROOT.RooAbsReal.NumEvent))
    roopdf_SS.plotOn(frame, 
                     ROOT.RooFit.LineColor(ROOT.kOrange), 
                     ROOT.RooFit.Normalization(nbkg_mc_SS, ROOT.RooAbsReal.NumEvent))
    
    frame.Draw()

    pad_ratio_samecharge.cd()
    ratio_graph_sc = ROOT.TGraphErrors(NBINS,
                            array("d", [50.5+(i*int(80/NBINS)) for i in range(NBINS)]), array("d", array_ratio_OS_samecharge), 
                            array("d", [0]*NBINS), array("d", err_norm_ratio_OS_samecharge))
    ratio_graph_sc.GetXaxis().SetRangeUser(50, 130)
    ratio_graph_sc.SetTitle("")
    ratio_graph_sc.SetMarkerStyle(20)
    ratio_graph_sc.SetMarkerSize(0.5)
    ratio_graph_sc.GetYaxis().SetTitle("(OS/SameCharge - 1)/#sigma")
    #ratio_graph.GetYaxis().SetTitleSize(0.5)
    ratio_graph_sc.Draw("ZAP")
    hline_sc = ROOT.TLine(50, 0, 130, 0)
    hline_sc.SetLineStyle(2)
    hline_sc.Draw()

    pad_ratio_SS.cd()
    ratio_graph_SS = ROOT.TGraphErrors(NBINS,
                            array("d", [50.5+(i*int(80/NBINS)) for i in range(NBINS)]), array("d", array_ratio_OS_SS), 
                            array("d", [0]*NBINS), array("d", err_norm_ratio_OS_SS))
    ratio_graph_SS.GetXaxis().SetRangeUser(50, 130)
    ratio_graph_SS.SetTitle("")
    ratio_graph_SS.SetMarkerStyle(20)
    ratio_graph_SS.SetMarkerSize(0.5)
    ratio_graph_SS.GetYaxis().SetTitle("(OS/SS - 1)/#sigma")
    #ratio_graph.GetYaxis().SetTitleSize(0.5)
    ratio_graph_SS.Draw("ZAP")
    hline = ROOT.TLine(50, 0, 130, 0)
    hline.SetLineStyle(2)
    hline.Draw()



    c.SaveAs(f"wjets_tests/bin_plots/{b_key}_wjets_ratio.pdf")


    if bin_eta%48==0: 
        bin_eta = 1
        bin_pt += 1
    else:
        bin_eta += 1

    

ratio_mcbkg_samecharge = np.array(ratio_mcbkg_samecharge)
ratio_OS_SS_wjets = np.array(ratio_OS_SS_wjets)
# pval_KS_samecharge = np.array(pval_KS_samecharge)
# pval_KS_wjets_SS = np.array(pval_KS_wjets_SS)


print("Ratio SameCharge/mcBkg:", ratio_mcbkg_samecharge.mean(), ratio_mcbkg_samecharge.max(), ratio_mcbkg_samecharge.min())
print("Ratio SS/OS mcBkg:", ratio_OS_SS_wjets.mean(), ratio_OS_SS_wjets.max(), ratio_OS_SS_wjets.min()) 
# print("Pval KS_test SameCharge/mcBkg (OS):", pval_KS_samecharge.mean(), pval_KS_samecharge.max(), pval_KS_samecharge.min())
# print("Pval KS_test SS/OS mcBkg:", pval_KS_wjets_SS.mean(), pval_KS_wjets_SS.max(), pval_KS_wjets_SS.min())


file_out = ROOT.TFile("wjets_tests/histograms.root", "RECREATE")
file_out.cd()
h2d_ratio_OS_SS.Write()
h2d_ratio_OS_SameCharge.Write()
# h2d_KSpval_OS_SS.Write()
# h2d_KSpval_OS_SameCharge.Write()
file_out.Close()
