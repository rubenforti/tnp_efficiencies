"""
"""

import ROOT
from utilities.base_library import bkg_lumi_scales, binning, bin_dict, get_new_binning
from utilities.dataset_utils import ws_init
from array import array
import sys


def initialize_h2d_bkg(bkg_categories, binning_pt, binning_eta):

    histos_pass = {}
    histos_fail = {}

    if "WW" in bkg_categories:
        h2d_WW_pass = ROOT.TH2D(
            "WW_pass", "WW_pass", len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta)
        h2d_WW_fail = ROOT.TH2D(
            "WW_fail", "WW_fail", len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta)
        histos_pass.update({"WW" : h2d_WW_pass})
        histos_fail.update({"WW" : h2d_WW_fail})

    if "WZ" in bkg_categories:
        h2d_WZ_pass = ROOT.TH2D(
            "WZ_pass", "WZ_pass", len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta)
        h2d_WZ_fail = ROOT.TH2D(
            "WZ_fail", "WZ_fail", len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta)
        histos_pass.update({"WZ" : h2d_WZ_pass})
        histos_fail.update({"WZ" : h2d_WZ_fail})

    if "ZZ" in bkg_categories:
        h2d_ZZ_pass = ROOT.TH2D(
            "ZZ_pass", "ZZ_pass", len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta)
        h2d_ZZ_fail = ROOT.TH2D(
            "ZZ_fail", "ZZ_fail", len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta)
        histos_pass.update({"ZZ" : h2d_ZZ_pass})
        histos_fail.update({"ZZ" : h2d_ZZ_fail})

    if "TTSemileptonic" in bkg_categories:
        h2d_TTsemileptonic_pass = ROOT.TH2D("TTsemileptonic_pass", "TTsemileptonic_pass", 
                                            len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta)   
        h2d_TTsemileptonic_fail = ROOT.TH2D("TTsemileptonic_fail", "TTsemileptonic_fail", 
                                            len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta) 
        histos_pass.update({"TTSemileptonic" : h2d_TTsemileptonic_pass})
        histos_fail.update({"TTSemileptonic" : h2d_TTsemileptonic_fail})
    
    if "Ztautau" in bkg_categories:
        h2d_Ztautau_pass = ROOT.TH2D("Ztautau_pass", "Ztautau_pass", 
                                     len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta)
        h2d_Ztautau_fail = ROOT.TH2D("Ztautau_fail", "Ztautau_fail", 
                                     len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta)
        histos_pass.update({"Ztautau" : h2d_Ztautau_pass})
        histos_fail.update({"Ztautau" : h2d_Ztautau_fail})
    
    if "SameCharge" in bkg_categories:
        h2d_samecharge_pass = ROOT.TH2D("SameCharge_pass", "SameCharge_pass", 
                                     len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta)
        h2d_samecharge_fail = ROOT.TH2D("SameCharge_fail", "SameCharge_fail", 
                                     len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta)
        histos_pass.update({"SameCharge" : h2d_samecharge_pass})
        histos_fail.update({"SameCharge" : h2d_samecharge_fail})

    return histos_pass, histos_fail

###############################################################################

def plot_bkg_on_data(axis, histo_data, histos_bkg, bin_keys):

    c = ROOT.TCanvas("c", "c", 1200, 900)
    c.cd()

    frame = axis.frame(ROOT.RooFit.Title("Bkg on data"))
    frame.GetXaxis().SetTitle("M_{TP} [GeV]")

    legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)

    legend.SetFillColor(ROOT.kWhite)
    legend.SetLineColor(ROOT.kWhite)

    histo_data = ROOT.RooDataHist("data", "data", axis, histo_data)

    colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen, ROOT.kOrange, ROOT.kMagenta, ROOT.kCyan, ROOT.kYellow]
    col_idx = 0
    histo_data.plotOn(frame, ROOT.RooFit.Name("data"))

    histo_bkg = histos_bkg["WW"]
    print(type(histo_bkg))
    tot_bkg = ROOT.RooDataHist("tot_bkg", "tot_bkg", axis, histo_bkg)

    '''
    for bkg in histos_bkg:

        histos_bkg[bkg] = ROOT.RooDataHist("bkg", "bkg", axis, histos_bkg[bkg])
        histos_bkg[bkg].plotOn(frame,
                               ROOT.RooFit.Name(bkg), 
                               ROOT.RooFit.MarkerSize(0),
                               ROOT.RooFit.LineColor(colors[col_idx]),
                               ROOT.RooFit.DataError(3))
        tot_bkg.add(histos_bkg[bkg])
        col_idx += 1
    '''

    
    tot_bkg.plotOn(frame, 
                   ROOT.RooFit.Name("tot_bkg"),
                   # ROOT.RooFit.MarkerSize(0), 
                   ROOT.RooFit.LineColor(ROOT.kBlack),
                   ROOT.RooFit.DataError(3))
    
    frame.Draw()
    c.SaveAs("plots/bkg_on_data.pdf")

###############################################################################

def draw_bkg_distrib_2d(workspace, bkg_categories, lumi_scales, bin_dict, binning_pt, binning_eta, 
                        divide_for_data=False, save_root_file=False):


    histos_pass, histos_fail = initialize_h2d_bkg(bkg_categories, binning_pt, binning_eta)

    nbins_pt = len(binning_pt) - 1
    nbins_eta = len(binning_eta) - 1

    rootfile_distrib = ROOT.TFile("root_files/bkg_2d_distributions.root", "RECREATE")

    for cat in bkg_categories:

        if divide_for_data:
            histos_pass[cat].SetTitle(f"{histos_pass[cat].GetTitle()} norm on data")
            histos_fail[cat].SetTitle(f"{histos_fail[cat].GetTitle()} norm on data")

        print(len(binning_pt), len(binning_eta))
        for bin_key in bin_dict:

            _, bin_pt, bin_eta = bin_dict[bin_key]

            print(bin_pt)
            print(bin_eta)

            h_pass = workspace.data(f"Minv_bkg_pass_{bin_key}_{cat}")
            h_fail = workspace.data(f"Minv_bkg_fail_{bin_key}_{cat}")
            n_pass = h_pass.sumEntries()
            n_fail = h_fail.sumEntries()
                
            if divide_for_data:
                h_data_pass = workspace.data(f"Minv_data_pass_{bin_key}")
                h_data_fail = workspace.data(f"Minv_data_fail_{bin_key}")
                # print(n_pass, h_data_pass.sumEntries())
                n_pass, n_fail = n_pass/h_data_pass.sumEntries(), n_fail/h_data_fail.sumEntries()

            # Bin transformation needed in case the bins are merged
            if type(bin_pt) is list:
                bin_pt = int(1+(nbins_pt*(bin_pt[0]-1)/15.))
            if type(bin_eta) is list:
                bin_eta = int(1+(nbins_eta*(bin_eta[0]-1)/48.))

            histos_pass[cat].SetBinContent(bin_pt, bin_eta, n_pass*1.0)
            # The histograms are weighted with the lumi-scale, so the error on the 
            # total number of events is the sqrt of the resulting number 
            # (= sqrt(lumi_scale * num_init)) multiplied with sqrt(lumi_scale), to obtain the 
            # correct value of lumi_scale*sqrt(num_init)
            if n_pass<0 or n_fail<0 or lumi_scales[cat]<0:
                print(bin_pt, bin_eta, cat)
                print("ERROR: negative number of events")
                sys.exit()
            histos_pass[cat].SetBinError(bin_pt, bin_eta, (n_pass*lumi_scales[cat])**0.5)  
            histos_fail[cat].SetBinContent(bin_pt, bin_eta, n_fail*1.0)
            histos_fail[cat].SetBinError(bin_pt, bin_eta, (n_fail*lumi_scales[cat])**0.5)
        
        c = ROOT.TCanvas("c", "c", 1200, 1200)
        c.Divide(1,2)

        c.cd(1)
        ROOT.gStyle.SetOptStat("")
        histos_pass[cat].Draw("COLZ")
        c.cd(2)
        ROOT.gStyle.SetOptStat("")
        histos_fail[cat].Draw("COLZ")

        if divide_for_data:
            c.SaveAs(f"figs/backgrounds/{cat}_bkg_h2d_norm.pdf")
        else:
            c.SaveAs(f"figs/backgrounds/{cat}_bkg_h2d.pdf")

        if save_root_file:
            rootfile_distrib.cd()
            histos_pass[cat].Write()
            histos_fail[cat].Write()
    
    if save_root_file:
        rootfile_distrib.Close()

    # return histos_pass, histos_fail
       
###############################################################################
###############################################################################


if __name__ == "__main__":


    bkg_categories = ["WW", "WZ", "ZZ", "TTSemileptonic", "Ztautau"]
    # bkg_categories = ["Ztautau"]

    type_eff = ("sa", "global", "ID", "iso", "trigger", "veto")
    t = type_eff[3]

    types_analysis = ["indep", "sim"]
    an = types_analysis[0]

    lumi_scales = bkg_lumi_scales(t, bkg_categories)
    # lumi_scales = {"WW":1, "WZ":1, "ZZ":1, "TTSemileptonic":1, "Ztautau":1}

    filename_data = "/scratchnvme/wmass/Steve_root_files/Standard_SF_files/tnp_iso_data_vertexWeights1_oscharge1.root"
    filename_mc = "/scratchnvme/wmass/Steve_root_files/Standard_SF_files/tnp_iso_mc_vertexWeights1_oscharge1.root"
    dirname_bkg = "/scratchnvme/rajarshi/Bkg_TNP_3D_Histograms/OS"

    bkg_filepaths = {}
    [bkg_filepaths.update({cat : 
        f"{dirname_bkg}/tnp_{t}_{cat}_vertexWeights1_oscharge1.root"}) for cat in bkg_categories]


    import_dictionary = {
        "data" : filename_data,
        # "mc" : filename_mc,
        "bkg" : {
            "filepaths" : bkg_filepaths,
            "lumi_scales" : lumi_scales
        }
    }
    
    binning_eta_key = "eta_8bins"
    binning_pt_key = "pt"
    binning_mass = binning("mass_60_120")

    new_binning_pt = binning(binning_pt_key)
    new_binning_eta = binning(binning_eta_key)

    # bin_set = bin_dict()
    bin_set = get_new_binning(new_binning_pt, new_binning_eta)

    ws = ws_init(import_dictionary, an, bin_set, binning_mass)

    draw_bkg_distrib_2d(ws, bkg_categories, lumi_scales, bin_set, binning(binning_pt_key), binning(binning_eta_key), 
                        divide_for_data=False, save_root_file=True)


    '''
    h_pass_gen = ROOT.TH2D(h_pass[0])
    h_pass_gen.Sumw2()
    print(h_pass[0].GetBinContent(1,1), h_pass[0].GetBinError(1,1))

    h_fail_gen = ROOT.TH2D(h_fail[0])
    h_fail_gen.Sumw2()

    for i in range(1, 5):
        print(h_pass[i].GetBinContent(1,1), h_pass[i].GetBinError(1,1))
        h_pass_gen.Add(h_pass[i])
        h_fail_gen.Add(h_fail[i])
    

    c = ROOT.TCanvas()
    c.Divide(1,2)
    c.cd(1)
    ROOT.gStyle.SetOptStat("en")
    h_pass_gen.Draw("COLZ")
    c.cd(2)
    ROOT.gStyle.SetOptStat("en")
    h_fail_gen.Draw("COLZ")
    '''

    '''
    res = results_manager("indep")
    res.Open("results/benchmark_iso/results_iso_indep_benchmark.pkl")
    res_dict = res.dictionary()

    print("*********")
    # print(h_pass_gen.GetBinContent(1,1), h_pass_gen.GetBinError(1,1))
    print(res_dict["1,1"]["pars_pass"].find("nbkg_pass_(1|1)"))
    print(res_dict["1,1"]["pars_fail"].find("nbkg_fail_(1|1)"))
    '''













