"""
"""

import ROOT
from workspace_config import get_roohist
from results_utils import results_manager
from array import array
import sys


def initialize_bkg_h2d(bkg_categories, binning_pt, binning_eta):

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


BR_TAUToMU = 0.1739
BR_TAUToE = 0.1782
Z_TAU_TO_LEP_RATIO = (1.-(1. - BR_TAUToMU - BR_TAUToE)**2)
xsec_ZmmPostVFP = 2001.9

cross_sections_bkg = {
    # Unit = pb
    "WW" : 12.6,
    "WZ" : 5.4341,
    "ZZ" : 0.60,
    "TTSemileptonic" : 366.34,
    "Ztautau" : xsec_ZmmPostVFP*Z_TAU_TO_LEP_RATIO
}


def plot_bkg_mass():

    pass

def get_bkg_lumi_scale(type_eff, bkg_categories):
    """
    "Lumi scale" defined as alpha that satisifies alpha*lumi_bkg = lumi_data
    """

    lumi_data = 16.8  # fb^-1

    lumi_scale = {}

    for cat in bkg_categories:
        file = ROOT.TFile(f"/scratchnvme/rajarshi/Bkg_TNP_3D_Histograms/OS/tnp_{type_eff}_{cat}_vertexWeights1_oscharge1.root")
        
        wsum_histo = file.Get("weightSum")
        num_init = wsum_histo.Integral()
        xsection = cross_sections_bkg[cat]*1000
        lumi_bkg = num_init/xsection

        scale = lumi_data/lumi_bkg

        lumi_scale.update({cat : scale})
    
    return lumi_scale



def plot_bkg_on_data(ws, ws_bkg, merge_bins=False, new_binning_pt=[], new_binning_eta=[]):


    pass
    



def generate_bkg_datasets(type_eff, bkg_categories, lumi_scales, binning_pt, binning_eta, 
                          merge_bins=False, new_binning_pt=[], new_binning_eta=[]):

    axis = ROOT.RooRealVar("x_bkg_pass", "x", 50.0, 130.0)  # If we are going to convolve the bkg distributions, we will need a specific axis for every histogram in every bin

    ws = ROOT.RooWorkspace("w")
    ws.Import(axis)
    
    for cat in bkg_categories:
        
        file = ROOT.TFile(f"/scratchnvme/rajarshi/Bkg_TNP_3D_Histograms/OS/tnp_{type_eff}_{cat}_vertexWeights1_oscharge1.root")


        if merge_bins:
            
            new_pt_bins_idx = [-1]*len(binning_pt)
            idx_pt = 0
            for bound_pt in binning_pt:
                if bound_pt in new_binning_pt:
                    new_pt_bins_idx[idx_pt] = bound_pt
                idx_pt +=1
        
            new_eta_bins_idx = [-1]*len(binning_eta)
            idx_eta = 0
            for bound_eta in binning_eta:
                if bound_eta in new_binning_eta:
                    new_eta_bins_idx[idx_eta] = bound_eta
                idx_eta += 1
        else:
            new_pt_bins_idx = binning_pt
            new_eta_bins_idx = binning_eta
    
        print(new_binning_pt)
        print(new_binning_eta)

        for i in range(1, len(binning_pt)):
            
            if new_pt_bins_idx[i-1] == -1:
                continue
            ii = i
            while new_pt_bins_idx[ii] == -1:
                ii += 1
            bin_pt_idx_list = [i, ii]

            for j in range(1, len(binning_eta)):

                if new_eta_bins_idx[j-1] == -1:
                    continue
                jj = j
                while new_eta_bins_idx[jj] == -1:
                    jj += 1
                bin_eta_idx_list = [j, jj]

                h_pass = get_roohist(file, "bkg", "pass", axis, bin_pt_idx_list, bin_eta_idx_list, global_scale=lumi_scales[cat])
                h_pass.SetName(f"{h_pass.GetName()}_{cat}")
                ws.Import(h_pass)

                h_fail = get_roohist(file, "bkg", "fail", axis, [i], [j], global_scale=lumi_scales[cat])
                h_fail.SetName(f"{h_fail.GetName()}_{cat}")
                ws.Import(h_fail)
            
    if not merge_bins:
        ws.writeToFile(f"root_files/ws/ws_backgrounds_{type_eff}.root")
    else:
        ws.writeToFile(f"root_files/ws/ws_backgrounds_{type_eff}_mergebins.root")

    return ws





def draw_bkg_distrib_2d(workspace_bkg, bkg_categories, lumi_scales, binning_pt, binning_eta,
                        divide_for_data=False, file_ws_data='', save_root_file=False):

    histos_pass, histos_fail = initialize_bkg_h2d(bkg_categories, binning_pt, binning_eta)

    if save_root_file:
        rootfile_distrib = ROOT.TFile("root_files/bkg_2d_distributions.root", "RECREATE")

    for cat in bkg_categories:

        if divide_for_data:
            histos_pass[cat].SetTitle(f"{histos_pass[cat].GetTitle()} norm on data")
            histos_fail[cat].SetTitle(f"{histos_fail[cat].GetTitle()} norm on data")
            file_data = ROOT.TFile(file_ws_data)
            ws_data = file_data.Get("w")

        for i in range(1, len(binning_pt)):
            for j in range(1, len(binning_eta)):
                h_pass = workspace_bkg.data(f"Minv_bkg_pass_({i}|{j})_{cat}")
                h_fail = workspace_bkg.data(f"Minv_bkg_fail_({i}|{j})_{cat}")
                n_pass, n_fail = h_pass.sumEntries(), h_fail.sumEntries()
                
                if divide_for_data:
                    h_data_pass = ws_data.data(f"Minv_data_pass_({i}|{j})")
                    h_data_fail = ws_data.data(f"Minv_data_fail_({i}|{j})")
                    print(n_pass, h_data_pass.sumEntries())
                    n_pass, n_fail = n_pass/h_data_pass.sumEntries(), n_fail/h_data_fail.sumEntries()

                histos_pass[cat].SetBinContent(i, j, n_pass)
                histos_pass[cat].SetBinError(i, j, (n_pass**0.5)*(lumi_scales[cat]**0.5))  # The histograms are weighted with the lumi-scale, so the error on the total number of events is the sqrt of the resulting number (= sqrt(lumi_scale * num_init)) multiplied with sqrt(lumi_scale), to obtain the correct value of lumi_scale*sqrt(num_init)

                histos_fail[cat].SetBinContent(i, j, n_fail)
                histos_fail[cat].SetBinError(i, j, (n_fail**0.5)*(lumi_scales[cat]**0.5))
        
        c = ROOT.TCanvas()
        c.Divide(1,2)

        c.cd(1)
        ROOT.gStyle.SetOptStat("en")
        histos_pass[cat].Draw("COLZ")
        c.cd(2)
        ROOT.gStyle.SetOptStat("en")
        histos_fail[cat].Draw("COLZ")

        if divide_for_data:
            c.SaveAs(f"figs/backgrounds/relative_distributions/{cat}_bkg_distribution.pdf")
        else:
            c.SaveAs(f"figs/backgrounds/distributions/{cat}_bkg_distribution.pdf")

        if save_root_file:
            rootfile_distrib.cd()
            histos_pass[cat].Write()
            histos_fail[cat].Write()
    
    if save_root_file:
        rootfile_distrib.Close()

    
    return histos_pass, histos_fail
       

if __name__ == "__main__":

    # bkg_categories = ["WW", "WZ", "ZZ", "TTSemileptonic", "Ztautau"]
    bkg_categories = ["ZZ"]
    binning_pt = array('d', [24., 26., 28., 30., 32., 34., 36., 38., 40., 42., 44., 47., 50., 55., 60., 65.])
    binning_eta = array('d', [round(-2.4 + i*0.1, 2) for i in range(49)])

    new_binning_pt = array('d', [24., 28., 32., 36., 40., 44., 47., 50., 55., 60., 65.])
    scale_bin_eta = 4
    new_binning_eta = array('d', [round(-2.4 + i*0.1*scale_bin_eta, 2) for i in range(int(48/scale_bin_eta)+1)])
    

    lumi_scales = get_bkg_lumi_scale("iso", bkg_categories)  # To be multiplied w/ the number of events
    
    ws = generate_bkg_datasets("iso", bkg_categories, lumi_scales, binning_pt, binning_eta, 
                               merge_bins=True, new_binning_pt=binning_pt, new_binning_eta=new_binning_eta)

    
    # file = ROOT.TFile("root_files/ws/ws_backgrounds_iso.root")
    # ws = file.Get("w")
    
    '''
    draw_bkg_distrib_2d(ws, bkg_categories, lumi_scales, binning_pt, binning_eta, divide_for_data=False,  
                        file_ws_data="results/benchmark_iso/ws_iso_indep_benchmark.root",
                        save_root_file=True)
    '''

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

    res = results_manager("indep")
    res.Open("results/benchmark_iso/results_iso_indep_benchmark.pkl")
    res_dict = res.dictionary()

    print("*********")
    # print(h_pass_gen.GetBinContent(1,1), h_pass_gen.GetBinError(1,1))
    print(res_dict["1,1"]["pars_pass"].find("nbkg_pass_(1|1)"))
    print(res_dict["1,1"]["pars_fail"].find("nbkg_fail_(1|1)"))















