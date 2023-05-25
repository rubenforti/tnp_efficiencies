"""
"""

import ROOT
from workspace_config import get_roohist
from results_utils import results_manager
from array import array


def initialize_bkg_h2d(categories, binning_pt, binning_eta):

    histos_pass = []
    histos_fail = []

    if "WW" in categories:
        h2d_WW_pass = ROOT.TH2D(
            "WW_pass", "WW_pass", len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta)
        h2d_WW_fail = ROOT.TH2D(
            "WW_fail", "WW_fail", len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta)
        histos_pass.append(h2d_WW_pass)
        histos_fail.append(h2d_WW_fail)

    if "WZ" in categories:
        h2d_WZ_pass = ROOT.TH2D(
            "WZ_pass", "WZ_pass", len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta)
        h2d_WZ_fail = ROOT.TH2D(
            "WZ_fail", "WZ_fail", len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta)
        histos_pass.append(h2d_WZ_pass)
        histos_fail.append(h2d_WZ_fail)

    if "ZZ" in categories:
        h2d_ZZ_pass = ROOT.TH2D(
            "ZZ_pass", "ZZ_pass", len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta)
        h2d_ZZ_fail = ROOT.TH2D(
            "ZZ_fail", "ZZ_fail", len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta)
        histos_pass.append(h2d_ZZ_pass)
        histos_fail.append(h2d_ZZ_fail)

    if "TTSemileptonic" in categories:
        h2d_TTsemileptonic_pass = ROOT.TH2D("TTsemileptonic_pass", "TTsemileptonic_pass", 
                                            len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta)   
        h2d_TTsemileptonic_fail = ROOT.TH2D("TTsemileptonic_fail", "TTsemileptonic_fail", 
                                            len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta) 
        histos_pass.append(h2d_TTsemileptonic_pass)
        histos_fail.append(h2d_TTsemileptonic_fail)
    
    if "Ztautau" in categories:
        h2d_Ztautau_pass = ROOT.TH2D("Ztautau_pass", "Ztautau_pass", 
                                     len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta)
        h2d_Ztautau_fail = ROOT.TH2D("Ztautau_fail", "Ztautau_fail", 
                                     len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta)
        histos_pass.append(h2d_Ztautau_pass)
        histos_fail.append(h2d_Ztautau_fail)
    
    if "SameCharge" in categories:
        h2d_samecharge_pass = ROOT.TH2D("SameCharge_pass", "SameCharge_pass", 
                                     len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta)
        h2d_samecharge_fail = ROOT.TH2D("SameCharge_fail", "SameCharge_fail", 
                                     len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta)
        histos_pass.append(h2d_samecharge_pass)
        histos_fail.append(h2d_samecharge_fail)



    return histos_fail, histos_pass


BR_TAUToMU = 0.1739
BR_TAUToE = 0.1782
Z_TAU_TO_LEP_RATIO = (1.-(1. - BR_TAUToMU - BR_TAUToE)**2)
xsec_ZmmPostVFP = 2001.9

cross_sections_bkg = {
    # Unit = pb
    "WW" : 12.6,
    "WZ" : 5.4341,
    "ZZ" : 0.60+5.1,
    "TTSemileptonic" : 366.34,
    "Ztautau" : xsec_ZmmPostVFP*Z_TAU_TO_LEP_RATIO
}
    

def get_bkg_lumi_scale(type_eff, categories):
    """
    "Lumi scale" defined as alpha that satisifies alpha*lumi_bkg = lumi_data
    """

    lumi_data = 16.8  # fb^-1

    lumi_scale = {}

    for cat in categories:
        file = ROOT.TFile(f"/scratchnvme/rajarshi/Bkg_TNP_3D_Histograms/OS/tnp_{type_eff}_{cat}_vertexWeights1_oscharge1.root")
        
        wsum_histo = file.Get("weightSum")
        num_init = wsum_histo.Integral()
        xsection = cross_sections_bkg[cat]*1000
        lumi_bkg = num_init/xsection

        scale = lumi_data/lumi_bkg


        lumi_scale.update({cat : scale})
    
    return lumi_scale



def plot_bkg_on_data(ws, ws_bkg, ):
    pass



def generate_bkg_roohists(type_eff, categories, binning_pt, binning_eta, make_2d_histos=True):

    if make_2d_histos:
        histos_fail, histos_pass = initialize_bkg_h2d(categories, binning_pt, binning_eta)

    axis = ROOT.RooRealVar("x_bkg_pass", "x", 50, 130)  # If we are going to convolve the bkg distributions, we will need a specific axis for every histogram in every bin


    ws = ROOT.RooWorkspace("w")
    ws.Import(axis)

    lumi_scales = get_bkg_lumi_scale(type_eff, categories)  # To be multiplied w/ the number of events
    
    for idx in range(len(categories)):

        cat = categories[idx]

        
        file = ROOT.TFile(f"/scratchnvme/rajarshi/Bkg_TNP_3D_Histograms/OS/tnp_{type_eff}_{categories[idx]}_vertexWeights1_oscharge1.root")

        print(file.GetName())
        for i in range(1, len(binning_pt)):
            for j in range(1, len(binning_eta)):
                
                h_pass = get_roohist(file, "bkg", "pass", axis, i, j, global_scale=lumi_scales[cat])
                print("AAAAA")
                h_pass.SetName(f"{h_pass.GetName()}_bkg_{cat}")
                n_pass = h_pass.sumEntries()
                histos_pass[idx].SetBinContent(i, j, n_pass)
                histos_pass[idx].SetBinError(i, j, n_pass*(lumi_scales[cat]**0.5))
                ws.Import(h_pass)

                h_fail = get_roohist(file, "bkg", "fail", axis, i, j, global_scale=lumi_scales[cat])
                h_fail.SetName(f"{h_fail.GetName()}_bkg_{cat}")
                n_fail = h_fail.sumEntries()
                histos_fail[idx].SetBinContent(i, j, n_fail)
                histos_fail[idx].SetBinError(i, j, n_fail*(lumi_scales[cat]**0.5))
                ws.Import(h_fail)

        file.Close()
            
        c = ROOT.TCanvas()
        c.Divide(1,2)

        c.cd(1)
        ROOT.gStyle.SetOptStat("en")
        histos_pass[idx].Draw("COLZ")
        c.cd(2)
        ROOT.gStyle.SetOptStat("en")
        histos_fail[idx].Draw("COLZ")

        c.SaveAs(f"figs/backgrounds/{categories[idx]}_bkg_distribution.pdf")
    
    ws.writeToFile("root_files/ws/ws_backgrounds.root")

    return histos_pass, histos_fail




if __name__ == "__main__":

    categories = ["WW", "WZ", "ZZ", "TTSemileptonic", "Ztautau"]
    binning_pt = array('d', [24., 26., 28., 30., 32., 34., 36., 38., 40., 42., 44., 47., 50., 55., 60., 65.])
    binning_eta = array('d', [round(-2.4 + i*0.1, 2) for i in range(49)])


    h_pass, h_fail = generate_bkg_roohists("iso", categories, binning_pt, binning_eta)

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

    res = results_manager("sim")
    res.Open("results/benchmark_iso/results_iso_indep_benchmark.pkl")
    res_dict = res.dictionary()

    print("*********")
    print(h_pass_gen.GetBinContent(1,1), h_pass_gen.GetBinError(1,1))
    print(res_dict["1,1"]["pars_pass"].Print("v"))















