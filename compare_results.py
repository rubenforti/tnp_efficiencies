"""
"""
import ROOT
from array import array
from utilities.base_library import binning, bin_dictionary, eval_efficiency, sumw2_error
from utilities.results_utils import results_manager, init_results_histos


bins_delta_eff = array("d", [round(-8e-4 + (1.6e-3/75)*i, 5) for i in range(75+1)])
bins_delta_deff = array("d", [round(-1e-4 + (2e-4/75)*i, 6) for i in range(75+1)])
bins_pull = array("d", [round(-1 + (2/75)*i, 5) for i in range(75+1)])
bins_ratio = array("d", [round(0.995 + (0.01/50)*i, 5) for i in range(50+1)])


def compare_efficiency(ws_bmark, ws_new, binning_pt, binning_eta, file_output, auxiliary_res={}):
    """
    Compare the efficiencies and their error between two results files. The first file is 
    considered as benchmark
    """
    
    for t_an in ["indep", "sim"]:
        an_bmark = t_an if t_an in ws_bmark.GetName() else ""
        an_new = t_an if t_an in ws_new.GetName() else ""

    
    res_benchmark = results_manager("indep", binning_pt, binning_eta, import_ws=ws_bmark)
    res_new = results_manager("indep", binning_pt, binning_eta, import_ws=ws_new)

    '''
    if auxiliary_res["filename"]!="":
        aux_file = ROOT.TFile(auxiliary_res.pop("filename"), "READ")
        aux_ws = aux_file.Get("w")
        for key_missing in auxiliary_res.keys():
            res_pass = aux_ws.obj(f"results_pass_{auxiliary_res[key_missing]}")
            res_fail = aux_ws.obj(f"results_fail_{auxiliary_res[key_missing]}")
            res_new.add_result({"pass":res_pass, "fail":res_fail}, key_missing)
    '''

    bins_pt, bins_eta = binning(binning_pt), binning(binning_eta)
    nbins_pt, nbins_eta = len(bins_pt)-1, len(bins_eta)-1

    bin_dict = bin_dictionary(binning_pt, binning_eta)
    
    bmark_dict = res_benchmark.dictionary()
            
 
    histos = {}
    histos.update(init_results_histos("delta_eff", "Delta efficiency", 
                                      bins_delta_eff, bins_pt, bins_eta))
    histos.update(init_results_histos("delta_error_eff", "Delta error on efficiency",
                                      bins_delta_deff, bins_pt, bins_eta))
    histos.update(init_results_histos("pull_eff", "Pull efficiency", 
                                      bins_pull, bins_pt, bins_eta))

    cnt_mergedpt = 0

    for bin_key in bmark_dict.keys():

        _, bin_pt, bin_eta = bin_dict[bin_key]

       # Bin transformation needed in case the bins are merged
        if type(bin_eta) is list:
            bin_eta = int(1+(nbins_eta*(bin_eta[0]-1)/48.))
        if type(bin_pt) is list:
            bin_pt_list = bin_pt
            bin_pt = int(bin_pt_list[0] - cnt_mergedpt)
            cnt_mergedpt += bin_pt_list[-1]-bin_pt_list[0] if bin_eta==nbins_eta else 0

        eff_1, deff_1 = res_benchmark.getEff(bin_key)
        eff_2, deff_2 = res_new.getEff(bin_key)

        delta_eff = eff_2-eff_1
        delta_deff = deff_2-deff_1

        histos["delta_eff"].Fill(delta_eff)
        histos["delta_eff_2d"].SetBinContent(bin_pt, bin_eta, delta_eff)

        histos["delta_error_eff"].Fill(delta_deff)
        histos["delta_error_eff_2d"].SetBinContent(bin_pt, bin_eta, delta_deff)

        histos["pull_eff"].Fill(delta_eff/deff_2)
        histos["pull_eff_2d"].SetBinContent(bin_pt, bin_eta, delta_eff/deff_1)

    file_out = ROOT.TFile(file_output, "RECREATE")
    file_out.cd()
    [histo.Write() for histo in histos.values()]
    file_out.Close()

###############################################################################

def compare_with_benchmark(ws, type_analysis, ref_txt, file_output):

    with open(ref_txt, "r") as file:
        row_list = file.readlines()
    
    bins_pt, bins_eta = binning("pt"), binning("eta")
    bin_dict = bin_dictionary("pt","eta")
            
    histos = {}
    histos.update(init_results_histos("delta_eff", "Delta efficiency", 
                                      bins_delta_eff, bins_pt, bins_eta))
    histos.update(init_results_histos("delta_error_eff", "Delta error on efficiency",
                                      bins_delta_deff, bins_pt, bins_eta))
    histos.update(init_results_histos("pull_eff", "Pull efficiency", 
                                      bins_pull, bins_pt, bins_eta))

    idx_list = 3

    print(row_list[0])
    print(row_list[1])
    print(row_list[2])
    print(row_list[3])

    results = results_manager(type_analysis, "pt", "eta", import_ws=ws)

    for bin_key in bin_dict.keys():
        
        _, bin_pt, bin_eta = bin_dict[bin_key]

        eff, deff = results.getEff(bin_key)
        elements = row_list[idx_list].split('\t')

        histos["delta_eff"].Fill(eff-float(elements[4]))
        histos["delta_eff_2d"].SetBinContent(bin_pt, bin_eta, eff-float(elements[4]))

        histos["delta_error_eff"].Fill(deff-float(elements[5]))
        histos["delta_error_eff_2d"].SetBinContent(bin_pt, bin_eta, deff-float(elements[5]))

        histos["pull_eff"].Fill((eff-float(elements[4]))/float(elements[5]))
        histos["pull_eff_2d"].SetBinContent(bin_pt, bin_eta, (eff-float(elements[4]))/float(elements[5]))

        idx_list += 1

    file_out = ROOT.TFile(file_output, "RECREATE")
    file_out.cd()
    [histo.Write() for histo in histos.values()]
    file_out.Close()
    '''
    c0 = ROOT.TCanvas("", "", 1200, 900)
    c0.Divide(2)
    c0.cd(1)
    ROOT.gStyle.SetOptStat("menr")
    ROOT.gPad.SetLogy()
    h_delta_eff.Draw()
    c0.cd(2)
    ROOT.gPad.SetLogy()
    h_delta_deff.Draw()
    c0.SaveAs("figs/delta_eff.pdf")

    c1 = ROOT.TCanvas("delta_eff", "delta_eff", 1200, 900)
    c1.cd()
    ROOT.gStyle.SetOptStat("en")
    ROOT.gPad.SetRightMargin(0.15)
    ROOT.gStyle.SetPalette(57)
    h2d_delta_eff.Draw("COLZ")
    h2d_delta_eff.SetContour(25)
    c1.SaveAs("figs/delta_eff_2d.pdf")

    c2 = ROOT.TCanvas("pull_delta_eff", "pull_delta_eff", 1200, 900)
    c2.cd()
    ROOT.gStyle.SetOptStat("en")
    ROOT.gPad.SetRightMargin(0.15)
    ROOT.gStyle.SetPalette(57)
    h2d_pull.Draw("COLZ")
    h2d_pull.SetContour(25)
    c2.SaveAs("figs/pull_delta_eff_2d.pdf")
    '''


   

###############################################################################

def compare_eff_pseudodata(ws, binning_pt, binning_eta, file_output):
    """
    """

    bins_pt, bins_eta = binning(binning_pt), binning(binning_eta)
    bin_dict = bin_dictionary(binning_pt, binning_eta)

    nbins_pt, nbins_eta = len(bins_pt)-1, len(bins_eta)-1

    #serializzare
    histos = {}

    histos.update(init_results_histos("pull", "Pull efficiency", bins_pull, bins_pt, bins_eta))
    histos.update(init_results_histos("ratio", "Ratio efficiency", bins_ratio, bins_pt, bins_eta))
    histos.update(init_results_histos("ratio_errors", "Ratio errors on efficiency", 
                                      bins_ratio, bins_pt, bins_eta))

    results = results_manager("indep", binning_pt, binning_eta, import_ws=ws)

    cnt_mergedpt = 0
    
    for bin_key in bin_dict.keys():

        _, bin_pt, bin_eta = bin_dict[bin_key]

       # Bin transformation needed in case the bins are merged
        if type(bin_eta) is list:
            bin_eta = int(1+(nbins_eta*(bin_eta[0]-1)/48.))
        if type(bin_pt) is list:
            bin_pt_list = bin_pt
            bin_pt = int(bin_pt_list[0] - cnt_mergedpt)
            cnt_mergedpt += bin_pt_list[-1]-bin_pt_list[0] if bin_eta==nbins_eta else 0
        
        eff, d_eff = results.getEff(bin_key)

        eff_mc, d_eff_mc = eval_efficiency(ws.data(f"Minv_mc_pass_{bin_key}").sumEntries(), 
                                           ws.data(f"Minv_mc_fail_{bin_key}").sumEntries(),
                                           sumw2_error(ws.data(f"Minv_mc_pass_{bin_key}")),
                                           sumw2_error(ws.data(f"Minv_mc_fail_{bin_key}")))

        pull = (eff-eff_mc)/d_eff_mc
        ratio = eff / eff_mc
        ratio_errors  = d_eff/d_eff_mc

        histos["ratio"].Fill(ratio)
        histos["ratio_2d"].SetBinContent(bin_pt, bin_eta, ratio)
        histos["pull"].Fill(pull)
        histos["pull_2d"].SetBinContent(bin_pt, bin_eta, pull)
        histos["ratio_errors"].Fill(ratio_errors)
        histos["ratio_errors_2d"].SetBinContent(bin_pt, bin_eta, ratio_errors)

    resfile = ROOT.TFile(file_output, "RECREATE")
    resfile.cd()
    [histo.Write() for histo in histos.values()]             
    resfile.Close()

###############################################################################
###############################################################################

if __name__ == '__main__':

    
    file = ROOT.TFile.Open("results/benchmark_iso/ws_iso_indep_benchmark.root")
    ws = file.Get("w")
    compare_with_benchmark(ws, "indep", "results/benchmark_iso/old_results.txt", 
                           "results/benchmark_iso/hres_cmp_iso_indep_benchmark.root")
    

    '''
    file_bmark = ROOT.TFile.Open("results/bmark_iso_2gev/ws_iso_indep_bmark_2gev.root")
    ws_bmark = file_bmark.Get("w")

    file_test = ROOT.TFile.Open("root_files/ws_iso_indep_mcbkg_2gev.root")
    ws_test = file_test.Get("w")

    compare_efficiency(ws_bmark, ws_test, "pt", "eta", "bkg_results/hres_mcbkg_2gev_cmp.root")
    '''

    '''
    file_pseudodata = ROOT.TFile.Open("bkg_results/ws_bkg_pseudodata.root")
    ws_pseudodata = file_pseudodata.Get("w")
    compare_eff_pseudodata(ws_pseudodata, "pt", "eta_8bins", "bkg_results/hres_pseudodata_cmp.root")
    '''

    

