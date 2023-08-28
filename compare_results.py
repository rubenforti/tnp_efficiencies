"""
"""
import ROOT
from array import array
from utilities.base_library import binning, bin_dictionary, eval_efficiency, sumw2_error
from utilities.results_utils import results_manager, init_results_histos, fill_res_histograms

NBINS = 75

eff_min = 0.85
rel_err_eff_min = 1e-4
rel_err_eff_max = 0.01
sf_min = 0.99
sf_max = 1.01


delta_min = -8e-3
delta_error_min = -1e-3
pull_min = -3
rm1_min = -0.01
ratio_error_min = 0.5

res_array_dict = {
    "eff" : array("d", [round(eff_min + (1-eff_min)*(i/NBINS), 4) for i in range(NBINS+1)]),
    "rel_err_eff" : array("d", [round(rel_err_eff_min + (rel_err_eff_max-rel_err_eff_min)*(i/NBINS), 6) 
                                for i in range(NBINS+1)]),
    "sf" : array("d", [round(sf_min + (sf_max-sf_min)*(i/NBINS), 5) for i in range(NBINS+1)]),

}

cmpres_array_dict = {
    "delta" : array("d", [round(delta_min + (-2*delta_min/NBINS)*i, 6) for i in range(NBINS+1)]),
    "delta_error" : array("d", [round(delta_error_min + (-2*delta_error_min/NBINS)*i, 6) for i in range(NBINS+1)]),
    "pull" : array("d", [round(pull_min + (-2*pull_min/NBINS)*i, 6) for i in range(NBINS+1)]),
    "rm1" : array("d", [round(rm1_min + (-2*rm1_min/NBINS)*i, 6) for i in range(NBINS+1)]),
    "ratio_error" : array("d", [round(ratio_error_min + (2*(1-ratio_error_min)/NBINS)*i, 6) for i in range(NBINS+1)])
}

###############################################################################

def save_eff_results(ws_name, type_analysis, binning_pt, binning_eta):
    """
    """
    file_in = ROOT.TFile.Open(ws_name, "UPDATE")
    ws = file_in.Get("w")

    bins_pt, bins_eta = binning(binning_pt), binning(binning_eta)

    bin_dict = bin_dictionary(binning_pt, binning_eta)
    
    results = results_manager(type_analysis, binning_pt, binning_eta, import_ws=ws)


    histos ={}
    
    [histos.update(init_results_histos(res, res, res_array_dict[res], bins_pt, bins_eta)) 
     for res in ["efficiency", "eff_rel_error", "sf"]]

    
    for bin_key in bin_dict.keys():
        
        _, bin_pt, bin_eta = bin_dict[bin_key]

        eff, d_eff = results.getEff(bin_key)
        eff_mc, d_eff_mc = eval_efficiency(ws.data(f"Minv_mc_pass_{bin_key}").sumEntries(), 
                                           ws.data(f"Minv_mc_fail_{bin_key}").sumEntries(),
                                           sumw2_error(ws.data(f"Minv_mc_pass_{bin_key}")),
                                           sumw2_error(ws.data(f"Minv_mc_fail_{bin_key}")))

        histos["efficiency"].Fill(eff)
        histos["efficiency_2d"].SetBinContent(bin_pt, bin_eta, eff)
        histos["eff_rel_error"].Fill(d_eff/eff)
        histos["eff_rel_error_2d"].SetBinContent(bin_pt, bin_eta, d_eff/eff)
        histos["sf"].Fill(eff/eff_mc)
        histos["sf_2d"].SetBinContent(bin_pt, bin_eta, eff/eff_mc)
    

    file_out = ROOT.TFile(ws_name.replace("ws", "hres"), "RECREATE")
    file_out.cd()
    [histo.Write() for histo in histos.values()]
    file_out.Close()

###############################################################################

def compare_efficiency(ws_txt_bmark_filename, ws_new_filename, binning_pt, binning_eta, res_list):
    """
    Compare the efficiencies and their error between two results files. The first file is 
    considered as benchmark
    """

    if ".root" in ws_txt_bmark_filename:
        file_bmark = ROOT.TFile(ws_txt_bmark_filename, "READ")
        ws_bmark = file_bmark.Get("w")
        print(type(ws_bmark))
        for t_an in ["indep", "sim"]: 
            if t_an in ws_txt_bmark_filename:
                res_benchmark = results_manager(t_an, binning_pt, binning_eta, import_ws=ws_bmark)
    else:
        with open(ws_txt_bmark_filename, "r") as file_bmark:
            row_list = file_bmark.readlines()
        print(type(row_list))
        res_benchmark = results_manager("indep", "pt", "eta", import_txt=row_list)
    
    file_new = ROOT.TFile(ws_new_filename, "READ")
    ws_new = file_new.Get("w")

    print("indep" in ws_new_filename)
    for t_an in ["indep", "sim"]:
        if t_an in ws_new_filename:
            res_new = results_manager(t_an, binning_pt, binning_eta, import_ws=ws_new)

    bins_pt, bins_eta = binning(binning_pt), binning(binning_eta)
    bin_dict = bin_dictionary(binning_pt, binning_eta)

    dict_new = res_new._dict_results
    print(type(dict_new))
    print("N keys", len(dict_new.keys()))
 
    histos = {}
    [histos.update(init_results_histos(res_key, res_key, cmpres_array_dict[res_key], bins_pt, bins_eta)) 
     for res_key in res_list]

    fill_res_histograms(res_benchmark, res_new, histos, bin_dict, len(bins_eta)-1) 

    histos_copy = histos.copy()
    for hist_key in histos_copy.keys():
        if "2d" in hist_key:
            histos[hist_key].Sumw2()
            histo_pt = histos[hist_key].ProjectionX(hist_key.replace("2d", "pt"), 1, len(bins_eta)-1)
            histo_pt.Scale(1/(len(bins_eta)-1))
            histo_eta = histos[hist_key].ProjectionY(hist_key.replace("2d", "eta"), 1, len(bins_pt)-1)
            histo_eta.Scale(1/(len(bins_pt)-1))
            histos.update({histo_pt.GetName() : histo_pt, histo_eta.GetName() : histo_eta})
       

    file_out = ROOT.TFile(ws_new_filename.replace("ws", "hres_cmp"), "RECREATE")
    file_out.cd()
    [histo.Write() for histo in histos.values()]
    file_out.Close()

###############################################################################

def compare_eff_pseudodata(ws_filename, binning_pt, binning_eta, res_list, file_output):
    """
    """

    file_pseudodata = ROOT.TFile.Open(ws_filename)
    ws = file_pseudodata.Get("w")

    bins_pt, bins_eta = binning(binning_pt), binning(binning_eta)
    bin_dict = bin_dictionary(binning_pt, binning_eta)

    nbins_pt, nbins_eta = len(bins_pt)-1, len(bins_eta)-1

    histos = {}
    [histos.update(init_results_histos(res_key, res_key, cmpres_array_dict[res_key], bins_pt, bins_eta)) 
     for res_key in res_list]

    results = results_manager("indep", binning_pt, binning_eta, import_ws=ws)

    cnt_mergedpt = 0

    print(histos.keys())
    
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

        if "delta" in histos.keys():
            histos["delta"].Fill(eff-eff_mc)
            histos["delta_2d"].SetBinContent(bin_pt, bin_eta, eff-eff_mc)
        if "delta_error" in histos.keys():
            histos["delta_error"].Fill(d_eff-d_eff_mc)
            histos["delta_error_2d"].SetBinContent(bin_pt, bin_eta, d_eff-d_eff_mc)
        if "pull" in histos.keys():
            histos["pull"].Fill((eff-eff_mc)/d_eff_mc)
            histos["pull_2d"].SetBinContent(bin_pt, bin_eta, (eff-eff_mc)/d_eff_mc)
        if "rm1" in histos.keys():
            histos["rm1"].Fill((eff/eff_mc)-1)
            histos["rm1_2d"].SetBinContent(bin_pt, bin_eta, (eff/eff_mc)-1)
        if "ratio_error" in histos.keys():
            histos["ratio_error"].Fill(d_eff/d_eff_mc)
            histos["ratio_error_2d"].SetBinContent(bin_pt, bin_eta, d_eff/d_eff_mc)

    resfile = ROOT.TFile(file_output, "RECREATE")
    resfile.cd()
    [histo.Write() for histo in histos.values()]             
    resfile.Close()

###############################################################################
###############################################################################

if __name__ == '__main__':

    res_list = cmpres_array_dict.keys()

    '''
    ws_results = [["results/benchmark_trig_plus/ws_triggerplus_indep_bmark.root", "results/benchmark_trig_plus/old_results.txt"],
                  ["results/benchmark_trig_minus/ws_triggerminus_indep_bmark.root", "results/benchmark_trig_minus/old_results.txt"],]

    for ws_new, ws_benchmark in ws_results:
        compare_efficiency(ws_benchmark, ws_new, "pt", "eta", res_list)
    '''
    
    ws_filename = "results/pseudodata_trig_minus/ws_triggerminus_pseudodata.root"
    
    compare_eff_pseudodata(ws_filename, "pt", "eta", res_list, "results/pseudodata_trig_minus/hres_cmp_pseudodata.root")

    

