"""
"""
import ROOT
from array import array
from copy import copy
from utilities.base_lib import binning, bin_dictionary, eval_efficiency, sumw2_error
from utilities.res_tools import results_manager, init_results_histos, fill_res_histograms
from utilities.dataset_utils import import_pdf_library

NBINS = 21

eff_min = 0.895
rel_err_eff_min = 1e-4
rel_err_eff_max = 7e-3
sf_min = 0.99
sf_max = 1.03


delta_min = -5e-2
delta_error_min = -5e-3
pull_min = -5
rm1_min = -3e-2
ratio_error_min = -1

res_var_dict = {
    "efficiency" : {
        "title" : "Efficiency", 
        "array" : array("d", [round(eff_min + (1-eff_min)*(i/NBINS), 4) for i in range(NBINS+1)])},
    "rel_err_efficiency" : {
        "title" : "Relative error on efficiency", 
        "array" : array("d", [round(rel_err_eff_min + (rel_err_eff_max-rel_err_eff_min)*(i/NBINS), 6) for i in range(NBINS+1)])},
    "sf" : {
        "title" : "Scale Factor", 
        "array" : array("d", [round(sf_min + (sf_max-sf_min)*(i/NBINS), 5) for i in range(NBINS+1)])},
}
resCmp_var_dict = {
    "delta" : {
        "title" : "Delta efficiency", 
        "array" : array("d", [round(delta_min + (-2*delta_min/NBINS)*i, 6) for i in range(NBINS+1)])},
    "delta_error" : {
        "title" : "Delta error", 
        "array" : array("d", [round(delta_error_min + (-2*delta_error_min/NBINS)*i, 6) for i in range(NBINS+1)])},
    "pull" : {
        "title" : "Pull", 
        "array" : array("d", [round(pull_min + (-2*pull_min/NBINS)*i, 6) for i in range(NBINS+1)])},
    "pull_ref" : {
        "title" : "Pull (referred to bmark stat. error)", 
        "array" : array("d", [round(pull_min + (-2*pull_min/NBINS)*i, 6) for i in range(NBINS+1)])},
    "rm1" : {
        "title" : "Relative bias", 
        "array" : array("d", [round(rm1_min + (-2*rm1_min/NBINS)*i, 6) for i in range(NBINS+1)])},
    "ratio_error" : {
        "title" : "Ratio error minus 1", 
        "array" : array("d", [round(ratio_error_min + (-2*ratio_error_min/NBINS)*i, 6) for i in range(NBINS+1)])}
}



###############################################################################

def save_eff_results(ws_name, type_analysis, binning_pt, binning_eta):
    """
    """
    file_in = ROOT.TFile(ws_name, "READ")
    ws = file_in.Get("w")

    bins_pt, bins_eta = binning(binning_pt), binning(binning_eta)

    bin_dict = bin_dictionary(binning_pt, binning_eta)
    
    results = results_manager(type_analysis, binning_pt, binning_eta, import_ws=ws)

    histos ={}
    
    for res_key in ["efficiency", "rel_err_efficiency", "sf"]:
        histos.update(init_results_histos(res_key, res_var_dict[res_key]["title"], 
                                          res_var_dict[res_key]["array"], bins_pt, bins_eta))

    
    for bin_key, [_, bin_pt, bin_eta] in bin_dict.items():

        eff, d_eff = results.getEff(bin_key)
        eff_mc, d_eff_mc = eval_efficiency(ws.data(f"Minv_mc_pass_{bin_key}").sumEntries(), 
                                           ws.data(f"Minv_mc_fail_{bin_key}").sumEntries(),
                                           sumw2_error(ws.data(f"Minv_mc_pass_{bin_key}")),
                                           sumw2_error(ws.data(f"Minv_mc_fail_{bin_key}")))

        histos["efficiency"].Fill(eff)
        histos["efficiency_2d"].SetBinContent(bin_pt, bin_eta, eff)
        histos["rel_err_efficiency"].Fill(d_eff/eff)
        histos["rel_err_efficiency_2d"].SetBinContent(bin_pt, bin_eta, d_eff/eff)
        histos["sf"].Fill(eff/eff_mc)
        histos["sf_2d"].SetBinContent(bin_pt, bin_eta, eff/eff_mc)
    

    file_out = ROOT.TFile(ws_name.replace("ws", "res"), "RECREATE")
    file_out.cd()
    [histo.Write() for histo in histos.values()]
    file_out.Close()

###############################################################################


def compare_efficiency(ws_txt_bmark_filename, ws_new_filename, binning_pt, binning_eta, res_list,
                       eval_nobkg_effect=False, postfix_name=""):
    """
    Compare the efficiencies and their error between two results files. The first file is 
    considered as benchmark
    """

    if ".root" in ws_txt_bmark_filename:
        file_bmark = ROOT.TFile(ws_txt_bmark_filename, "READ")
        ws_bmark = file_bmark.Get("w")
        print(type(ws_bmark))
        res_benchmark = results_manager("indep", binning_pt, binning_eta, import_ws=ws_bmark)
    else:
        with open(ws_txt_bmark_filename, "r") as file_bmark:
            row_list = file_bmark.readlines()
        print(type(row_list))
        res_benchmark = results_manager("indep", "pt_tracking", "eta", import_txt=row_list)
    
    file_new = ROOT.TFile(ws_new_filename, "READ")
    ws_new = file_new.Get("w")


    res_new = results_manager("indep", binning_pt, binning_eta, import_ws=ws_new)

    bins_pt, bins_eta = binning(binning_pt), binning(binning_eta)
    bin_dict = bin_dictionary(binning_pt, binning_eta)

    bin_dict_original = copy(bin_dict)

    if eval_nobkg_effect is True:
        if t_an == "indep":
            for b_key in bin_dict_original.keys():
                pars_pass, pars_fail = [res_new.getPars(flag, b_key) for flag in ["pass", "fail"]]
                if pars_pass.getSize() == 5 or pars_fail.getSize() == 5:
                    bin_dict.pop(b_key)
        elif t_an == "sim":
            for b_key in bin_dict_original.keys():
                pars_sim = res_new.getPars("sim", b_key)
                if pars_sim.getSize() == 10:
                    bin_dict.pop(b_key)

        NBINS = len(bin_dict.keys())
    
 
    histos = {}
    [histos.update(init_results_histos(res_key, resCmp_var_dict[res_key]["title"], 
                                       resCmp_var_dict[res_key]["array"], bins_pt, bins_eta)) 
     for res_key in res_list]
    
    print(type(res_benchmark), type(res_new))

    fill_res_histograms(res_benchmark, res_new, histos, bin_dict) 

    '''
    histos_copy = histos.copy()
    for hist_key in histos_copy.keys():
        if "2d" in hist_key:
            histos[hist_key].Sumw2()
            histo_pt = histos[hist_key].ProjectionX(hist_key.replace("2d", "pt"), 1, len(bins_eta)-1)
            histo_pt.Scale(1/(len(bins_eta)-1))
            histo_eta = histos[hist_key].ProjectionY(hist_key.replace("2d", "eta"), 1, len(bins_pt)-1)
            histo_eta.Scale(1/(len(bins_pt)-1))
            histos.update({histo_pt.GetName() : histo_pt, histo_eta.GetName() : histo_eta})
    '''

    add_flag = "cmp" if ".root" in ws_txt_bmark_filename else "cmp-egm"

    if postfix_name!="": add_flag += f"-{postfix_name}"

    if eval_nobkg_effect: add_flag += "_noBkgFits"

    file_out = ROOT.TFile(ws_new_filename.replace("ws", f"res_{add_flag}"), "RECREATE")
    file_out.cd()
    [histo.Write() for histo in histos.values()]
    file_out.Close()

    print(file_out.GetName())


###############################################################################

def compare_eff_pseudodata(ws_filename, binning_pt, binning_eta, res_list, isolate_effect=""):
    """
    """

    file_pseudodata = ROOT.TFile.Open(ws_filename)
    ws = file_pseudodata.Get("w")

    bins_pt, bins_eta = binning(binning_pt), binning(binning_eta)
    bin_dict = bin_dictionary(binning_pt, binning_eta)

    nbins_pt, nbins_eta = len(bins_pt)-1, len(bins_eta)-1

    histos = {}
    [histos.update(init_results_histos(res_key, resCmp_var_dict[res_key]["title"], 
                                       resCmp_var_dict[res_key]["array"], bins_pt, bins_eta)) 
     for res_key in res_list]

    results = results_manager("indep", binning_pt, binning_eta, import_ws=ws)
    
    for bin_key, [_, bin_pt, bin_eta] in bin_dict.items():
        
        eff, d_eff = results.getEff(bin_key)

        if isolate_effect=="pass":
            nfail = ws.data(f"Minv_mc_fail_{bin_key}").sumEntries()
            d_nfail = sumw2_error(ws.data(f"Minv_mc_fail_{bin_key}"))   
            pars_fail = ws.obj(f"results_pass_{bin_key}").floatParsFinal()
            npass, d_npass = pars_fail.find(f"nsig_pass_{bin_key}").getVal(), pars_fail.find(f"nsig_pass_{bin_key}").getError()
            eff, d_eff = eval_efficiency(npass, nfail, d_npass, d_nfail)

        elif isolate_effect=="fail":
            npass = ws.data(f"Minv_mc_pass_{bin_key}").sumEntries()
            d_npass = sumw2_error(ws.data(f"Minv_mc_pass_{bin_key}"))

            pars_fail = ws.obj(f"results_fail_{bin_key}").floatParsFinal()
            nfail, d_nfail = pars_fail.find(f"nsig_fail_{bin_key}").getVal(), pars_fail.find(f"nsig_fail_{bin_key}").getError()

            eff, d_eff = eval_efficiency(npass, nfail, d_npass, d_nfail)
        else:
            pass

        eff_mc, d_eff_mc = eval_efficiency(ws.data(f"Minv_mc_pass_{bin_key}").sumEntries(), 
                                           ws.data(f"Minv_mc_fail_{bin_key}").sumEntries(),
                                           sumw2_error(ws.data(f"Minv_mc_pass_{bin_key}")),
                                           sumw2_error(ws.data(f"Minv_mc_fail_{bin_key}")))
        
        '''
        if bin_key in ["[55.0to65.0][-0.3to-0.2]", "[55.0to65.0][0.2to0.3]",
                       "[55.0to65.0][-0.2to-0.1]", "[55.0to65.0][0.5to0.6]", "[55.0to65.0][0.7to0.8]"]:
            
            histos["delta_2d"].SetBinContent(bin_pt, bin_eta, -999)
            histos["delta_error_2d"].SetBinContent(bin_pt, bin_eta, -999)
            histos["pull_2d"].SetBinContent(bin_pt, bin_eta, -999)
            histos["rm1_2d"].SetBinContent(bin_pt, bin_eta, -999)
            histos["ratio_error_2d"].SetBinContent(bin_pt, bin_eta, -999)
            continue
        '''

        if "delta" in histos.keys():
            histos["delta"].Fill(eff-eff_mc)
            histos["delta_2d"].SetBinContent(bin_pt, bin_eta, eff-eff_mc)
        if "delta_error" in histos.keys():
            histos["delta_error"].Fill(d_eff-d_eff_mc)
            histos["delta_error_2d"].SetBinContent(bin_pt, bin_eta, d_eff-d_eff_mc)
        if "pull" in histos.keys():
            histos["pull"].Fill((eff-eff_mc)/d_eff)
            histos["pull_2d"].SetBinContent(bin_pt, bin_eta, (eff-eff_mc)/d_eff)
        if "rm1" in histos.keys():
            histos["rm1"].Fill((eff/eff_mc)-1)
            histos["rm1_2d"].SetBinContent(bin_pt, bin_eta, (eff/eff_mc)-1)
        if "ratio_error" in histos.keys():
            histos["ratio_error"].Fill((d_eff/d_eff_mc)-1)
            histos["ratio_error_2d"].SetBinContent(bin_pt, bin_eta, (d_eff/d_eff_mc)-1)

        print(bin_key, (eff-eff_mc)/d_eff)


    resfile = ROOT.TFile(ws_filename.replace("ws", f"res_cmpMC"), "RECREATE")
    resfile.cd()
    [histo.Write() for histo in histos.values()]             
    resfile.Close()

###############################################################################

def eval_minos(ws_hesse_filename, ws_filename, binning_pt, binning_eta):
    """
    """
    file = ROOT.TFile.Open(ws_filename)
    ws = file.Get("w")

    file_hesse = ROOT.TFile.Open(ws_hesse_filename)
    ws_hesse = file_hesse.Get("w")

    bins_pt, bins_eta = binning(binning_pt), binning(binning_eta)
    bin_dict = bin_dictionary(binning_pt, binning_eta)

    asymm_binning = array("d", [-0.5 + (i/75.) for i in range(76)])
    ratio_errors_binning = array("d", [-0.3 + (0.6*i/75.) for i in range(76)])

    histos ={}
    
    
    histos.update(init_results_histos("asymm_minos", "MINOS error asymmetry",
                                      asymm_binning, bins_pt, bins_eta))
    histos.update(init_results_histos("ratio_errors", "Ratio errors",
                                      ratio_errors_binning, bins_pt, bins_eta))
    

    for bin_key, [_, bin_pt, bin_eta] in bin_dict.items():

        res, res_hesse = ws.obj(f"results_sim_{bin_key}"), ws_hesse.obj(f"results_sim_{bin_key}")
        
        pars, pars_hesse = res.floatParsFinal(), res_hesse.floatParsFinal()

        eff, eff_hesse = pars.find(f"efficiency_{bin_key}"), pars_hesse.find(f"efficiency_{bin_key}")

        asymm_minos = (eff.getErrorHi() - abs(eff.getErrorLo()))/(eff.getErrorHi() + abs(eff.getErrorLo()))

        minos_width = (eff.getErrorHi() + abs(eff.getErrorLo()))
        hesse_width = eff_hesse.getError()*2

        print(minos_width/hesse_width)
        histos["asymm_minos"].Fill(asymm_minos)
        histos["asymm_minos_2d"].SetBinContent(bin_pt, bin_eta, asymm_minos)
        histos["ratio_errors"].Fill((minos_width/hesse_width)-1)
        histos["ratio_errors_2d"].SetBinContent(bin_pt, bin_eta, (minos_width/hesse_width)-1)

    
    file_out = ROOT.TFile(ws_filename.replace("ws", "asym"), "RECREATE")
    file_out.cd()
    [histo.Write() for histo in histos.values()]
    file_out.Close()



###############################################################################
###############################################################################

if __name__ == '__main__':

    resCmp_list = resCmp_var_dict.keys()
    print(resCmp_list)

    import_pdf_library("RooCMSShape")

    base_folder = "/scratch/rforti/tnp_efficiencies_results/tracking"

    ws_filename = base_folder+"/pseudodata/ws_tracking_pseudodata.root"

    # benchmark_res = base_folder+"/legacy_fit_onlyFailSA_allMC/ws_tracking.root"


    # save_eff_results(benchmark_res, "indep", "pt_tracking", "eta")
    # compare_efficiency(benchmark_res, ws_filename, "pt_tracking", "eta", resCmp_list, postfix_name="legacy-onlyFailSA")

    # eval_minos("results/iso_sim/ws_iso_sim.root", "results/iso_sim_minos/ws_iso_sim_minos_eff.root", "pt", "eta")

    compare_eff_pseudodata(ws_filename, "pt_tracking", "eta", resCmp_list, isolate_effect="")