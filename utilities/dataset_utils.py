"""
"""

import sys
import os
import ROOT
from utilities.base_library import binning, bin_dictionary, lumi_factors, get_idx_from_bounds, bin_global_idx_dict


bkg_datasets_sel = {
    "WW" : 1,
    "WZ" : 1,
    "ZZ" : 1,
    "TTFullyleptonic" : 1,
    "TTSemileptonic" : 0,
    "Ztautau" : 1,
    "SameCharge" : 1,
    "WJets" : 0
}


def import_pdf_library(*functions):
    """
    """
    current_path = os.path.dirname(__file__)
    import_path = os.path.join(current_path, '..', 'libCpp')
    
    for function in functions:
        ctrl_head = ROOT.gInterpreter.Declare(f' #include "{import_path}/{function}.h"')
        ctrl_source = ROOT.gSystem.CompileMacro(f"{import_path}/{function}.cc", opt="ks")

        if ctrl_head is not True:
            print("ERROR in header loading")
            sys.exit()
        if ctrl_source != 1:
            print("ERROR in sourcefile compiling and loading")
            sys.exit()


###############################################################################

def import_totbkg_hist(ws, bin_key, bkg_categories):
    """
    """
    for flag in ["pass", "fail"]:
        axis = ws.var(f"x_{flag}_{bin_key}")
        bkg_total_histo = ROOT.RooDataHist(f"Minv_bkg_{flag}_{bin_key}_total", "bkg_total_histo", 
                                       ROOT.RooArgSet(axis), "plot_binning")

        for bkg_cat in bkg_categories:
            bkg_dataset = ws.data(f"Minv_bkg_{flag}_{bin_key}_{bkg_cat}")
            if bkg_datasets_sel[bkg_cat] == 1: bkg_total_histo.add(bkg_dataset)

        ws.Import(bkg_total_histo)


###############################################################################

def get_roohist(file, type_set, flag, axis, bin_key, bin_pt, bin_eta, global_scale=-1.):
    """
    Returns a RooDataHists of the variable TP_invmass, in a single (pt, eta) bin
    """
    if type_set=="data" or type_set=="data_SC":
        type_suffix = "RunGtoH"
    elif type_set=="mc" or type_set=="mc_w" or type_set=="bkg":
        type_suffix = "DY_postVFP"
    else:
        print("ERROR: invalid category type given!")
        sys.exit()

    histo3d = file.Get(f"{flag}_mu_{type_suffix}")

    if type(bin_pt) is int:  bin_pt = [bin_pt, bin_pt]
    if type(bin_eta) is int: bin_eta = [bin_eta, bin_eta]    

    th1_histo = histo3d.ProjectionX(f"Histo_{type_set}_{flag}", bin_pt[0], bin_pt[1], bin_eta[0], bin_eta[1], "e")  # Option "e" is specified to calculate the bin errors in the new histogram for generic selection of bin_pt and bin_eta. Without it, it all works well ONLY IF the projection is done on one single bin of (pt, eta)
    
    # print("Th1 entries before scaling: ", th1_histo.Integral())
    # print("Under/overflow: ", th1_histo.GetBinContent(0), th1_histo.GetBinContent(th1_histo.GetNbinsX()+1))

    if global_scale > 0:
        th1_histo.Scale(global_scale)

    print(f"\nTH1 integral = {th1_histo.Integral()}")
    if th1_histo.Integral() < 0:
        print("ERROR: negative entries in TH1")
        print(f"{bin_pt}, {bin_eta}")
        # sys.exit()
    
    # print(type_set, global_scale)

    numBins = axis.getBinning().numBins()
    th1_histo.Rebin(int(th1_histo.GetNbinsX()/numBins))

    if type_set=="data_SC": type_set = "bkg"
    
    roohisto = ROOT.RooDataHist(f"Minv_{type_set}_{flag}_{bin_key}", f"Minv_{type_set}_{flag}_{bin_key}",
                                ROOT.RooArgList(axis), th1_histo)
    

    print(bin_pt, bin_eta)

    print(axis.getBinning().numBins())
    if axis.getBinning().numBins() == th1_histo.GetNbinsX():
        print(f"TH1 binning: {th1_histo.GetNbinsX()}")
        for i in range(th1_histo.GetNbinsX()):
            if roohisto.weight(i) != th1_histo.GetBinContent(i+1):
                print("ERROR")
                print(f"{bin_pt}, {bin_eta}")
                sys.exit()
        print(f"RooDataHist integral = {roohisto.sumEntries()}\n")

    return roohisto

###############################################################################

def ws_init(import_datasets, type_analysis, binning_pt, binning_eta, binning_mass, 
            import_existing_ws=False, existing_ws_filename="", altBinning_bkg=False):
    """
    Initializes a RooWorkspace with the datasets corresponding to a given efficiency step. The objects 
    stored are RooDataHist of TP_invmass in every (pt,eta) bin. The datasets are of the type ("data, "mc",
    "bkg") specified in the import_sets dictionary that must contain, as values: the path of the file for
    "data" and "mc types; two other dictionaries in the the "bkg" case, containing the paths of the files 
    and the dictionary of the luminosity scales.
    """

    if altBinning_bkg is False:
        bins = bin_dictionary(binning_pt, binning_eta)
    else:
        bin_idx_dict = bin_global_idx_dict(binning_pt, binning_eta)
        bins = bin_dictionary()

    bins_mass = binning(binning_mass)

    if import_existing_ws is True:
        ws_file = ROOT.TFile(existing_ws_filename)
        ws = ws_file.Get("w")
    else:
        ws = ROOT.RooWorkspace("w")
        x = ROOT.RooRealVar("x", "TP M_inv", bins_mass[0], bins_mass[-1], unit="GeV/c^2")
        x_binning = ROOT.RooUniformBinning(bins_mass[0], bins_mass[-1], len(bins_mass)-1, "x_binning")
        ws.Import(x)

    global_counter = 0


    for dataset_type in import_datasets:       

        if dataset_type=="data":
            file = ROOT.TFile(import_datasets[dataset_type], "READ")
        elif dataset_type=="mc":
            file = ROOT.TFile(import_datasets["mc"]["filename"], "READ")
            sig_lumi_scale = import_datasets["mc"]["lumi_scale"]
        elif dataset_type=="bkg":
            file_set={}
            for cat in import_datasets[dataset_type]["filenames"]:
                file_set.update({cat : ROOT.TFile(import_datasets["bkg"]["filenames"][cat], "READ")})
                print(file_set[cat].GetName())
            bkg_lumi_scales = import_datasets["bkg"]["lumi_scales"]
        else:
            print("****\nINVALID DATASET TYPE\n****")
            sys.exit()
        

        for bin_key in bins:

            gl_idx, bin_pt, bin_eta = bins[bin_key]

            if import_existing_ws is True:
                if type_analysis == "indep":
                    x_pass = ws.var(f"x_pass_{bin_key}")
                    x_fail = ws.var(f"x_fail_{bin_key}")
                    axis = (x_fail, x_pass)
                elif type_analysis == "sim":
                    x_sim = ws.var(f"x_sim_{bin_key}")
                    axis = (x_sim, x_sim)
            else:
                if type_analysis == 'indep':
                    x_pass = ROOT.RooRealVar(f"x_pass_{bin_key}", "TP M_inv", 
                                            bins_mass[0], bins_mass[-1], unit="GeV/c^2")
                    x_pass.setBinning(x_binning)
                    x_fail = ROOT.RooRealVar(f"x_fail_{bin_key}", "TP M_inv",
                                            bins_mass[0], bins_mass[-1], unit="GeV/c^2")
                    x_fail.setBinning(x_binning)
                    axis = (x_fail, x_pass)
                elif type_analysis == 'sim':
                    x_sim = ROOT.RooRealVar(f"x_sim_{bin_key}", "TP M_inv",
                                            bins_mass[0], bins_mass[-1], unit="GeV/c^2")
                    x_sim.setBinning(x_binning)
                    axis = (x_sim, x_sim)

            if dataset_type == "data":
                histo_pass = get_roohist(file, dataset_type, "pass", axis[1], bin_key, bin_pt, bin_eta)
                histo_fail = get_roohist(file, dataset_type, "fail", axis[0], bin_key, bin_pt, bin_eta)
                ws.Import(histo_pass)
                ws.Import(histo_fail)
                global_counter += 2

            if dataset_type == "mc":
                print(sig_lumi_scale)
                histo_pass = get_roohist(file, dataset_type, "pass", axis[1], 
                                         bin_key, bin_pt, bin_eta, global_scale=sig_lumi_scale)
                histo_fail = get_roohist(file, dataset_type, "fail", axis[0], 
                                         bin_key, bin_pt, bin_eta, global_scale=sig_lumi_scale)
                ws.Import(histo_pass)
                ws.Import(histo_fail)
                global_counter += 2

            elif dataset_type == "bkg":
                if altBinning_bkg is True:
                    gl_idx_key = str(gl_idx)
                    bin_pt, bin_eta = bin_idx_dict[gl_idx_key]
                for cat in bkg_lumi_scales:
                    dset_type = "data_SC" if cat == "SameCharge" else "bkg"
                    histo_pass = get_roohist(file_set[cat], dset_type, "pass", axis[1], 
                                            bin_key, bin_pt, bin_eta, global_scale=bkg_lumi_scales[cat])
                    histo_pass.SetName(f"{histo_pass.GetName()}_{cat}")
                    histo_fail = get_roohist(file_set[cat], dset_type, "fail", axis[0], 
                                            bin_key, bin_pt, bin_eta, global_scale=bkg_lumi_scales[cat])
                    histo_fail.SetName(f"{histo_fail.GetName()}_{cat}")
                    ws.Import(histo_pass)
                    ws.Import(histo_fail)
                    global_counter += 2
            else:
                pass

    print(global_counter)
            
    return ws


###############################################################################


'''
# NOT USED ANYWHERE. Probably created to add the bkg datasets in merged bins to a workspace with data in
# reco ones, but evidently it couldn't perform well. This feature is now implemented as an accessory path
# in ws_init.

def extend_merged_datasets(ws, merged_bin_key, bkg_categories):
    """
    """
    bin_dict = bin_dictionary("pt", "eta")

    bkg_datasets = {}
    for flag in ["pass", "fail"]:
        [bkg_datasets.update({f"{cat}_{flag}" : ws.data(f"Minv_bkg_{flag}_{merged_bin_key}_{cat}")}) for cat in bkg_categories]
        bkg_datasets.update({f"axis_{flag}" : ws.var(f"x_{flag}_{merged_bin_key}")})

    string_pt, string_eta = merged_bin_key.split("][")
    pt_min, pt_max = string_pt[1:].split("to")
    eta_min, eta_max = string_eta[:-1].split("to")

    pt_min, pt_max = float(pt_min), float(pt_max)
    eta_min, eta_max = float(eta_min), float(eta_max)

    global_merged_bins, _, _ = get_idx_from_bounds([pt_min, pt_max], [eta_min, eta_max])

    for bin_key in bin_dict.keys():

        gl_bin, _, _ = bin_dict[bin_key]

        if gl_bin in global_merged_bins:

            for flag in ["pass", "fail"]:
                new_axis = bkg_datasets[f"axis_{flag}"].Clone(f"x_{flag}_{bin_key}")
                ws.Import(new_axis)
                for cat in bkg_categories:
                    new_data = bkg_datasets[f"{cat}_{flag}"].Clone(f"Minv_bkg_{flag}_{bin_key}_{cat}")
                    ws.Import(new_data)
'''

      
###############################################################################     
###############################################################################


if __name__ == '__main__':

    bkg_categories = ["WW", "WZ", "ZZ", "TTSemileptonic", "Ztautau"]
    # bkg_categories = ["Ztautau"]

    type_eff = ("sa", "global", "ID", "iso", "trigger", "veto")
    t = type_eff[3]

    types_analysis = ["indep", "sim"]
    an = types_analysis[0]


    filename_data = "root_files/datasets/tnp_iso_data_vertexWeights1_oscharge1.root"
    filename_mc = "root_files/datasets/tnp_iso_mc_vertexWeights1_oscharge1.root"
    dirname_bkg = "root_files/datasets"
    
    
    bkg_filenames = {}
    [bkg_filenames.update({cat : 
        f"{dirname_bkg}/tnp_{t}_{cat}_vertexWeights1_oscharge1.root"}) for cat in bkg_categories]
    
    print(bkg_filenames)
    

    lumi_scales = lumi_factors(t, bkg_categories)

    print(lumi_scales)


    lumi_scale_signal = lumi_scales.pop("Zmumu")

    import_dictionary = {
        "data" : filename_data,
        "mc" : {
            "filename": filename_mc,
            "lumi_scale" : lumi_scale_signal
        },
        "bkg" : {
            "filenames" : bkg_filenames,
            "lumi_scales" : lumi_scales
        } 
    }