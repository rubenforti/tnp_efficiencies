"""
"""

import sys
import os
import ROOT
from utilities.base_library import binning, bin_dictionary, lumi_factors, get_idx_from_bounds


def import_pdf_library(*functions):
    """
    """
    current_path = os.path.dirname(__file__)
    
    for function in functions:
        header_incl = ' #include "libCpp/'+function+'.h"'
        sourcefile = os.path.join(current_path, 'libCpp', f'{function}.cc')
        # print(sourcefile)
        ctrl_head = ROOT.gInterpreter.Declare(header_incl)
        ctrl_source = ROOT.gSystem.CompileMacro(sourcefile, opt="ks")

        if ctrl_head is not True:
            print("ERROR in header loading")
            sys.exit()
        if ctrl_source != 1:
            print("ERROR in sourcefile compiling and loading")
            sys.exit()

###############################################################################

def get_roohist(file, type_set, flag, axis, bin_key, bin_pt, bin_eta, global_scale=-1.):
    """
    Returns a RooDataHists of the variable TP_invmass, in a single (pt, eta) bin
    """

    if type_set=="data":
        type_suffix = "RunGtoH"
    elif type_set=="mc" or type_set=="mc_w" or type_set=="bkg":
        type_suffix = "DY_postVFP"
    else:
        print("ERROR: invalid category type given!")
        sys.exit()

    print(f"{flag}_mu_{type_suffix}")

    histo3d = file.Get(f"{flag}_mu_{type_suffix}")


    if type(bin_pt) is int:
        bin_pt = [bin_pt, bin_pt]
    if type(bin_eta) is int:
        bin_eta = [bin_eta, bin_eta]    

    th1_histo = histo3d.ProjectionX(f"Histo_{type_set}_{flag}", bin_pt[0], bin_pt[1], bin_eta[0], bin_eta[1], "e")  # Option "e" is specified to calculate the bin errors in the new histogram for generic selection of bin_pt and bin_eta. Without it, it all works well ONLY IF the projection is done on one single bin of (pt, eta)
    
    if global_scale > 0:
        th1_histo.Scale(global_scale)

    print(f"TH1 integral = {th1_histo.Integral()}")
    if th1_histo.Integral() < 0:
        print("ERROR: negative entries in TH1")
        print(f"{bin_pt}, {bin_eta}")
        # sys.exit()

    numBins = axis.getBinning().numBins()
    th1_histo.Rebin(int(th1_histo.GetNbinsX()/numBins))

    
    roohisto = ROOT.RooDataHist(f"Minv_{type_set}_{flag}_{bin_key}", f"Minv_{type_set}_{flag}_{bin_key}",
                                ROOT.RooArgList(axis), th1_histo)
    
    return roohisto

###############################################################################

def ws_init(import_datasets, type_analysis, bins, bins_mass):
    """
    Initializes a RooWorkspace with the datasets corresponding to a given efficiency step. The objects 
    stored are RooDataHist of TP_invmass in every (pt,eta) bin. The datasets are of the type ("data, "mc",
    "bkg") specified in the import_sets dictionary that must contain, as values: the path of the file for
    "data" and "mc types; two other dictionaries in the the "bkg" case, containing the paths of the files 
    and the dictionary of the luminosity scales.
    """
    w = ROOT.RooWorkspace("w")
    x = ROOT.RooRealVar("x", "TP M_inv", bins_mass[0], bins_mass[-1], unit="GeV/c^2")
    w.Import(x)

    x_binning = ROOT.RooUniformBinning(bins_mass[0], bins_mass[-1], len(bins_mass)-1, "x_binning")

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

            _, bin_pt, bin_eta = bins[bin_key]

            if type_analysis == 'indep':
                x_pass = ROOT.RooRealVar(f"x_pass_{bin_key}", "TP M_inv", 
                                        bins_mass[0], bins_mass[-1], unit="GeV/c^2")
                x_pass.setBinning(x_binning)
                x_fail = ROOT.RooRealVar(f"x_fail_{bin_key}", "TP M_inv",
                                        bins_mass[0], bins_mass[-1], unit="GeV/c^2")
                x_fail.setBinning(x_binning)
                axis = (x_fail, x_pass)
                # w.Import(x_pass)
                # w.Import(x_fail)
            elif type_analysis == 'sim':
                x_sim = ROOT.RooRealVar(f"x_sim_{bin_key}", "TP M_inv",
                                        bins_mass[0], bins_mass[-1], unit="GeV/c^2")
                x_sim.setBinning(x_binning)
                axis = (x_sim, x_sim)
            else:
                print("****\nINVALID ANALYSIS TYPE\n****")
                sys.exit()
            
            if dataset_type == "data":
                histo_pass = get_roohist(file, dataset_type, "pass", axis[1], bin_key, bin_pt, bin_eta)
                histo_fail = get_roohist(file, dataset_type, "fail", axis[0], bin_key, bin_pt, bin_eta)
                w.Import(histo_pass)
                w.Import(histo_fail)
                global_counter += 2
            if dataset_type == "mc":
                histo_pass = get_roohist(file, dataset_type, "pass", axis[1], 
                                         bin_key, bin_pt, bin_eta, global_scale=sig_lumi_scale)
                histo_fail = get_roohist(file, dataset_type, "fail", axis[0], 
                                         bin_key, bin_pt, bin_eta, global_scale=sig_lumi_scale)
                w.Import(histo_pass)
                w.Import(histo_fail)
                global_counter += 2
            elif dataset_type == "bkg":
                for cat in bkg_lumi_scales:
                    histo_pass = get_roohist(file_set[cat], dataset_type, "pass", axis[1], 
                                             bin_key, bin_pt, bin_eta, global_scale=bkg_lumi_scales[cat])
                    histo_pass.SetName(f"{histo_pass.GetName()}_{cat}")
                    histo_fail = get_roohist(file_set[cat], dataset_type, "fail", axis[0], 
                                             bin_key, bin_pt, bin_eta, global_scale=bkg_lumi_scales[cat])
                    histo_fail.SetName(f"{histo_fail.GetName()}_{cat}")
                    w.Import(histo_pass)
                    w.Import(histo_fail)
                    global_counter += 2
            else:
                pass

    print(global_counter)
            
    return w

###############################################################################

def show_negweighted_bins(ws, bkg_categories, binning_pt, binning_eta):
    """
    """

    bins_pt, bins_eta = binning(binning_pt), binning(binning_eta)
    nbins_pt, nbins_eta = len(bins_pt)-1, len(bins_eta)-1

    bin_dict = bin_dictionary(binning_pt, binning_eta)

    negweight_single_bkg = ROOT.TH2D("negweight_single_bkg", "negweight_single_bkg", 
                                     nbins_pt, bins_pt, nbins_eta, bins_eta)
    negweight_total_bkg = ROOT.TH2D("negweight_total_bkg", "negweight_total_bkg", 
                                    nbins_pt, bins_pt, nbins_eta, bins_eta)

    cnt_mergedpt=0

    for bin_key in bin_dict.keys():

        cnt_negweight_single = 0

        _, bin_pt, bin_eta = bin_dict[bin_key]

        # Bin transformation needed in case the bins are merged
        if type(bin_eta) is list:
            bin_eta = int(1+(nbins_eta*(bin_eta[0]-1)/48.))
        if type(bin_pt) is list:
            bin_pt_list = bin_pt
            bin_pt = int(bin_pt_list[0] - cnt_mergedpt)
            cnt_mergedpt += bin_pt_list[-1]-bin_pt_list[0] if bin_eta==nbins_eta else 0

        n_ev_pass = 0
        n_ev_fail = 0
        for cat in bkg_categories:
            histo_pass = ws.data(f"Minv_bkg_pass_{bin_key}_{cat}")
            histo_fail = ws.data(f"Minv_bkg_fail_{bin_key}_{cat}")
            n_ev_pass += histo_pass.sumEntries()
            n_ev_fail += histo_fail.sumEntries()

            if histo_pass.sumEntries()<0 or histo_fail.sumEntries()<0:
                cnt_negweight_single += 1
        
        if cnt_negweight_single>0:
            negweight_single_bkg.SetBinContent(bin_pt, bin_eta, cnt_negweight_single)
                
        if n_ev_pass<0 or n_ev_fail<0:
            negweight_single_bkg.SetBinContent(bin_pt, bin_eta, 1)
    
    c = ROOT.TCanvas()
    c.Divide(1,2)
    c.cd(1)
    negweight_single_bkg.Draw("colz")
    c.cd(2)
    negweight_total_bkg.Draw("colz")
    c.SaveAs(f"prob_bins_pt{nbins_pt}bins_eta{nbins_eta}bins.png")

###############################################################################


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

    binning_mass = binning("mass_60_120")

    binnings_list = [# ["pt", "eta"], ["pt", "eta_24bins"], ["pt", "eta_16bins"], 
                     # ["pt", "eta_8bins"], ["pt_12bins", "eta"], ["pt_10bins", "eta"],
                     # ["pt_8bins", "eta"], ["pt_12bins", "eta_24bins"]]
                     ["pt_10bins", "eta_16bins"]]
    
    # bin_dict = bin_dictionary("pt_12bins", "eta_24bins")

    for binning_pt, binning_eta in binnings_list:
        bin_set = bin_dictionary(binning_pt, binning_eta)

        w = ws_init(import_dictionary, an, bin_set, binning_mass)
        
        show_negweighted_bins(w, bkg_categories, binning_pt, binning_eta)




    # w = ws_init_backgrounds(dirname_bkg, bkg_categories, t, lumi_scales, bin_set, binning_mass)
    # w.Print()

   




   #w.writeToFile(f"root_files/ws/ws_data_mc_bkg.root")