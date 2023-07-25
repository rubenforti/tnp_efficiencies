"""
"""

import sys
import os
import ROOT
from utilities.base_library import binning, bin_dictionary, lumi_factors


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
        #sys.exit()

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

def show_empty_bins(ws, bkg_categories, binning, nbins_pt, nbins_eta):
    """
    """

    prob_bins_histo = ROOT.TH2D("prob_bins_histo", "prob_bins_histo", 
                                nbins_pt, 1, nbins_pt, nbins_eta, 1, nbins_eta)

    for bin_key in binning:
        _, bin_pt, bin_eta = binning[bin_key]

        n_ev_pass = 0
        n_ev_fail = 0
        for cat in bkg_categories:
            histo_pass = ws.data(f"Minv_bkg_pass_{bin_key}_{cat}")
            histo_fail = ws.data(f"Minv_bkg_fail_{bin_key}_{cat}")
            n_ev_pass += histo_pass.sumEntries()
            n_ev_fail += histo_fail.sumEntries()
        
        if n_ev_pass<1 or n_ev_fail<1:
            prob_bins_histo.Fill(bin_pt[0], bin_eta[0])
    
    c = ROOT.TCanvas()
    c.cd()
    prob_bins_histo.Draw("colz")
    c.SaveAs("../prob_bins_Ztautau.png")
        
###############################################################################     
###############################################################################


if __name__ == '__main__':

    bkg_categories = ["WW", "WZ", "ZZ", "TTSemileptonic", "Ztautau"]
    # bkg_categories = ["Ztautau"]

    type_eff = ("sa", "global", "ID", "iso", "trigger", "veto")
    t = type_eff[3]

    types_analysis = ["indep", "sim"]
    an = types_analysis[0]


    lumi_scales = lumi_factors(t, bkg_categories)
    # lumi_scales = {"WW":1, "WZ":1, "ZZ":1, "TTSemileptonic":1, "Ztautau":1}


    filename_data = "/scratchnvme/wmass/Steve_root_files/Standard_SF_files/tnp_iso_data_vertexWeights1_oscharge1.root"
    filename_mc = "/scratchnvme/wmass/Steve_root_files/Standard_SF_files/tnp_iso_mc_vertexWeights1_oscharge1.root"
    dirname_bkg = "/scratchnvme/rajarshi/Bkg_TNP_3D_Histograms/OS"


    bkg_filepaths = {}
    [bkg_filepaths.update({cat : 
        f"{dirname_bkg}/tnp_{t}_{cat}_vertexWeights1_oscharge1.root"}) for cat in bkg_categories]



    import_dictionary = {
        "data" : filename_data,
        "mc" : filename_mc,
        "bkg" : {
            "filepaths" : bkg_filepaths,
            "lumi_scales" : lumi_scales
        }
    }


    binning_mass = binning("mass_60_120")

    bin_set = bin_dictionary("pt_6bins", "eta_4bins")

    w = ws_init(import_dictionary, an, bin_set, binning_mass)

    # w = ws_init_backgrounds(dirname_bkg, bkg_categories, t, lumi_scales, bin_set, binning_mass)
    # w.Print()

    # show_problematic_bins(w, bkg_categories, bin_set, 15, 48)




    w.writeToFile(f"root_files/ws/ws_data_mc_bkg.root")