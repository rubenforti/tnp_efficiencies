"""
"""

import sys
import os
import ROOT
from array import array


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


def get_roohist(file, type_set, flag, axis, bin_key, bin_pt, bin_eta, global_scale=-1.):
    """
    Returns a RooDataHists of the variable TP_invmass, in a single (pt, eta) bin
    """

    if type_set == "data":
        type_suffix = "RunGtoH"
    elif type_set == "mc" or type_set == "bkg":
        type_suffix = "DY_postVFP"
    else:
        print("ERROR: invalid category type given!")
        sys.exit()

    print(f"{flag}_mu_{type_suffix}")

    histo3d = file.Get(f"{flag}_mu_{type_suffix}")

    if len(bin_pt) == 1:
        bin_pt.append(bin_pt[0])
    
    if len(bin_eta) == 1:
        bin_eta.append(bin_eta[0])


    th1_histo = histo3d.ProjectionX(f"Histo_data_{flag}", bin_pt[0], bin_pt[1], bin_eta[0], bin_eta[1], "e")  # Option "e" is specified to calculate the bin errors in the new histogram for generic selection of bin_pt and bin_eta. Without it, it all works well ONLY IF the projection is done on one single bin of (pt, eta)
    
    if global_scale > 0:
        th1_histo.Scale(global_scale)

    print(f"Num TH1 entries = {th1_histo.Integral()}")

    roohistogram = ROOT.RooDataHist(f"Minv_{type_set}_{flag}_{bin_key}", 
                                    f"Minv_{type_set}_{flag}_{bin_key}",
                                    ROOT.RooArgList(axis), th1_histo)
    
    return roohistogram


def ws_init(file_data, file_mc, type_eff, type_analysis, bins_pt, bins_eta, bins_mass):
    """
    Initializes a RooWorkspace with the dataset corresponding to a given efficiency step. The objects 
    stored are RooDataHist of TP_invmass in every (pt,eta) bin
    """
    
    # Import of the 3D histograms


    w = ROOT.RooWorkspace("w")

    x = ROOT.RooRealVar(
        "x", "TP M_inv", bins_mass[0], bins_mass[-1], unit="GeV/c^2")
    w.Import(x)

    for i in range(1, len(bins_pt)):
        for j in range(1, len(bins_eta)):

            if type_analysis == 'indep':
                x_pass = ROOT.RooRealVar(f"x_pass_({i}|{j})", "TP M_inv",
                                         bins_mass[0], bins_mass[-1], unit="GeV/c^2")

                x_fail = ROOT.RooRealVar(f"x_fail_({i}|{j})", "TP M_inv",
                                         bins_mass[0], bins_mass[-1], unit="GeV/c^2")

                axis = (x_fail, x_pass)
                w.Import(x_pass)
                w.Import(x_fail)

            elif type_analysis == 'sim':
                x_sim = ROOT.RooRealVar(f"x_sim_({i}|{j})", "TP M_inv",
                                        bins_mass[0], bins_mass[-1], unit="GeV/c^2")
                axis = (x_sim, x_sim)

            else:
                print("INVALID ANALYSIS TYPE")
                sys.exit()

            histo_data_pass = get_roohist(file_data, "data", "pass", axis[1], [i], [j])
            histo_mc_pass = get_roohist(file_mc, "mc", "pass", axis[1], [i], [j])

            histo_data_fail = get_roohist(file_data, "data", "fail", axis[0], [i], [j])
            histo_mc_fail = get_roohist(file_mc, "mc", "fail", axis[0], [i], [j])

            w.Import(histo_data_pass)
            w.Import(histo_mc_pass)
            w.Import(histo_data_fail)
            w.Import(histo_mc_fail)

    return w


def eval_efficiency(npass, nfail, sigma_npass, sigma_nfail):
    """
    """
    eff = npass/(npass+nfail)
    var1 = (nfail**2)*(sigma_npass**2)
    var2 = (npass**2)*(sigma_nfail**2)
    sigma_eff = ROOT.TMath.Sqrt(var1+var2)/((npass+nfail)**2)

    return eff, sigma_eff