"""
"""

import ROOT
import sys
from array import array


def get_roohist(file, type_set, flag, axis, bin_pt, bin_eta, global_scale=-1.):
    """
    Returns two RooDataHists (data and mc) of the variable TP_mass, in a
    (pt, eta) bin, selected from the TH3 given as input.
    """

    if len(bin_pt) == 1:
        bin_pt.append(bin_pt[0])
    
    if len(bin_eta) == 1:
        bin_eta.append(bin_eta[0])

    if type_set == "data":
        type_suffix = "RunGtoH"
    elif type_set == "mc" or type_set == "bkg":
        type_suffix = "DY_postVFP"
    else:
        print("ERROR: invalid category type given!")
        sys.exit()

    print(f"{flag}_mu_{type_suffix}")

    
    histo3d = file.Get(f"{flag}_mu_{type_suffix}")

    th1_histo = histo3d.ProjectionX(f"Histo_data_{flag}", bin_pt[0], bin_pt[1], bin_eta[0], bin_eta[1], "e")  # Option "e" is specified to calculate the bin errors in the new histogram for generic selection of bin_pt and bin_eta. Without it, it all works well ONLY IF the projection is done on one single bin of (pt, eta)
    
    if global_scale > 0:
        th1_histo.Scale(global_scale)

    print(f"Num TH1 entries = {th1_histo.Integral()}")

    roohistogram = ROOT.RooDataHist(f"Minv_{type_set}_{flag}_({bin_pt}|{bin_eta})",
                                    f"Minv_{type_set}_{flag}_({bin_pt}|{bin_eta})",
                                    ROOT.RooArgList(axis), th1_histo)

    return roohistogram


def ws_init(file_data, file_mc, type_eff, type_analysis, bins_pt, bins_eta, bins_mass):
    """
    Initializes a RooWorkspace with the dataset corresponding to a given
    efficiency step. The objects stored are RooDataHist of TP inv. mass in
    every (pt,eta) bin
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



if __name__ == '__main__':

    type_eff = ("sa", "global", "ID", "iso", "trigger", "veto")
    t = type_eff[3]

    types_analysis = ["indep", "sim"]
    an = types_analysis[1]

    binning_pt = array('d', [24., 26., 28., 30., 32., 34., 36., 38., 40., 42., 44., 47., 50., 55., 60., 65.])

    binning_eta = array('d', [round(-2.4 + i*0.1, 2) for i in range(49)])

    binning_mass = array('d', [60 + i for i in range(61)])

    filename_data = "/scratchnvme/wmass/Steve_root_files/Standard_SF_files/tnp_iso_data_vertexWeights1_oscharge1.root"
    filename_mc = "/scratchnvme/wmass/Steve_root_files/Standard_SF_files/tnp_iso_mc_vertexWeights1_oscharge1.root"

    f_data = ROOT.TFile(filename_data)
    f_mc = ROOT.TFile(filename_mc)
    print(type(f_data))
    print(type(f_mc))

    w = ws_init(f_data, f_mc, t, an, [1,1], [1,1], binning_mass)

    w.writeToFile(f"root_files/ws/ws_{t}_{an}_new.root")


