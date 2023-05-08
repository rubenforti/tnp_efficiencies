"""
"""

import ROOT
import sys
from array import array


def get_roohist(histos, axis, bin_pt, bin_eta, flag):
    """
    Returns two RooDataHists (data and mc) of the variable TP_mass, in a
    (pt, eta) bin, selected from the TH3 given as input.
    """

    # Option "e" has to be activated? Not clear how "errors are computed"

    if flag == 'pass':
        h_data = histos[0].pass_mu_RunGtoH
        h_mc = histos[1].pass_mu_DY_postVFP
    elif flag == 'fail':
        h_data = histos[0].fail_mu_RunGtoH
        h_mc = histos[1].fail_mu_DY_postVFP

    th1_data = h_data.ProjectionX(
        f"Histo_data_{flag}", bin_pt, bin_pt, bin_eta, bin_eta)
    th1_mc = h_mc.ProjectionX(
        f"Histo_mc_{flag}", bin_pt, bin_pt, bin_eta, bin_eta)

    print(f"Num TH1 entries = {th1_data.Integral()}")

    roohist_data = ROOT.RooDataHist(f"Minv_data_{flag}_({bin_pt}|{bin_eta})",
                                    f"Minv_data_{flag}({bin_pt}|{bin_eta})",
                                    ROOT.RooArgList(axis), th1_data)
    roohist_mc = ROOT.RooDataHist(f"Minv_mc_{flag}_({bin_pt}|{bin_eta})",
                                  f"Minv_mc_{flag}_({bin_pt}|{bin_eta})",
                                  ROOT.RooArgList(axis), th1_mc)

    return (roohist_data, roohist_mc)


def ws_init(filenames, type_eff, type_analysis, bins_pt, bins_eta, bins_mass):
    """
    Initializes a RooWorkspace with the dataset corresponding to a given
    efficiency step. The objects stored are RooDataHist of TP inv. mass in
    every (pt,eta) bin
    """

    w = ROOT.RooWorkspace("w")

    x = ROOT.RooRealVar(
        "x", "TP M_inv", bins_mass[0], bins_mass[-1], unit="GeV/c^2")
    w.Import(x)

    # Import of the 3D histograms
    f_data = ROOT.TFile(filenames[0])
    f_mc = ROOT.TFile(filenames[1])

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

            histos_pass = get_roohist((f_data, f_mc), axis[1], i, j, 'pass')
            histos_fail = get_roohist((f_data, f_mc), axis[0], i, j, 'fail')

            # Datasets are written in this order: data_pass, mc_pass, data_fail, mc_fail
            w.Import(histos_pass[0])
            w.Import(histos_pass[1])
            w.Import(histos_fail[0])
            w.Import(histos_fail[1])


    return w



if __name__ == '__main__':

    type_eff = ("sa", "global", "ID", "iso", "trigger", "veto")
    t = type_eff[3]

    types_analysis = ["indep", "sim"]
    an = types_analysis[0]

    binning_pt = array('d', [24., 26., 28., 30., 32., 34., 36., 38., 40., 42., 44., 47., 50., 55., 60., 65.])

    binning_eta = array('d', [round(-2.4 + i*0.1, 2) for i in range(49)])

    binning_mass = array('d', [60 + i for i in range(61)])

    filename_data = "/scratchnvme/wmass/Steve_root_files/Standard_SF_files/tnp_iso_data_vertexWeights1_oscharge1.root"
    filename_mc = "/scratchnvme/wmass/Steve_root_files/Standard_SF_files/tnp_iso_mc_vertexWeights1_oscharge1.root"


    w = ws_init((filename_data, filename_mc), t, an, binning_pt, binning_eta, binning_mass)

    w.writeToFile(f"root_files/ws/ws_{t}_{an}.root")

