"""
"""

import ROOT
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

    roohist_data = ROOT.RooDataHist(f"Minv_data_{flag}_({bin_pt},{bin_eta})",
                                    f"Minv_data_{flag}({bin_pt},{bin_eta})",
                                    [axis], Import=th1_data)
    roohist_mc = ROOT.RooDataHist(f"Minv_mc_{flag}_({bin_pt},{bin_eta})",
                                  f"Minv_mc_{flag}_({bin_pt},{bin_eta})",
                                  [axis], Import=th1_mc)

    return (roohist_data, roohist_mc)


def ws_init(type_eff, bins_pt, bins_eta, bins_mass):
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
    f_data = ROOT.TFile(f"root_files/tnp_{type_eff}_data.root")
    f_mc = ROOT.TFile(f"root_files/tnp_{type_eff}_mc.root")

    for i in range(1, len(bins_pt)):
        for j in range(1, len(bins_eta)):
            histos_pass = get_roohist((f_data, f_mc), x, i, j, 'pass')
            histos_fail = get_roohist((f_data, f_mc), x, i, j, 'fail')
            # Datasets are written in this order: data_pass, mc_pass, data_fail, mc_fail
            w.Import(histos_pass[0])
            w.Import(histos_pass[1])
            w.Import(histos_fail[0])
            w.Import(histos_fail[1])

    w.writeToFile(f"root_files/{type_eff}_workspace.root")

    return w


def ws_init_std_pdf(workspace):
    """
    """

    axis = workspace["x"]

    # Gaussian smearing
    # -----------------
    mean = ROOT.RooRealVar("mean", "mean", 0, -2, 2)
    sigma = ROOT.RooRealVar("sigma", "sigma", 0.5, 0.001, 2)
    gaussian = ROOT.RooGaussian(
        "gaus_smearing", "gaussian smearing", axis, mean, sigma)
    workspace.Import(gaussian)

    # Exponential bkg
    # ---------------
    tau = ROOT.RooRealVar("tau", "tau", -10, 0)
    expo_bkg = ROOT.RooExponential(
        "expo_bkg", "exponential background", axis, tau)
    workspace.Import(expo_bkg)


if __name__ == '__main__':

    binning_pt = array('d', [24., 26., 28., 30., 32., 34.,
                       36., 38., 40., 42., 44., 47., 50., 55., 60., 65.])

    binning_eta = array('d', [round(-2.4 + i*0.1, 2) for i in range(49)])

    binning_mass = array('d', [50 + i for i in range(81)])

    w = ws_init('iso', [1, 1], binning_eta, binning_mass)
    ws_init_std_pdf(w)
    w.writeToFile(f"root_files/iso_workspace.root")

    w.Print()

    '''
    f = ROOT.TFile("root_files/iso_workspace.root")
    w = f.Get("w")
    '''