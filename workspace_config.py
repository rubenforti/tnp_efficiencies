"""
"""

import ROOT
from array import array


def get_roohist(histos, axis, bin_pt, bin_eta, flag):
    """
    Returns a RooDataHist of the variable TP_mass, in a (pt, eta) bin, selected
    from the TH3 given as input.
    """

    # Option "e" has to be activated? Not clear how "errors are computed"

    if flag == 'data':
        h_pass = histos.pass_mu_RunGtoH
        h_fail = histos.fail_mu_RunGtoH
    elif flag == 'mc':
        h_pass = histos.pass_mu_DY_postVFP
        h_fail = histos.fail_mu_DY_postVFP

    th1_pass = h_pass.ProjectionX(
        f"Histo_{flag}_pass", bin_pt, bin_pt, bin_eta, bin_eta)
    th1_fail = h_fail.ProjectionX(
        f"Histo_{flag}_fail", bin_pt, bin_pt, bin_eta, bin_eta)

    print(f"Num TH1 entries = {th1_pass.GetEntries()}")

    roohist_pass = ROOT.RooDataHist(f"Minv_{flag}_pass_({bin_pt},{bin_eta})",
                                    f"Minv_{flag}_pass_({bin_pt},{bin_eta})",
                                    [axis], Import=th1_pass)
    roohist_fail = ROOT.RooDataHist(f"Minv_{flag}_fail_({bin_pt},{bin_eta})",
                                    f"Minv_{flag}_fail_({bin_pt},{bin_eta})",
                                    [axis], Import=th1_fail)

    return (roohist_pass, roohist_fail)


def ws_init(type_eff, bins_pt, bins_eta):

    w = ROOT.RooWorkspace("w")

    x = ROOT.RooRealVar("x", "TP M_{inv}", 50, 130, unit="GeV/c^2")
    w.Import(x)

    # Import of the 3D histograms
    f_data = ROOT.TFile(f"root_files/tnp_{type_eff}_data.root")
    f_mc = ROOT.TFile(f"root_files/tnp_{type_eff}_mc.root")

    for i in range(1, len(bins_pt)):
        for j in range(1, len(bins_eta)):
            histos_data = get_roohist(f_data, x, i, j, 'data')
            histos_mc = get_roohist(f_mc, x, i, j, 'mc')
            w.Import(histos_data[0])
            w.Import(histos_data[1])
            w.Import(histos_mc[0])
            w.Import(histos_mc[1])

    mean = ROOT.RooRealVar("mean", "mean", 0, -2, 2)
    sigma = ROOT.RooRealVar("sigma", "sigma", 0.5, 0.001, 2)
    smearing = ROOT.RooGaussian("smearing", "smearing", x, mean, sigma)
    w.Import(smearing)

    tau = ROOT.RooRealVar("tau", "tau", -10, 0)
    expo_bkg = ROOT.RooExponential("expo", "expo", x, tau)
    w.Import(expo_bkg)

    w.writeToFile(f"root_files/{type_eff}_workspace.root")

    '''
    nevts = ((nev_fail_data, nev_pass_data), (nev_fail_mc, nev_pass_mc))
    histos_data = (h_fail_data, h_pass_data)
    histos_mc = (h_fail_mc, h_pass_mc)
    '''

    return w


if __name__ == '__main__':

    binning_pt = array('d', [24., 26., 28., 30., 32., 34.,
                       36., 38., 40., 42., 44., 47., 50., 55., 60., 65.])

    binning_eta = array('d', [round(-2.4 + i*0.1, 2) for i in range(49)])

    w = ws_init('iso', [1,1], [1,1])

    '''
    f = ROOT.TFile("root_files/iso_workspace.root")
    w = f.Get("w")
    '''

    w.Print()

    aaa = w["Minv_data_fail_(1,1)"]
    print(aaa.sumEntries())
    print(aaa.numEntries())
    '''
    for i in range(1, 16):
        for j in range(1, 49):
            w.removeSet(f"Minv_data_fail_({i},{j})")
            w.removeSet(f"Minv_mc_fail_({i},{j})")
    '''
    # w.Print()
