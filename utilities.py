"""
"""

import sys
import pickle
import ROOT


def import_pdf_library(*functions):
    """
    """
    for function in functions:
        header_incl = ' #include "libCpp/'+function+'.h"'
        sourcefile = 'libCpp/'+function+'.cc'

        ctrl_head = ROOT.gInterpreter.Declare(header_incl)
        ctrl_source = ROOT.gSystem.CompileMacro(sourcefile, opt="ks")

        if ctrl_head is not True:
            print("ERROR in header loading")
            sys.exit()
        if ctrl_source != 1:
            print("ERROR in sourcefile compiling and loading")
            sys.exit()


def profile_histo(histo_th3, axis, bin_pt, bin_eta, flag):
    """
    Returns a RooDataHist of the variable TP_mass, in a (pt, eta) bin, selected
    from the TH3 given as input.
    """

    if len(bin_pt) == 1:
        bin_pt.append(bin_pt[0])

    if len(bin_eta) == 1:
        bin_eta.append(bin_eta[0])

    # Option "e" has to be activated? Not clear how "errors are computed"
    histo_th1 = histo_th3.ProjectionX(
        f"Histo_mass_{flag}", bin_pt[0], bin_pt[1], bin_eta[0], bin_eta[1])

    nbins = histo_th1.GetNbinsX()

    xAxis = histo_th1.GetXaxis()
    x = ROOT.RooRealVar("x", "x", xAxis.GetXmin(),
                        xAxis.GetXmax(), unit="GeV/c^2")
    print(f"Num TH1 entries = {histo_th1.GetEntries()}")
    roohist = ROOT.RooDataHist(
        f"roohist_{flag}", f"roohist_{flag}", [x], Import=histo_th1)
    if roohist.numEntries() != nbins:
        print('**********************************')
        print('ERRORE NEL NUMERO DI BIN IMPORTATI')
        print('**********************************')

    return roohist, histo_th1.Integral(1, nbins)


def th3_checks(histo_th3):
    """
    """
    Nbinsx = histo_th3.GetNbinsX()
    Nbinsy = histo_th3.GetNbinsY()
    Nbinsz = histo_th3.GetNbinsZ()
    print(Nbinsx, Nbinsy, Nbinsz)
    histo_x = histo_th3.ProjectionX("histo_mass", 0, -1, 0, -1)
    print(f"Number of events (with uf/of): {histo_x.GetEntries()}")
    print(
        f"Number of effective events (with uf/of): {histo_x.GetEffectiveEntries()}")
    histo_x_o = histo_th3.ProjectionX("histo_mass_2", 1, 1, 1, 1)
    print(f"Number of events (without uf/of): {histo_x_o.GetEntries()}")
    print(
        f"Number of effective events (without uf/of): {histo_x_o.GetEffectiveEntries()}")
    print(f"Integral: {histo_x_o.Integral()}")
    print(f"Weights: {histo_x_o.GetSumOfWeights()}")


def import_Steve_histos(type_eff, bin_pt, bin_eta):

    x = ROOT.RooRealVar("x", "TP M_{inv}", 50, 130, unit="GeV/c^2")

    # Import of the 3D histograms
    f_data = ROOT.TFile(f"root_files/tnp_{type_eff}_data.root")
    f_mc = ROOT.TFile(f"root_files/tnp_{type_eff}_mc.root")

    h_pass_data, nev_pass_data = profile_histo(
        f_data.pass_mu_RunGtoH, x, bin_pt, bin_eta, 1)
    h_pass_data.SetNameTitle(
        f"Events {type_eff} pass", f"Events {type_eff} pass")
    h_fail_data, nev_fail_data = profile_histo(
        f_data.fail_mu_RunGtoH, x, bin_pt, bin_eta, 2)
    h_fail_data.SetNameTitle(
        f"Events {type_eff} fail", f"Events {type_eff} fail")
    h_pass_mc, nev_pass_mc = profile_histo(
        f_mc.pass_mu_DY_postVFP, x, bin_pt, bin_eta, 3)
    h_fail_mc, nev_fail_mc = profile_histo(
        f_mc.fail_mu_DY_postVFP, x, bin_pt, bin_eta, 4)

    nevts = ((nev_fail_data, nev_pass_data), (nev_fail_mc, nev_pass_mc))
    histos_data = (h_fail_data, h_pass_data)
    histos_mc = (h_fail_mc, h_pass_mc)

    return histos_data, histos_mc, nevts, x


def eval_efficiency(npass, nfail, sigma_npass, sigma_nfail):
    """
    """
    eff = npass/(npass+nfail)
    var1 = ((1-npass)**2)*(sigma_npass**2)
    var2 = (npass**2)*(sigma_nfail**2)
    sigma_eff = ROOT.TMath.Sqrt(var1+var2)/((npass+nfail)**2)

    return eff, sigma_eff


def add_result(dict_results, res_pass, res_fail, eff, bin_pt, bin_eta):
    """
    """
    dict_results.update({
       f"{bin_pt},{bin_eta}": {
            "efficiency": eff,
            "fit_stat_pass": {
                "migrad_status": res_pass.status(),
                "parameters": res_pass.floatParsFinal(),
                "cov_matrix": res_pass.covarianceMatrix(),
                "cov_matrix_quality": res_pass.covQual(),
                "global_correlation": res_pass.globalCorr(),
                "EDM": res_pass.edm()
                },
            "fit_stat_fail": {
                "migrad_status": res_fail.status(),
                "parameters": res_fail.floatParsFinal(),
                "cov_matrix": res_fail.covarianceMatrix(),
                "cov_matrix_quality": res_fail.covQual(),
                "global_correlation": res_fail.globalCorr(),
                "EDM": res_fail.edm()
                }
            }
        })
    return dict_results


if __name__ == '__main__':
    '''
    file = ROOT.TFile("root_files/tnp_iso_mc.root")
    histo = file.pass_mu_DY_postVFP
    th3_checks(histo)
    '''
