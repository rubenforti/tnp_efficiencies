"""
"""

import sys
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

    # Option "e" has to be activated? Not clear how "errors are computed"
    histo_th1 = histo_th3.ProjectionX(
        f"histo_mass_{flag}", bin_pt[0], bin_pt[1], bin_eta[0], bin_eta[1])

    xAxis = histo_th1.GetXaxis()
    x = ROOT.RooRealVar("x", "x", xAxis.GetXmin(),
                        xAxis.GetXmax(), unit="GeV/c^2")
    print(f"Num TH1 entries = {histo_th1.GetEntries()}")
    roohist = ROOT.RooDataHist(
        f"roohist_{flag}", f"roohist_{flag}", [x], Import=histo_th1)
    print(f"Num RooDataHist entries = {roohist.numEntries()}")
    return roohist, histo_th1.GetEffectiveEntries()


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
    
    if len(bin_pt) == 1:
        bin_pt.append(bin_pt[0])

    if len(bin_eta) == 1:
        bin_eta.append(bin_eta[0])

    x = ROOT.RooRealVar("x", "TP M_{inv}", 50, 130, unit="GeV/c^2")

    # Import of the 3D histograms
    f_data = ROOT.TFile(f"root_files/tnp_{type_eff}_data.root")
    f_mc = ROOT.TFile(f"root_files/tnp_{type_eff}_mc.root")

    h_pass_data, nev_pass_data = profile_histo(
        f_data.pass_mu_RunGtoH, x, bin_pt, bin_eta, 1)
    h_pass_data.SetNameTitle(f"Events {type_eff} pass", f"Events {type_eff} pass")
    h_fail_data, nev_fail_data = profile_histo(
        f_data.fail_mu_RunGtoH, x, bin_pt, bin_eta, 2)
    h_fail_data.SetNameTitle(f"Events {type_eff} fail", f"Events {type_eff} fail")
    h_pass_mc, nev_pass_mc = profile_histo(
        f_mc.pass_mu_DY_postVFP, x, bin_pt, bin_eta, 3)
    h_fail_mc, nev_fail_mc = profile_histo(
        f_mc.fail_mu_DY_postVFP, x, bin_pt, bin_eta, 4)

    nevts = ((nev_fail_data, nev_pass_data), (nev_fail_mc, nev_pass_mc))
    histos_data = (h_fail_data, h_pass_data)
    histos_mc = (h_fail_mc, h_pass_mc)

    return histos_data, histos_mc, nevts, x


def llr_test_bkg(histo, pdf):
    """
    """
    
    pars_set = pdf.getParameters(histo)
    
    profiled_var = ROOT.RooRealVar()

    for idx in range(pars_set.getSize()):
        if pars_set[idx].GetName() == 'nbkg':
            profiled_var = pars_set[idx]

    null_profiled_var = ROOT.RooRealVar("nbkg", "nbkg", 0, 0, 0)
    
    llr_obj = ROOT.RooStats.ProfileLikelihoodCalculator(
            histo, pdf, ROOT.RooArgSet(profiled_var), 0.05, ROOT.RooArgSet(null_profiled_var))
    

    test_res = llr_obj.GetHypoTest()

    pval = test_res.NullPValue()

    print(f"p-value for null hypo is: {pval}")

    null = True if pval>0.05 else False

    return null




def makeGaussianHisto():
    hh = ROOT.TH1D("hh", "hh", 100, 50, 130)
    for i in range(100000):
        hh.Fill(ROOT.gRandom.Gaus(91, 2.5))
    return hh


def makeAndSavePlot(axis, data, function, name='prova.png', title="Histo", pull=False):

    c = ROOT.TCanvas()
    if pull is True:
        c.Divide(2)

    c.cd(1)
    ROOT.gPad.SetLeftMargin(0.15)
    frame = axis.frame(Title=title+' '+str(axis))
    data.plotOn(frame)
    for comp in function.getComponents():
        print(comp.GetName())
        function.plotOn(frame, Components=comp, LineStyle=':')
    function.plotOn(frame)
    frame.Draw()

    if pull:
        c.cd(2)
        ROOT.gPad.SetLeftMargin(0.15)
        hpull = frame.pullHist()
        frame2 = axis.frame(Title="Residual Distribution")
        frame2.addPlotable(hpull, "P")
        frame2.Draw()
    c.SaveAs(name)


def eval_efficiency(npass, nfail, sigma_npass, sigma_nfail):
    """
    """
    eff = npass/(npass+nfail)
    var1 = ((1-npass)**2)*(sigma_npass**2)
    var2 = (npass**2)*(sigma_nfail**2)
    sigma_eff = ROOT.TMath.Sqrt(var1+var2)/((npass+nfail)**2)

    return eff, sigma_eff


def pearson_chi2_eval(histo, pdf, nbins, res):
    """
    """
    npars = nbins - res.floatParsFinal().getSize()
    print(npars)
    chi2_sqrtvar = (2*npars)**(1/2.)
    print(f"Expected chi2 pars: mu={npars}, sqrt(var)={chi2_sqrtvar}")

    chi2_obj = ROOT.RooChi2Var("chi2", "chi2", pdf, histo, Verbose=False)
    print(f"Measured chi2 = {chi2_obj.getVal()}")
    print(f"Distance in sigma = {(chi2_obj.getVal()-npars)/chi2_sqrtvar}")


def add_result(dict_results, res_pass, res_fail, eff, bin_pt, bin_eta):
    """
    """
    
    dict_results.update({
       f"[{bin_pt}, {bin_eta}]" : {
            "efficiency" : eff,
            "fit_stat_pass" : {
                "migrad_status" : res_pass.status(), 
                "parameters" : res_pass.floatParsFinal(),
                "cov_matrix" : res_pass.covarianceMatrix(),
                "cov_matrix_quality" : res_pass.covQual(),
                "global_correlation" : res_pass.globalCorr(),
                "EDM" : res_pass.edm()
                },
            "fit_stat_fail" : {
                "migrad_status" : res_fail.status(),
                "parameters" : res_fail.floatParsFinal(),
                "cov_matrix" : res_fail.covarianceMatrix(),
                "cov_matrix_quality" : res_fail.covQual(),
                "global_correlation" : res_fail.globalCorr(),
                "EDM" : res_fail.edm()
                }
            }
        })

    return dict_results


if __name__ == '__main__':
    file = ROOT.TFile("root_files/tnp_iso_mc.root")
    histo = file.pass_mu_DY_postVFP
    th3_checks(histo)
