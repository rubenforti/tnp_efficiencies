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

    histo_th1 = histo_th3.ProjectionX(f"histo_mass_{flag}", bin_pt[0], bin_pt[1], bin_eta[0], bin_eta[1])
    
    xAxis = histo_th1.GetXaxis()
    x = ROOT.RooRealVar("x", "x", xAxis.GetXmin(), xAxis.GetXmax(), unit="GeV/c^2")
    print(f"Num TH1 entries = {histo_th1.GetEntries()}")
    roohist = ROOT.RooDataHist(f"roohist_{flag}", f"roohist_{flag}", [x], Import=histo_th1)
    print(f"Num RooDataHist entries = {roohist.numEntries()}") 
    return roohist, histo_th1.GetEntries()


def import_Steve_histos(type_eff, bin_pt, bin_eta):

    if len(bin_pt) == 1:
        bin_pt.append(bin_pt[0])
        bin_pt[0] -= 1
    if len(bin_eta) == 1:
        bin_eta.append(bin_eta[0])
        bin_eta[0] -= 1
    
    x = ROOT.RooRealVar("x", "TP M_{inv}", 50, 130, unit="GeV/c^2") 

    # Import of the 3D histograms
    f_data = ROOT.TFile(f"root_files/tnp_{type_eff}_data.root")
    f_mc = ROOT.TFile(f"root_files/tnp_{type_eff}_mc.root")

    h_pass_data, nev_pass_data = profile_histo(f_data.pass_mu_RunGtoH, x, bin_pt, bin_eta, 1)
    h_fail_data, nev_fail_data = profile_histo(f_data.fail_mu_RunGtoH, x, bin_pt, bin_eta, 2)
    h_pass_mc, nev_pass_mc = profile_histo(f_mc.pass_mu_DY_postVFP, x, bin_pt, bin_eta, 3)
    h_fail_mc, nev_fail_mc = profile_histo(f_mc.fail_mu_DY_postVFP, x, bin_pt, bin_eta, 4)
    
    nevts = ((nev_fail_data, nev_pass_data), (nev_fail_mc, nev_pass_mc))
    histos_data = (h_fail_data, h_pass_data)
    histos_mc = (h_fail_mc, h_pass_mc)

    return histos_data, histos_mc, nevts, x


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
    

'''
def init_Gaussian(x, mean=91., sigma=1.2, name="gaussian", title="gaussian"):
    """
    """
    mu = ROOT.RooRealVar("mu", "mu", mean, 80, 100)
    s = ROOT.RooRealVar("sigma", "sigma", sigma, 0.1, 5)
    return ROOT.RooGaussian(name, title, x, mu, s)
'''

'''
def init_CrystalBall(x, mean=91., sigma=1.2, n=1, alpha=1.5, name="CB", title="CB"):
    """
    """
    mu = roodouble(mean, 80, 100)
    s = roodouble(sigma, 0.1, 5)
    nn = roodouble(n, 0, 4)
    al = roodouble(alpha, -5, 5)
    return ROOT.RooCBShape(name, title, x, mu, s, nn, al)
'''

if __name__ == '__main__':

    '''
    f = ROOT.TFile("root_files/tnp_iso_data.root")
    histo_pass = f.pass_mu_RunGtoH
    histo_fail = f.fail_mu_RunGtoH

    ymin = histo_pass.GetYaxis().GetXmin()
    ymax = histo_pass.GetYaxis().GetNbins()

    zmin = histo_pass.GetZaxis().GetXmin()
    zmax = histo_pass.GetZaxis().GetNbins()

    c = ROOT.TCanvas()
    c.cd()

    x = ROOT.RooRealVar("x", "x", 91, 50, 130)
    frame = x.frame(Title='prova')
    histo = makeGaussianHisto()

    dh = ROOT.RooDataHist("dh", "dh", [x], Import=histo)

    mean = ROOT.RooRealVar("mean", "mean", 91, 85, 97)
    sigma = ROOT.RooRealVar("sigma", "sigma", 0.1, 5)
    gauss_1 = ROOT.RooGaussian("gauss", "gauss", x, mean, sigma)
    '''

    i,b,c,d=import_Steve_histos("iso", [1], [1])
    print(c)
