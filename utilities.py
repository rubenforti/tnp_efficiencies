"""
"""

import ROOT


def import_pdf_library(*functions):
    """
    """
    for function in functions:
        header_incl = ' #include "libCpp/'+function+'.h"'
        sourcefile = 'libCpp/'+function+'.cc'

        ctrl_head = ROOT.gInterpreter.Declare(header_incl)
        ctrl_source = ROOT.gSystem.CompileMacro(sourcefile, opt="ks")

        if not ctrl_head:
            print("ERROR in header loading")
            quit()
        if not ctrl_source == 1:
            print("ERROR in sourcefile compiling and loading")
            quit()


def import_histo(histo_th3, bin_pt, bin_eta):
    """
    Returns a RooDataHist of the variable TP_mass, in a (pt, eta) bin, selected
    from the TH3 given as input.
    """
    if len(bin_pt) == 1:
        bin_pt.append(bin_pt[0])
        bin_pt[0] -= 1
    if len(bin_eta) == 1:
        bin_eta.append(bin_eta[0])
        bin_eta[0] -= 1

    histo_th1 = histo_th3.ProjectionX(
        "histo_mass", bin_pt[0], bin_pt[1], bin_eta[0], bin_eta[1])
    xAxis = histo_th1.GetXaxis()

    x = ROOT.RooRealVar("x", "x", xAxis.GetXmin(),
                        xAxis.GetXmax(), unit="GeV/c^2")

    roohist = ROOT.RooDataHist("roohist", "roohist", [x], Import=histo_th1)

    n_events = histo_th1.GetEntries()
    print(n_events)

    return roohist, x, n_events


def makeGaussianHisto():
    hh = ROOT.TH1D("hh", "hh", 100, 50, 130)
    for i in range(100000):
        hh.Fill(ROOT.gRandom.Gaus(91, 2.5))
    return hh


def makeAndSavePlot(axis, histo, function, name='prova.png', title="Histo"):
    c = ROOT.TCanvas()
    c.cd()
    frame = axis.frame(Title=title+' '+str(axis))
    for comp in function.getComponents():
        print(comp.GetName())
        function.plotOn(frame, Components={comp}, LineStyle=':')
    histo.plotOn(frame)
    function.plotOn(frame)
    frame.Draw()
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

    makeAndSavePlot(x, dh, gauss_1)
