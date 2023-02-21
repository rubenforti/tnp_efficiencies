"""
"""

import ROOT


def makeGaussianHisto():
    hh = ROOT.TH1D("hh", "hh", 100, 50, 130)
    for i in range(100000):
        hh.Fill(ROOT.gRandom.Gaus(91, 2.5))
    return hh


def makeAndSavePlot(axis, histo, func, name='prova.png', title="Histo"):
    c = ROOT.TCanvas()
    c.cd()
    frame = axis.frame(Title=title+' '+str(axis))
    dh.plotOn(frame)
    gauss.plotOn(frame)
    frame.Draw()
    c.SaveAs(name)


if __name__ == '__main__':

    f = ROOT.TFile("root_files/tnp_iso_data.root")
    histo_pass = f.pass_mu_RunGtoH
    histo_fail = f.fail_mu_RunGtoH

    ymin = histo_pass.GetYaxis().GetXmin()
    ymax = histo_pass.GetYaxis().GetNbins()

    zmin = histo_pass.GetZaxis().GetXmin()
    zmax = histo_pass.GetZaxis().GetNbins()

    histo = histo_pass.ProjectionX("histo1", 0, ymax, 0, zmax)
    # histo = makeGaussianHisto()

    print(ymin, zmax, histo_pass.GetXaxis().GetNbins())

    x = ROOT.RooRealVar("x", "x", histo.GetXaxis().GetXmin(),
                        histo.GetXaxis().GetXmax())

    dh = ROOT.RooDataHist("dh", "dh", [x], Import=histo)

    mean = ROOT.RooRealVar("mean", "mean", 91, 85, 97)
    sigma = ROOT.RooRealVar("sigma", "sigma", 0.1, 5)
    gauss = ROOT.RooGaussian("gauss", "gauss", x, mean, sigma)

    gauss.fitTo(dh)

    makeAndSavePlot(x, dh, gauss)
