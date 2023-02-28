"""
"""

import ROOT


def makeGaussianHisto():
    hh = ROOT.TH1D("hh", "hh", 100, 50, 130)
    for i in range(100000):
        hh.Fill(ROOT.gRandom.Gaus(91, 2.5))
    return hh


def makeAndSavePlot(axis, histo, *functions, name='prova.png', title="Histo"):
    c = ROOT.TCanvas()
    c.cd()
    func_list = []
    [func_list.append(func) for func in functions]
    sum = ROOT.RooAddPdf("pdf", "pdf", func_list, [0.5])
    frame = axis.frame(Title=title+' '+str(axis))
    histo.plotOn(frame)
    for func in functions:
        func.plotOn(frame)
    sum.plotOn(frame)
    frame.Draw()
    c.SaveAs(name)


def roodouble(x, xmin, xmax):
    return ROOT.RooRealVar(str(x), str(x), x, xmin, xmax)

def init_Gaussian(x, mean=[], sigma=[], name="gaussian", title="gaussian"):
    """
    """
    if not mean:
        mean = [91., 75., 105.]
    if not sigma:
        sigma = [2.5, 0.1, 5.]
    mu = ROOT.RooRealVar("mu", "mu", mean[0], mean[1], mean[2])
    s = ROOT.RooRealVar("sigma", "sigma", sigma[0], sigma[1], sigma[2])
    output = ROOT.RooGaussian(name, title, x, mu, s)
    return output


if __name__ == '__main__':

    f = ROOT.TFile("root_files/tnp_iso_data.root")
    histo_pass = f.pass_mu_RunGtoH
    histo_fail = f.fail_mu_RunGtoH

    ymin = histo_pass.GetYaxis().GetXmin()
    ymax = histo_pass.GetYaxis().GetNbins()

    zmin = histo_pass.GetZaxis().GetXmin()
    zmax = histo_pass.GetZaxis().GetNbins()

    x = roodouble(91, 50, 130)

    histo = makeGaussianHisto()

    dh = ROOT.RooDataHist("dh", "dh", [x], Import=histo)

    mean = ROOT.RooRealVar("mean", "mean", 91, 85, 97)
    sigma = ROOT.RooRealVar("sigma", "sigma", 0.1, 5)
    gauss = ROOT.RooGaussian("gauss", "gauss", x, mean, sigma)
    
    # gauss = init_Gaussian(x)    
    gauss.fitTo(dh)

    makeAndSavePlot(x, dh, gauss)
