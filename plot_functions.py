"""
"""

import sys
import pickle
import ROOT


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


def differential_eff_plot(file, binning_pt=(), binning_eta=()):
    """
    """
    if len(binning_pt) == 0:
        binning_pt = tuple([24., 26., 28., 30., 32., 34.,
                           36., 38., 40., 42., 44., 47., 50., 55., 60., 65.])
    if len(binning_eta) == 0:
        binning_eta = tuple([round(-2.4 + i*0.1, 2) for i in range(49)])

    histo = ROOT.TH2D("eff_histo", "Efficiency (pt,eta)", len(binning_pt)-1, binning_pt[0], binning_pt[-1],
                      len(binning_eta)-1, binning_eta[0], binning_eta[-1])

    with open("indep_eff_results.pkl", "rb") as f:
        results = pickle.load(f)

    idx_pt, idx_eta = 0, 0

    for key in results.keys():
        idx_pt = int(key.split(",")[0])
        idx_eta = int(key.split(",")[1])
        eff_val = results[key]["efficiency"][0]
        histo.Fill((binning_pt[idx_pt-1]+binning_pt[idx_pt])/2,
                   (binning_eta[idx_eta-1]+binning_eta[idx_eta])/2, eff_val)

    c = ROOT.TCanvas()
    c.cd()
    histo.Draw("COLZ")
    c.SaveAs("differential_eff.png")


if __name__ == '__main__':
    with open("indep_eff_results.pkl", "rb") as f:
        results = pickle.load(f)
    differential_eff_plot(results)
