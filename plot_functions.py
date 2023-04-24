"""
"""
import ROOT
import sys
# import pickle

from array import array
from results_utilities import res_manager_indep


def makeAndSavePlot(axis, data, function, bkg_name='expo',
                    name='prova.png', title="Histo", pull=False):

    c = ROOT.TCanvas()
    if pull is True:
        c.Divide(2)

    c.cd(1)
    ROOT.gPad.SetLeftMargin(0.15)
    frame = axis.frame(Title=title+' '+str(axis))
    data.plotOn(frame)
    function.plotOn(frame, LineColor='kBlue')
    #[function.plotOn(frame, Components=comp, LineColor='kRed')
     #for comp in function.getComponents() if comp.GetName() == bkg_name]

    function.plotOn(frame, LineColor='kBlue')
    frame.Draw()

    if pull:
        c.cd(2)
        ROOT.gPad.SetLeftMargin(0.15)
        hpull = frame.pullHist()
        frame2 = axis.frame(Title="Residual Distribution")
        frame2.addPlotable(hpull, "P")
        frame2.Draw()
    c.SaveAs(name)


def differential_efficiency(filename, binning_pt=(), binning_eta=()):
    """
    """

    file = ROOT.TFile.Open("root_files/iso_indep_eff.root", "RECREATE")

    if len(binning_pt) == 0:
        binning_pt = array('d', [24., 26., 28., 30., 32., 34.,
                           36., 38., 40., 42., 44., 47., 50., 55., 60., 65.])
    if len(binning_eta) == 0:
        binning_eta = array('d', [round(-2.4 + i*0.1, 2) for i in range(49)])

    ROOT.gStyle.SetOptStat("n")

    results = res_manager_indep()
    results.open(filename)
    res_dict = results.dictionary()

    idx_pt, idx_eta = 0, 0

    histo = ROOT.TH2D("eff_histo", "Efficiency (pt,eta)", len(binning_pt)-1, binning_pt,
                      len(binning_eta)-1, binning_eta)
    histo_err = ROOT.TH2D(
        "eff_error_histo", "Relative error on efficiency (pt,eta)",
        len(binning_pt)-1, binning_pt, len(binning_eta)-1, binning_eta)

    for key in res_dict.keys():
        idx_pt = int(key.split(",")[0])
        idx_eta = int(key.split(",")[1])
        eff = res_dict[key]["efficiency"]
        # print(eff[0], eff[1])
        histo.SetBinContent(idx_pt, idx_eta, eff[0])
        histo.SetBinError(idx_pt, idx_eta, eff[1])
        histo_err.SetBinContent(idx_pt, idx_eta, eff[1]/eff[0])

    c1 = ROOT.TCanvas()
    c1.cd()
    histo.Draw("COLZ")
    c1.SaveAs("differential_eff_indep.png")

    c2 = ROOT.TCanvas()
    c2.cd()
    histo_err.Draw("COLZ")
    c2.SaveAs("differential_eff_indep_errors.png")

    file.Write("eff_histo")
    # file.Write("eff_error_histo")
    file.Close()


if __name__ == '__main__':
    print("EO")

    differential_efficiency("results/indep_eff_results.pkl")
