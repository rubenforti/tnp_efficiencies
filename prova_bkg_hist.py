
import ROOT
from array import array
from workspace_config import get_roohist


type_eff = 'iso'


f_data = ROOT.TFile(f"root_files/tnp_{type_eff}_data.root")
f_bkg = ROOT.TFile(
    "/scratchnvme/rajarshi/Bkg_TNP_3D_Histograms/OS/tnp_reco_TTSemileptonic_vertexWeights1_oscharge1.root")


bins_mass = array('d', [50 + i for i in range(81)])

x = ROOT.RooRealVar(
    "x", "TP M_inv", bins_mass[0], bins_mass[-1], unit="GeV/c^2")

histos_pass = get_roohist((f_data, f_bkg), x, 6, 25, 'fail')


c = ROOT.TCanvas()
c.cd()
ROOT.gPad.SetLeftMargin(0.15)
frame = x.frame()
histos_pass[1].plotOn(frame)

frame.Draw()
c.SaveAs("provaBkg_fail.png")
