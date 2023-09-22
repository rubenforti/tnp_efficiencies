"""
"""

import ROOT
from utilities.plot_utils import style_settings
from utilities.CMS_lumi import CMS_lumi

ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

bkg_cat = "TTSemileptonic"


file = ROOT.TFile("bkg_studies/iso/bkg_2d_distrib.root", "READ")

hist_pass = file.Get(f"{bkg_cat}_pass_2d")
hist_fail = file.Get(f"{bkg_cat}_fail_2d")

c = ROOT.TCanvas("c", "c", 1600, 900)
c.cd()

style_settings()

pad_title = ROOT.TPad("pad_title", "pad_title", 0, 0.94, 1, 1)
pad_subtitle_pass = ROOT.TPad("pad_subtitle_pass", "pad_subtitle_pass", 0, 0.9, 0.5, 0.94)
pad_subtitle_fail = ROOT.TPad("pad_subtitle_fail", "pad_subtitle_fail", 0.5, 0.9, 1, 0.94)
pad_pass = ROOT.TPad("pad_pass", "pad_plot", 0, 0, 0.5, 0.9)
pad_fail = ROOT.TPad("pad_fail", "pad_plot", 0.5, 0, 1, 0.9)

pad_title.SetMargin(0.1, 0.1, 0.1, 0.1), pad_title.Draw()
pad_subtitle_pass.SetMargin(0.1, 0.1, 0.1, 0.1), pad_subtitle_pass.Draw()
pad_subtitle_fail.SetMargin(0.1, 0.1, 0.1, 0.1), pad_subtitle_fail.Draw()
pad_pass.SetMargin(0.1, 0.135, 0.12, 0.05), pad_pass.Draw()
pad_fail.SetMargin(0.1, 0.135, 0.12, 0.05), pad_fail.Draw()


pad_title.cd()
titlebox = ROOT.TPaveText(0, 0, 1, 1, "NDC NB")
titlebox.SetFillColor(0)
# titlebox.SetTextFont(42)
titlebox.SetTextSize(0.6)
titlebox.AddText(f"{bkg_cat}  2D distribution")
titlebox.Draw()
c.Update()

pad_subtitle_pass.cd()
title_pass = ROOT.TPaveText(0, 0, 1, 1, "NDC NB")
title_pass.SetFillColor(0)
title_pass.SetTextSize(0.6)
title_pass.AddText(f"Passing probes")
title_pass.Draw()
c.Update()

pad_subtitle_fail.cd()
title_fail = ROOT.TPaveText(0, 0, 1, 1, "NDC NB")
title_fail.SetFillColor(0)
title_fail.SetTextSize(0.6)
title_fail.AddText(f"Failing probes")
title_fail.Draw()
c.Update()

pad_pass.cd()
hist_pass.Draw("colz")
hist_pass.SetTitle("")
hist_pass.SetTitleSize(0)
Yaxis_pass = hist_pass.GetYaxis()
Yaxis_pass.SetTitle("#eta^{#mu}")
Yaxis_pass.SetTitleSize(0.035)
Yaxis_pass.SetTitleOffset(1.2)
Xaxis_pass = hist_pass.GetXaxis()
Xaxis_pass.SetTitle("p_{T}^{#mu} [GeV]")
Xaxis_pass.SetTitleSize(0.035)
Xaxis_pass.SetTitleOffset(1.2)
CMS_lumi(pad_pass, 5, 0, simulation=True)
pad_pass.Update()


pad_fail.cd()
hist_fail.Draw("colz")
hist_fail.SetTitle("")
hist_fail.SetTitleSize(0)
Yaxis_fail = hist_fail.GetYaxis()
Yaxis_fail.SetTitle("#eta^{#mu}")
Yaxis_fail.SetTitleSize(0.035)
Yaxis_fail.SetTitleOffset(1.2)
Xaxis_fail = hist_fail.GetXaxis()
Xaxis_fail.SetTitle("p_{T}^{#mu} [GeV]")
Xaxis_fail.SetTitleSize(0.035)
Xaxis_fail.SetTitleOffset(1.2)
CMS_lumi(pad_fail, 5, 0, simulation=True)
pad_fail.Update()


c.SaveAs(f"bkg_studies/iso/{bkg_cat}_distrib_2d.pdf")

