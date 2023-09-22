import ROOT
from utilities.base_library import bin_dictionary
import matplotlib.pyplot as plt
import numpy as np


file = ROOT.TFile("bkg_studies/iso/ws_iso_backgrounds.root")

ws = file.Get("w")


bin_key = "[24.0to26.0][-2.4to-2.3]"

bin_dict = bin_dictionary("pt", "eta")

bkg_categories = ["WW", "WZ", "ZZ", "TTSemileptonic", "Ztautau", "SameCharge"]


axis = ws.var(f"x_pass_{bin_key}")

data = ws.data(f"Minv_data_pass_{bin_key}")

mc = ws.data(f"Minv_mc_pass_{bin_key}")
mc_pdf = ROOT.RooHistPdf("mc_pdf", "mc_pdf", ROOT.RooArgSet(axis), mc, 0)

tot_bkg = ROOT.RooDataHist("tot_bkg", "tot_bkg", ROOT.RooArgList(axis), "")

tot_bkg_pdf = ROOT.RooHistPdf("tot_bkg_pdf", "tot_bkg_pdf", ROOT.RooArgSet(axis), tot_bkg, 0)

for cat in bkg_categories:
    bkg_cat = ws.data(f"Minv_bkg_pass_{bin_key}_{cat}")
    tot_bkg.add(bkg_cat)

c = ROOT.TCanvas()
c.cd()

frame = axis.frame()
mc_pdf.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Normalization(mc.sumEntries(), ROOT.RooAbsReal.NumEvent))
data.plotOn(frame, ROOT.RooFit.MarkerColor(ROOT.kBlack))
tot_bkg_pdf.plotOn(frame, ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Normalization(tot_bkg.sumEntries(), ROOT.RooAbsReal.NumEvent))
frame.Draw()


c.SaveAs("bkg_vs_data.pdf")





'''


bkg_frac_fit, bkg_frac_real = [], []

sig_ratio = []

for bin_key in bin_dict.keys():

    for flag in ["pass", "fail"]:

        if flag == "fail": continue

        res  = ws.obj(f"results_{flag}_{bin_key}")

        if type(res) is ROOT.RooFitResult:
            pars = res.floatParsFinal()

            Nsig = pars.find(f"nsig_{flag}_{bin_key}")
            Nbkg = pars.find(f"nbkg_{flag}_{bin_key}")
    
            if type(Nbkg) is ROOT.RooRealVar: bkg_frac_fit.append(Nbkg.getVal()/(Nsig.getVal() + Nbkg.getVal()))

        data_hist = ws.data(f"Minv_data_{flag}_{bin_key}")
        Nsig_data = data_hist.sumEntries()

        mc_hist = ws.data(f"Minv_mc_{flag}_{bin_key}")
        Nsig_real = mc_hist.sumEntries()

        sig_ratio.append(Nsig_data/(Nsig_real*0.66*0.66))


        Nbkg_real = 0
        for cat in bkg_categories:
            bkg_cat_hist = ws.data(f"Minv_bkg_{flag}_{bin_key}_{cat}")
            Nbkg_real += bkg_cat_hist.sumEntries()

        bkg_frac_real.append(Nbkg_real/(Nsig_real + Nbkg_real))


bkg_frac_fit = np.array(bkg_frac_fit)
bkg_frac_real = np.array(bkg_frac_real)

print(bkg_frac_fit.mean(), bkg_frac_real.mean())
print(bkg_frac_fit.std(), bkg_frac_real.std())

plt.plot()
plt.subplot(1, 3, 1)
plt.hist(bkg_frac_fit, bins=20, color="red")
plt.subplot(1, 3, 2)
plt.hist(bkg_frac_real, bins=20, color="blue")
plt.subplot(1, 3, 3)
plt.hist(sig_ratio, bins=20, color="green")
plt.show()
plt.savefig("bkg_frac_pass.png")
'''