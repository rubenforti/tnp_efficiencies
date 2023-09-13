import ROOT

file_original = ROOT.TFile("/scratchnvme/wmass/Steve_root_files/Standard_SF_files/tnp_iso_mc_vertexWeights1_oscharge1.root", "READ")

file_new = ROOT.TFile("/scratchnvme/rajarshi/Signal_TNP_3D_Histograms/OS/tnp_iso_mc_vertexWeights1_oscharge1.root", "READ")


histo_p_old = file_original.Get("pass_mu_DY_postVFP")
histo_p_new = file_new.Get("pass_mu_DY_postVFP")

print(histo_p_old.GetNbinsX(), histo_p_old.GetNbinsY(), histo_p_old.GetNbinsZ())
print(histo_p_new.GetNbinsX(), histo_p_new.GetNbinsY(), histo_p_new.GetNbinsZ())

print("")

n_errors = 0

differences = 0

max_diff = 0

for j in range(1, 16):
    for k in range(1, 49):

        th1_old = histo_p_old.ProjectionX("th1_old", j, j, k, k)
        th1_new = histo_p_new.ProjectionX("th1_new", j, j, k, k)

        if th1_old.Integral() != th1_new.Integral():
            if max_diff < abs(th1_old.Integral()-th1_new.Integral()):
                max_diff = abs(th1_old.Integral()-th1_new.Integral())
            differences += (th1_old.Integral()-th1_new.Integral())
            print("Bin", j, k, "has different content")
            n_errors += 1

print(differences/n_errors)
print(max_diff)