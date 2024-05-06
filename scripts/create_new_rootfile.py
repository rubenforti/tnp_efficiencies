import ROOT
import os
from copy import deepcopy as dcp

file_in = ROOT.TFile("allSmooth_GtoHout_vtxAgnIso.root", "READ")

hists = []

idx=0


'''
for hist_name_in in file_in.GetListOfKeys():

    #print(hist_name_in.GetName())

    #print(file_in.Get(hist_name_in.GetName()).GetName())
    
    if "plus" in file_in.Get(hist_name_in.GetName()).GetName(): 
        ch = "plus"
    elif "minus" in file_in.Get(hist_name_in.GetName()).GetName(): 
        ch = "minus"
    else: 
        ch = "all"
    
    if (file_in.Get(hist_name_in.GetName()).GetName() != f"SF_nomiAndAlt_GtoH_tracking_{ch}") or (ch=="all"):
        file_in.cd()
        h_in = file_in.Get(hist_name_in.GetName()).Clone()
        hists.append(h_in)
        # dprint("A", type(h_in), idx)
    else:
        file_in_new = ROOT.TFile(f"smoothLeptonScaleFactors/GtoH/mu_tracking_{ch}/smoothedSFandEffi_tracking_GtoH_{ch}.root", "READ")
        file_in_new.cd()
        h_in_new = file_in_new.Get(f"SF_nomiAndAlt_GtoH_tracking_{ch}")
        h_in_new.SetName(f"SF_nomiAndAlt_GtoH_tracking_{ch}")
        hists.append(h_in_new)
        # print("B", file_in.Get(hist_name_in.GetName()).GetName(), ch, idx)  
    
    idx+=1
'''

file_out = ROOT.TFile("allSmooth_GtoHout_vtxAgnIso_altBkg.root", "RECREATE")

for ch in ["plus", "minus"]:
    file_in_new = ROOT.TFile(f"smoothLeptonScaleFactors/GtoH/mu_tracking_{ch}/smoothedSFandEffi_tracking_GtoH_{ch}.root", "READ")
    #file_in_new.cd()
    h_in_new = file_in_new.Get(f"SF_nomiAndAlt_GtoH_tracking_{ch}").Clone()
    h_in_new.SetName(f"SF_nomiAndAlt_GtoH_tracking_{ch}_altBkg")
    hists.append(h_in_new)

    file_out.cd()
    h_in_new.Write()

idx_ctrl=0
for hist in hists:
    print(type(hist), idx_ctrl)
    idx_ctrl+=1
    #hist.Write()
    #print(hist.GetName())


file_out.Close()

