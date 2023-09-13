import ROOT

from utilities.base_library import binning, bin_dictionary

"""
Modulo di test che controlla se il numero totale di eventi di bkg fittato in un macro-bin (sui dati)
è uguale alla somma dlla stima di nkbg nei micro-bin corrispondenti
"""


file_old = ROOT.TFile("results/benchmark_iso/ws_iso_indep_benchmark.root", "READ")
ws_old = file_old.Get("w")

file_new = ROOT.TFile("root_files/ws_bkg_studies.root", "READ")
ws_new = file_new.Get("w")

bin_dict_old = bin_dictionary()
bin_dict_new = bin_dictionary("pt", "eta_8bins")

idx_pt = 1
idx_eta = 1

counter_1sigma = 0
counter_3sigma = 0
counter_over3sigma = 0

for bin_key in bin_dict_new:

    gl_bins, bins_pt, bins_eta = bin_dict_new[bin_key]
    
    res_n = ws_new.obj(f"results_pass_{bin_key}")
    pars_n = res_n.floatParsFinal()
    nbkg_n_obj = pars_n.find(f"nbkg_pass_{bin_key}")
    if type(nbkg_n_obj) is ROOT.RooRealVar:
        nbkg_n = [nbkg_n_obj.getVal(), nbkg_n_obj.getError()]
    else:
        nbkg_n = [0, 0]
        
    nbkg_o = [0,0]

    for i in range(bins_pt[0], bins_pt[-1]+1):
        for j in range(bins_eta[0], bins_eta[-1]+1):
            print(i,j)
            res_o = ws_old.obj(f"results_pass_({i}|{j})")
            pars_o = res_o.floatParsFinal()
            nbkg_o_obj = pars_o.find(f"nbkg_pass_({i}|{j})")
            if type(nbkg_o_obj) is ROOT.RooRealVar:
                nbkg_o[0] += nbkg_o_obj.getVal()
                nbkg_o[1] += nbkg_o_obj.getError()**2
            else:
                pass
    
    nbkg_o[1] = nbkg_o[1]**0.5

    print(f"New: {round(nbkg_n[0])}+-{round(nbkg_n[1])}", f"Old: {round(nbkg_o[0])}+-{round(nbkg_o[1])}")
    
    sigma = (nbkg_n[1]**2 + nbkg_o[1]**2)**0.5 


    nsigma = (nbkg_n[0] - nbkg_o[0])/sigma if sigma != 0 else 999

    if nsigma <= 1:
        counter_1sigma += 1
    
    if nsigma <= 3:
        counter_3sigma += 1
    elif nsigma > 3:
        counter_over3sigma += 1

    
print("*********")
print(counter_1sigma, counter_3sigma, counter_over3sigma)
print("*********")






