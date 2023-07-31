"""
"""

import sys
import os
import ROOT
from array import array


binnings = {
    "pt": array('d', [24., 26., 28., 30., 32., 34., 36., 38., 40., 42., 44., 47., 50., 55., 60., 65.]),
    "eta" : array('d', [round(-2.4 + i*0.1, 2) for i in range(49)]),
    "mass_50_130" : array('d', [50 + i for i in range(81)]),
    "mass_60_120" : array('d', [60 + i for i in range(61)]),
    "pt_12bins" : array('d', [24., 26., 28., 32., 36., 40., 42., 44., 47., 50., 55., 60., 65.]),
    "pt_10bins" : array('d', [24., 28., 32., 36., 40., 44., 47., 50., 55., 60., 65.]),
    "pt_8bins" : array('d', [24., 30., 36., 40., 44., 50., 55., 60., 65.]),
    "pt_6bins" : array('d', [24., 30., 36., 42., 47., 55., 65.]),
    "eta_24bins" : array('d', [round(-2.4 + i*0.2, 2) for i in range(25)]),
    "eta_16bins" : array('d', [round(-2.4 + i*0.3, 2) for i in range(17)]),
    "eta_8bins" : array('d', [round(-2.4 + i*0.6, 2) for i in range(9)]),
    "eta_4bins" : array('d', [round(-2.4 + i*1.2, 2) for i in range(5)]),
    "mass_2GeV" : array('d', [60 + 2*i for i in range(31)]),
    "mass_3GeV" : array('d', [60 + 3*i for i in range(21)]),
    "mass_4GeV" : array('d', [60 + 4*i for i in range(16)]),
    }

lumi_data = 16.8  # fb^-1

# sig_mc_repo = "/scratchnvme/rajarshi/Signal_TNP_3D_Histograms/OS"
# bkg_repo = "/scratchnvme/rajarshi/Bkg_TNP_3D_Histograms/OS"

sig_mc_repo = "root_files/datasets"
bkg_repo = "root_files/datasets"

BR_TAUToMU = 0.1739
BR_TAUToE = 0.1782
Z_TAU_TO_LEP_RATIO = (1.-(1. - BR_TAUToMU - BR_TAUToE)**2)
xsec_ZmmPostVFP = 2001.9

xsec_bkg = {
    # Unit = pb
    "WW" : 12.6,
    "WZ" : 5.4341,
    "ZZ" : 0.60,
    "TTSemileptonic" : 366.34, 
    "Ztautau" : xsec_ZmmPostVFP*Z_TAU_TO_LEP_RATIO
}

weightSum_signal = 50602137.0


###############################################################################

def binning(type):
    """
    Returns the binning array for the given type of variable. The type can be "pt", "eta", 
    "mass_50_130" or "mass_60_120"
    """
    return binnings[type]

###############################################################################

def get_idx_from_bounds(bounds_pt, bounds_eta):
    """
    Function that, given the bin bounds (that have to be present in the binning
    arrays!!!), returns the coordinates of the selected region. The coordinates 
    are given as list of global indexes and as bounds on pt/eta indexes (max 
    and min, since the region is rectangular)
    """
    initial_dict = bin_dictionary("pt", "eta")
    global_bins = []

    bounds_pt_idx = []
    bounds_eta_idx = []

    for el in initial_dict:

        gl_idx, idx_pt, idx_eta = initial_dict[el]

        string_pt, string_eta = el.split("][")
        pt_min, pt_max = string_pt[1:].split("to")
        eta_min, eta_max = string_eta[:-1].split("to")

        pt_min, pt_max = float(pt_min), float(pt_max)
        eta_min, eta_max = float(eta_min), float(eta_max)

        if(pt_min>=bounds_pt[0]) and (pt_max<=bounds_pt[1]) and \
          (eta_min>=bounds_eta[0]) and (eta_max<=bounds_eta[1]):
            global_bins.append(gl_idx)
            bounds_pt_idx.append(idx_pt)
            bounds_eta_idx.append(idx_eta)
            
    return global_bins, [min(bounds_pt_idx), max(bounds_pt_idx)], [min(bounds_eta_idx), max(bounds_eta_idx)]

###############################################################################

def bin_dictionary(binning_pt_name="pt", binning_eta_name="eta"):
    """
    Creates a dictionary that relates the bounds of every bin (keys) to the 
    indexes useful for the bin selection. The indexes are returned in a 
    3-element list containing the global index, the pt index and the eta index
    in this order.
    """
    index_dictionary = {}
    global_idx = 1

    binning_pt, binning_eta = binnings[binning_pt_name], binnings[binning_eta_name]

    for idx_pt in range(1, len(binning_pt)):
        for idx_eta in range(1, len(binning_eta)):

            bin_key = f"[{binning_pt[idx_pt-1]}to{binning_pt[idx_pt]}][{binning_eta[idx_eta-1]}to{binning_eta[idx_eta]}]"

            if binning_pt_name == "pt" and binning_eta_name == "eta":
                index_dictionary.update({bin_key : [global_idx, idx_pt, idx_eta]})
                global_idx +=1
            else:
                global_idx, bounds_idx_pt, bounds_idx_eta = get_idx_from_bounds(
                    [binning_pt[idx_pt-1], binning_pt[idx_pt]], [binning_eta[idx_eta-1], binning_eta[idx_eta]])
                index_dictionary.update({bin_key : [global_idx, bounds_idx_pt, bounds_idx_eta]})

    return index_dictionary

###############################################################################

def cross_section_bkg(bkg_process):
    """
    Returns the cross section of the given background process in pb
    """
    return xsec_bkg[bkg_process]

###############################################################################

def lumi_factors(type_eff, bkg_categories):
    """
    "Lumi scale" defined as alpha that satisifies lumi_data=alpha*lumi_bkg
    """
    lumi_scales = {}

    # file_sig = ROOT.TFile(f"{sig_mc_repo}/tnp_{type_eff}_mc_vertexWeights1_oscharge1.root")
    # wsum_histo_sig = file_sig.Get("weightSum")
    # num_init_sig = wsum_histo_sig.Integral()
    num_init_sig = weightSum_signal
    xsection_sig = xsec_ZmmPostVFP*1000  # has to be put in fb
    lumi_sig = num_init_sig/xsection_sig
    lumi_scales.update({"Zmumu" : lumi_data/lumi_sig})

    for cat in bkg_categories:
        file = ROOT.TFile(f"{bkg_repo}/tnp_{type_eff}_{cat}_vertexWeights1_oscharge1.root")
        
        wsum_histo = file.Get("weightSum")
        num_init = wsum_histo.Integral()
        xsection = cross_section_bkg(cat)*1000  # has to be put in fb
        lumi_bkg = num_init/xsection

        scale = lumi_data/lumi_bkg

        lumi_scales.update({cat : scale})
    
    return lumi_scales

###############################################################################

def eval_efficiency(npass, nfail, sigma_npass, sigma_nfail):
    """
    """
    eff = npass/(npass+nfail)
    var1 = (nfail**2)*(sigma_npass**2)
    var2 = (npass**2)*(sigma_nfail**2)
    sigma_eff = ROOT.TMath.Sqrt(var1+var2)/((npass+nfail)**2)

    return eff, sigma_eff

###############################################################################

def sumw2_error(histo):
    """
    Calculates the error on the integral of a RooDataHist as the (square root 
    of the) sum of the errors associated to each bin. The errors considered are
    the "SumW2", already stored in the RooDataHist.
    """

    variance = 0
    for i in range(0, histo.numEntries()):
        histo.get(i)
        variance += histo.weightError(ROOT.RooAbsData.SumW2)**2
    sum_error = variance**0.5

    return sum_error


###############################################################################
###############################################################################


if __name__ == '__main__':
    
    x = ROOT.RooRealVar("x", "x", 0, 100)
    mu = ROOT.RooRealVar("mu", "mu", 50, 0, 100)
    sigma = ROOT.RooRealVar("sigma", "sigma", 10, 0.5, 100)

    unif = ROOT.RooUniform("gaussian", "gaussian", x)

    data = unif.generateBinned(ROOT.RooArgSet(x), 100000)

    histo = data.createHistogram("histo", x, ROOT.RooFit.Binning(100, 0, 100))
    histo2 = data.createHistogram("histo", x, ROOT.RooFit.Binning(100, 0, 100))

 

    histo.Scale(1.5)
    histo2.Scale(2.7)


    histo.Add(histo2)

    roohisto = ROOT.RooDataHist("roohisto", "roohisto", ROOT.RooArgList(x), histo)

    print(roohisto.sumEntries())
    print(sumw2_error(roohisto))
