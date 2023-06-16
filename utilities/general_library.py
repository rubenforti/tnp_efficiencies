"""
"""

import sys
import os
import ROOT
from array import array


binning_pt = array('d', [24., 26., 28., 30., 32., 34., 36., 38., 40., 42., 44., 47., 50., 55., 60., 65.])
binning_eta = array('d', [round(-2.4 + i*0.1, 2) for i in range(49)])
binning_mass = array('d', [60 + i for i in range(61)])


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


def binning_dict():
    """
    Creates a dictionary that relates the bounds of every bin (keys) to the indexes useful for the bin 
    selection. The indexes are returned in a 3-element list containing the global index, the pt index
    and the eta index in this order.
    """

    index_dictionary = {}
    global_idx = 1
    for idx_pt in range(1, len(binning_pt)):
        for idx_eta in range(1, len(binning_eta)):

            indexes = [global_idx, idx_pt, idx_eta]
            bin_elem = {f"[{binning_pt[idx_pt-1]},{binning_pt[idx_pt]}][{binning_eta[idx_eta-1]},{binning_eta[idx_eta]}]" : indexes }
            index_dictionary.update(bin_elem)
            global_idx +=1
    
    return index_dictionary


def get_indexes_from_bounds(bounds_pt, bounds_eta):
    """
    Function that, given the bin bounds (that have to be present in the binning arrays!!!), returns the 
    coordinates of the selected region. The coordinates are given as list of global indexes and as bounds
    on pt/eta indexes (max and min, since the region is rectangular)
    """

    dict_idx = get_binning_dictionary()
    global_bins = []

    bounds_pt_idx = []
    bounds_eta_idx = []

    for el in dict_idx:

        gl_idx, idx_pt, idx_eta = dict_idx[el]



        string_pt, string_eta = el.split("][")
        pt_min, pt_max = string_pt[1:].split(",")
        eta_min, eta_max = string_eta[:-1].split(",")

        pt_min, pt_max = float(pt_min), float(pt_max)
        eta_min, eta_max = float(eta_min), float(eta_max)

        if(pt_min>=bounds_pt[0]) and (pt_max<=bounds_pt[1]) and \
          (eta_min>=bounds_eta[0]) and (eta_max<=bounds_eta[1]):
            global_bins.append(gl_idx)
            bounds_pt_idx.append(idx_pt)
            bounds_eta_idx.append(idx_eta)
            
    return global_bins, [min(bounds_pt_idx), max(bounds_pt_idx)], [min(bounds_eta_idx), max(bounds_eta_idx)]


def cross_section_bkg(bkg_process):
    """
    """
    return xsec_bkg[bkg_process]


def bkg_lumi_scales(type_eff, bkg_categories):
    """
    "Lumi scale" defined as alpha that satisifies lumi_data=alpha*lumi_bkg
    """

    lumi_data = 16.8  # fb^-1

    lumi_scales = {}

    for cat in bkg_categories:
        file = ROOT.TFile(f"/scratchnvme/rajarshi/Bkg_TNP_3D_Histograms/OS/tnp_{type_eff}_{cat}_vertexWeights1_oscharge1.root")
        
        wsum_histo = file.Get("weightSum")
        num_init = wsum_histo.Integral()
        xsection = cross_section_bkg(cat)*1000
        lumi_bkg = num_init/xsection

        scale = lumi_data/lumi_bkg

        lumi_scales.update({cat : scale})
    
    return lumi_scales





if __name__ == '__main__':
    
    a, b, c = get_indexes_from_bounds([24, 28], [-2.4, -1.3])

    print(len(a))
    print(b[0], b[1], c[0], c[1])