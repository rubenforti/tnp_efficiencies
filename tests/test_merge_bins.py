"""
"""

import ROOT
import unittest
from utilities.base_library import binning, get_new_binning, bin_dict


def create_flat_th3():
    std_binning_pt = binning("pt")
    std_binning_eta = binning("eta")
    mass_binning = binning("mass_60_120")

    histo = ROOT.TH3D("histo", "histo", len(mass_binning)-1, mass_binning,
                      len(std_binning_pt)-1, std_binning_pt, len(std_binning_eta)-1, std_binning_eta)
        
    nbins = (len(std_binning_pt)-1)*(len(std_binning_eta)-1)*(len(mass_binning)-1)
    
    for i in range(1, histo.GetNbinsX()+1):
        for j in range(1, histo.GetNbinsY()+1):
            for k in range(1, histo.GetNbinsZ()+1):
                histo.SetBinContent(i, j, k, 1)

    return histo


class TestMergeBins(unittest.TestCase):
    

    def test_coverage_merge_pt(self):

        newbins_pt = binning("pt_6bins")
        newbins_eta = binning("eta")

        new_bin_dictionary = get_new_binning(newbins_pt, newbins_eta)

        global_bins_check = []

        for bin_key in new_bin_dictionary:
            global_bins, _, _ = new_bin_dictionary[bin_key]

            for gl_bin in global_bins:
                global_bins_check.append(gl_bin) if gl_bin not in global_bins_check else self.assertEqual(1,0)
            
        global_bins_check.sort()
        ordered_bins = list(range(1, 721))

        self.assertTrue(global_bins_check == ordered_bins)

    
    def test_coverage_merge_eta(self):

        newbins_pt = binning("pt")
        newbins_eta = binning("eta_8bins")

        new_bin_dictionary = get_new_binning(newbins_pt, newbins_eta)

        global_bins_check = []

        for bin_key in new_bin_dictionary:
            global_bins, _, _ = new_bin_dictionary[bin_key]

            for gl_bin in global_bins:
                global_bins_check.append(gl_bin) if gl_bin not in global_bins_check else self.assertEqual(1,0)
            
        global_bins_check.sort()
        ordered_bins = list(range(1, 721))

        self.assertTrue(global_bins_check == ordered_bins)        
    

    def test_coverage_merge_both(self):

        newbins_pt = binning("pt_9bins")
        newbins_eta = binning("eta_16bins")

        new_bin_dictionary = get_new_binning(newbins_pt, newbins_eta)

        global_bins_check = []

        for bin_key in new_bin_dictionary:
            global_bins, _, _ = new_bin_dictionary[bin_key]

            for gl_bin in global_bins:
                global_bins_check.append(gl_bin) if gl_bin not in global_bins_check else self.assertEqual(1,0)
            
        global_bins_check.sort()
        ordered_bins = list(range(1, 721))

        self.assertTrue(global_bins_check == ordered_bins)


    def test_integral_merged_bins(self):

        histo = create_flat_th3()

        newbins_pt = binning("pt_6bins")
        newbins_eta = binning("eta_16bins")

        new_bin_dictionary = get_new_binning(newbins_pt, newbins_eta)

        axis = ROOT.RooRealVar("x", "x", 60, 120)
        axis.setRange("range", 60, 120)
        axis.setBins(60)

        for bin_key in new_bin_dictionary:

            _, bounds_pt_idx, bounds_eta_idx = new_bin_dictionary[bin_key]

            h1d = histo.ProjectionX("h1d", bounds_pt_idx[0], bounds_pt_idx[1], bounds_eta_idx[0], bounds_eta_idx[1], "e")

            roohisto = ROOT.RooDataHist("roohisto", "roohisto", ROOT.RooArgList(axis), h1d)

            self.assertEqual(roohisto.sumEntries(), 60*(bounds_pt_idx[1]-bounds_pt_idx[0]+1)*(bounds_eta_idx[1]-bounds_eta_idx[0]+1))  






        
 







if __name__ == "__main__":
    unittest.main()