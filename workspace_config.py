"""
"""

import ROOT
import sys
from array import array





if __name__ == '__main__':

    type_eff = ("sa", "global", "ID", "iso", "trigger", "veto")
    t = type_eff[3]

    types_analysis = ["indep", "sim"]
    an = types_analysis[1]

    binning_pt = array('d', [24., 26., 28., 30., 32., 34., 36., 38., 40., 42., 44., 47., 50., 55., 60., 65.])

    binning_eta = array('d', [round(-2.4 + i*0.1, 2) for i in range(49)])

    binning_mass = array('d', [60 + i for i in range(61)])

    filename_data = "/scratchnvme/wmass/Steve_root_files/Standard_SF_files/tnp_iso_data_vertexWeights1_oscharge1.root"
    filename_mc = "/scratchnvme/wmass/Steve_root_files/Standard_SF_files/tnp_iso_mc_vertexWeights1_oscharge1.root"

    f_data = ROOT.TFile(filename_data)
    f_mc = ROOT.TFile(filename_mc)
    print(type(f_data))
    print(type(f_mc))

    w = ws_init(f_data, f_mc, t, an, [1,1], [1,1], binning_mass)

    w.writeToFile(f"root_files/ws/ws_{t}_{an}_new.root")


