"""
"""
import utilities.base_lib as base_lib


def bin_dictionary(binning_pt_name="pt", binning_eta_name="eta", get_mergedbins_bounds=False):
    """
    Creates a dictionary that relates the physical bounds of every bin (keys)
    to the indexes useful for the bin selection. The indexes are returned in a 
    3-element list containing the global index, the pt index and the eta index
    in this order. If the bins are not merged (w.r.t. the original binning of a
    given step) these objects are integers i; otherwise, two cases can be chosen
    by the "get_mergedbins_bounds" flag:
      - get_mergedbins_bounds=True: the pt and eta indexes are returned as
        lists, containing the indexes of bins in the original binning that are
        merged in the given bin;
      - get_mergedbins_bounds=False: the pt and eta indexes are returned as
        integers, containing the index referred to the actual grid pt-eta.
    """
    index_dictionary = {}
    global_idx = 1

    binning_pt, binning_eta = base_lib.binnings[binning_pt_name], base_lib.binnings[binning_eta_name]

    if "tracking" in binning_pt_name: binning_pt_name = "pt"

    cnt_mergedpt = 0
    nbins_eta = len(binning_eta)-1

    for idx_pt in range(1, len(binning_pt)):
        for idx_eta in range(1, len(binning_eta)):

            bin_key = f"[{binning_pt[idx_pt-1]}to{binning_pt[idx_pt]}][{binning_eta[idx_eta-1]}to{binning_eta[idx_eta]}]"

            if binning_pt_name == "pt" and binning_eta_name == "eta":
                index_dictionary[bin_key] = [global_idx, idx_pt, idx_eta]
                global_idx +=1
            else:
                global_idx, bounds_idx_pt, bounds_idx_eta = get_idx_from_bounds([binning_pt[idx_pt-1], binning_pt[idx_pt]],
                                                                                [binning_eta[idx_eta-1], binning_eta[idx_eta]])
                if get_mergedbins_bounds is False:
                    eta_idx = int(1+(nbins_eta*(bounds_idx_eta[0]-1)/48.))
                    pt_idx = int(bounds_idx_pt[0] - cnt_mergedpt)
                    cnt_mergedpt += bounds_idx_pt[-1]-bounds_idx_pt[0] if eta_idx==nbins_eta else 0
                    index_dictionary[bin_key] = [global_idx, pt_idx, eta_idx]
                else:
                    index_dictionary[bin_key] = [global_idx, bounds_idx_pt, bounds_idx_eta]

    return index_dictionary

###############################################################################


def get_idx_from_bounds(bounds_pt, bounds_eta):
    """
    Function that, given the bin bounds (that have to be present in the binning
    arrays!!!), returns the coordinates of the selected region. The coordinates 
    are given as list of global indexes and as bounds on pt/eta indexes (max 
    and min, since the region is rectangular)
    """
    initial_dict = bin_dictionary()
    global_bins = []
    bounds_pt_idx, bounds_eta_idx = [], []

    for el, [gl_idx, idx_pt, idx_eta] in initial_dict.items():

        string_pt, string_eta = el.split("][")
        pt_min, pt_max = string_pt[1:].split("to")
        eta_min, eta_max = string_eta[:-1].split("to")

        pt_min, pt_max = float(pt_min), float(pt_max)
        eta_min, eta_max = float(eta_min), float(eta_max)

        if(pt_min>=bounds_pt[0]) and (pt_max<=bounds_pt[1]) and (eta_min>=bounds_eta[0]) and (eta_max<=bounds_eta[1]):
            global_bins.append(gl_idx)
            bounds_pt_idx.append(idx_pt)
            bounds_eta_idx.append(idx_eta)
            
    return global_bins, [min(bounds_pt_idx), max(bounds_pt_idx)], [min(bounds_eta_idx), max(bounds_eta_idx)]

###############################################################################

def bin_global_idx_dict(binning_pt, binning_eta):
    """
    Returns a dictionary that relates every global index to the pt-eta indexes
    in the given binning.
    """
    bin_dict = bin_dictionary(binning_pt, binning_eta, get_mergedbins_bounds=True)

    bin_idx_dict = {}

    for bin_key, [gl_idx, pt_idx, eta_idx] in bin_dict.items():
        if type(gl_idx) is list:
            [bin_idx_dict.update({str(gl) : [pt_idx, eta_idx]}) for gl in gl_idx]
        else:
            bin_idx_dict.update({str(gl_idx) : [pt_idx, eta_idx]})
            
    return bin_idx_dict