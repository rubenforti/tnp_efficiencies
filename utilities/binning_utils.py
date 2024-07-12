"""
"""
import sys
import utilities.base_lib as base_lib


def get_pt_binning_ref(type_eff):
    """
    """
    if type_eff in ["reco", "tracking"]:
        return f"pt_{type_eff}"
    elif type_eff in ["idip", "trigger", "iso"]:
        return "pt"
    else:
        sys.exit("Wrong efficiency step inserted. Exiting...")


def bin_dictionary(binning_pt_name, binning_eta_name, get_mergedbins_bounds=False, pt_binning_ref="pt"):
    """
    Creates a dictionary that relates the physical bounds of every bin (keys)
    to the indexes useful for the bin selection. The indexes are returned in a 
    3-element list containing the global index, the pt index and the eta index
    in this order. If the bins are not merged (w.r.t. the original binning of a
    given step) these objects are integers i; otherwise, two cases can be 
    chosen by the "get_mergedbins_bounds" flag:
      - get_mergedbins_bounds=True: the pt and eta indexes are returned as
        lists, containing the indexes of bins in the original binning that are
        merged in the given bin;
      - get_mergedbins_bounds=False: the pt and eta indexes are returned as
        integers, containing the index referred to the actual grid pt-eta.
    
    """
    index_dictionary = {}
    global_idx = 1

    binning_pt, binning_eta = base_lib.binnings[binning_pt_name], base_lib.binnings[binning_eta_name]

    # The special binnings for reco and tracking are not to be intended as
    # "merged" binnings: therefore, the binning name is changed to "pt" to
    # enter the regular loop that returns the simple indexes, instead of
    if "reco" in binning_pt_name or "tracking" in binning_pt_name: binning_pt_name = "pt"
        
    cnt_mergedpt = 0
    nbins_eta = len(binning_eta)-1    
    
    for idx_pt in range(1, len(binning_pt)):
        for idx_eta in range(1, len(binning_eta)):

            bin_key = f"[{binning_pt[idx_pt-1]}to{binning_pt[idx_pt]}][{binning_eta[idx_eta-1]}to{binning_eta[idx_eta]}]"

            if binning_pt_name == "pt" and binning_eta_name == "eta":
                index_dictionary[bin_key] = [global_idx, idx_pt, idx_eta]
                global_idx +=1
                # print(global_idx, idx_pt, idx_eta)
            else:
                global_idx, bounds_idx_pt, bounds_idx_eta = get_idx_from_bounds([binning_pt[idx_pt-1], binning_pt[idx_pt]],
                                                                                [binning_eta[idx_eta-1], binning_eta[idx_eta]],
                                                                                pt_binning_ref=pt_binning_ref)
                if get_mergedbins_bounds is False:
                    eta_idx = int(1+(nbins_eta*(bounds_idx_eta[0]-1)/48.))
                    pt_idx = int(bounds_idx_pt[0] - cnt_mergedpt)
                    cnt_mergedpt += bounds_idx_pt[-1]-bounds_idx_pt[0] if eta_idx==nbins_eta else 0
                    index_dictionary[bin_key] = [global_idx, pt_idx, eta_idx]
                else:
                    index_dictionary[bin_key] = [global_idx, bounds_idx_pt, bounds_idx_eta]

                # print(global_idx, bounds_idx_pt, bounds_idx_eta)

    return index_dictionary

###############################################################################


def get_idx_from_bounds(bounds_pt, bounds_eta, pt_binning_ref="pt"):
    """
    Function that, given the "physical" bin bounds, returns the coordinates of
    the corresponding region in the reference system of the standard binning. 
    The coordinates are given as lists indexes (note that every region is 
    rectangular). The "pt_binning_ref" parameter is used to specify the pt 
    binning to be used as standard, needed for 'reco' and 'tracking' steps.
    """

    initial_dict = bin_dictionary(binning_pt_name=pt_binning_ref, binning_eta_name="eta")
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
