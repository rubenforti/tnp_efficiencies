"""
"""

import sys
import os
import ROOT
from copy import deepcopy
from utilities.base_library import binning, bin_dictionary, lumi_factors, lumi_factor, get_idx_from_bounds, bin_global_idx_dict


bkg_sum_selector = {
    "WW" : 1,
    "WZ" : 1,
    "ZZ" : 1,
    "TTFullyleptonic" : 1,
    "TTSemileptonic" : 0,
    "Ztautau" : 1,
    "SameCharge" : 1,
    "WJets" : 0,
    "signalMC_SS" : -1.
}

def gen_import_dictionary(base_folder, type_eff, dset_names,
                          ch_set = [""],
                          scale_MC = False,
                          add_SS_mc = False,
                          add_SS_bkg = False):
    """
    """
    import_dictionary = {}

    for dset_name in dset_names:

        import_dictionary[dset_name] = {}
        import_dictionary[dset_name]["filenames"] = []
        
        lumi_scale = -1

        S_sel, c_sel = "OS", "1"

        flag_set = dset_name.replace("bkg_", "")  # If it is not a background, it simply copies the dataset name
        
        if "SameCharge" in dset_name: 
            flag_set = "data"
            S_sel, c_sel = "SS", "0"

        for ch in ch_set:
                
                filename = f"{base_folder}/{S_sel}/tnp_{type_eff}{ch}_{flag_set}_vertexWeights1_oscharge{c_sel}.root"
                import_dictionary[dset_name]["filenames"].append(filename)
        
        print(dset_name)
        if (dset_name=="mc" and scale_MC is True) or ("bkg" in dset_name and "SameCharge" not in dset_name): 
            lumi_scale = lumi_factor(filename, flag_set)
        
        import_dictionary[dset_name]["lumi_scale"] = lumi_scale

        if dset_name=="mc" and add_SS_mc is True:
            import_dictionary["mc_SS"] = deepcopy(import_dictionary["mc"])
            import_dictionary["mc_SS"]["filenames"] = [f.replace("OS", "SS").replace("oscharge1", "oscharge0") 
                                                       for f in import_dictionary["mc"]["filenames"]]
            import_dictionary["mc_SS"]["lumi_scale"] = import_dictionary["mc"]["lumi_scale"]
        
        if "bkg" in dset_name and "SameCharge" not in dset_name and add_SS_bkg is True:
            import_dictionary[f"bkg_{flag_set}_SS"] = deepcopy(import_dictionary[f"bkg_{flag_set}"])
            import_dictionary[f"bkg_{flag_set}_SS"]["filenames"] = [f.replace("OS", "SS").replace("oscharge1", "oscharge0")
                                                                    for f in import_dictionary[f"bkg_{flag_set}"]["filenames"]]
            import_dictionary[f"bkg_{flag_set}_SS"]["lumi_scale"] = import_dictionary[f"bkg_{flag_set}"]["lumi_scale"]

    return import_dictionary
                


def import_pdf_library(*functions):
    """
    Imports the C++ defined functions into the ROOT interpreter. 
    """
    current_path = os.path.dirname(__file__)
    import_path = os.path.join(current_path, '..', 'libCpp')
    
    for function in functions:
        ctrl_head = ROOT.gInterpreter.Declare(f' #include "{import_path}/{function}.h"')
        ctrl_source = ROOT.gSystem.CompileMacro(f"{import_path}/{function}.cc", opt="ks")

        if ctrl_head is not True:
            print("ERROR in header loading")
            sys.exit()
        if ctrl_source != 1:
            print("ERROR in sourcefile compiling and loading")
            sys.exit()


###############################################################################

def import_totbkg_hist(ws, bin_key, bkg_categories, bkg_selector=bkg_sum_selector):
    """
    Imports into a workspace a RooDataHist that represents the total 
    background, for a given pt-eta bin. The background sources need to be 
    correctly normalized before calling this function, and their impact in 
    the total bkg histogram is managed by the "bkg_selector" dictionary.
    """
    for flag in ["pass", "fail"]:
        axis = ws.var(f"x_{flag}_{bin_key}")
        bkg_total_histo = ROOT.RooDataHist(f"Minv_bkg_{flag}_{bin_key}_total", "bkg_total_histo", 
                                       ROOT.RooArgSet(axis), "x_binning")

        for bkg_cat in bkg_categories:
            bkg_dataset = ws.data(f"Minv_bkg_{flag}_{bin_key}_{bkg_cat}")
            if bkg_selector[bkg_cat] == 1: 
                bkg_total_histo.add(bkg_dataset)
            elif bkg_selector[bkg_cat] == -1:
                bkg_dataset_copy = bkg_dataset.Clone(f"Minv_bkg_{flag}_{bin_key}_{bkg_cat}_copy")
                for i in range(axis.getBinning("x_binning").numBins()):
                    bkg_dataset_copy.get(i)
                    bkg_dataset_copy.set(i, -1*bkg_dataset_copy.weight(i), bkg_dataset_copy.weightError(ROOT.RooAbsData.SumW2))
                bkg_total_histo.add(bkg_dataset_copy)
            else:
                pass

        ws.Import(bkg_total_histo)


###############################################################################

def get_roohist(files, type_set, flag, axis, bin_key, bin_pt, bin_eta, 
                global_scale=-1.):
    """
    Returns a RooDataHists of the variable TP_mass, in a single pt-eta bin.
    For MC datasets, the "global_scale" parameter is used to scale the
    histogram to the correct luminosity.
    """
    if type_set=="data" or type_set=="data_SC":
        type_suffix = "RunGtoH"
    elif type_set=="mc" or type_set=="mc_w" or type_set=="bkg":
        type_suffix = "DY_postVFP"
    else:
        sys.exit("ERROR: invalid category type given")

    if type(bin_pt) is int:  bin_pt = [bin_pt, bin_pt]
    if type(bin_eta) is int: bin_eta = [bin_eta, bin_eta]    

    numBins = axis.getBinning("x_binning").numBins()

    n_files = len(files) if type(files) is list else 1

    if type_set=="data_SC": type_set = "bkg"

    roohisto = ROOT.RooDataHist(f"Minv_{type_set}_{flag}_{bin_key}", f"Minv_{type_set}_{flag}_{bin_key}",
                                ROOT.RooArgSet(axis), "x_binning")
    
    print(type_suffix)

    tmp_histos, tmp_roohistos = [], []

    print(f"\nImporting histograms in bin {bin_key}")

    for it_file in range(n_files):
        
        file = files[it_file] if n_files > 1 else files

        print(f"{flag}_mu_{type_suffix}")
        histo3d = file.Get(f"{flag}_mu_{type_suffix}")

        # In the projection, option "e" is specified to calculate the bin 
        # content errors in the new histogram for generic selection of bin_pt
        # and bin_eta. Without it, all works as expected ONLY IF the projection
        # is done on a single bin of pt-eta
        th1_histo = histo3d.ProjectionX(f"Histo_{type_set}_{flag}", bin_pt[0], bin_pt[1], bin_eta[0], bin_eta[1], "e")

        if global_scale > 0: th1_histo.Scale(global_scale)


        if th1_histo.Integral() < 0:
            print(f"ERROR: negative entries in TH1 in bin {bin_pt}-{bin_eta}")
            #sys.exit()

        th1_histo.Rebin(int(th1_histo.GetNbinsX()/numBins))

        tmp_histos.append(th1_histo.Clone(f"Minv_{type_set}_{flag}_{bin_key}_tmp"))

        tmp_roohistos.append(ROOT.RooDataHist(f"Minv_{type_set}_{flag}_{bin_key}", f"Minv_{type_set}_{flag}_{bin_key}",
                                      ROOT.RooArgList(axis), th1_histo))

        roohisto.add(tmp_roohistos[it_file])

    print("Beginning consistency check...")
    th1_integral = 0
    for it_file in range(n_files):
        th1_integral += tmp_histos[it_file].Integral()
    if round(th1_integral, 5) != round(roohisto.sumEntries(), 5):
        sys.exit(f"ERROR: TH1 and RooDataHist have different integrals in bin {bin_key}")
    for i in range(numBins):
        roohisto.get(i)
        th1_bincontent = 0
        for it_file in range(n_files): 
            th1_bincontent += tmp_histos[it_file].GetBinContent(i+1)
        if round(roohisto.weight(i), 5) != round(th1_bincontent, 5):
            # print(roohisto.weight(i), th1_bincontent)
            sys.exit(f"ERROR: bin {i} has different content in TH1 and RooDataHist")
    print("Consistency check passed, RooDataHist imported correctly\n")
            
    return roohisto


###############################################################################

def ws_init(import_datasets, type_analysis, binning_pt, binning_eta, binning_mass, 
            import_existing_ws = False, 
            existing_ws_filename = "", 
            lightMode_bkg = False, 
            altBinning_bkg = False):
    """
    Initializes a RooWorkspace with the datasets corresponding to a given 
    efficiency step. The objects stored are RooDataHist of TP_mass in every 
    pt-eta bin requested. The import_sets dictionary has as keys the names of 
    the datasets to be imported ("data", "mc", and the various backgrounds) and
    must contain as values:
        - the paths of the file;
        - the luminosity scale factors (for MC datasets);
    The value of a MC dataset has to be a dictionary with the "filename" and 
    "lumi_scale" keys. For "mc" case (referred to the signal template), it is
    possible not to specify the luminosity scale factor, since it could be used
    as a mere template: furhter consistency checks could be necessary in this
    case, depending on the analysis requested. 
    """

    if altBinning_bkg is False:
        bins = bin_dictionary(binning_pt, binning_eta)
    else:
        bin_idx_dict = bin_global_idx_dict(binning_pt, binning_eta)
        bins = bin_dictionary()

    bins_mass = binning(binning_mass)

    if import_existing_ws is True:
        ws_file = ROOT.TFile(existing_ws_filename)
        ws = ws_file.Get("w")
    else:
        ws = ROOT.RooWorkspace("w")
        x = ROOT.RooRealVar("x", "TP M_inv", bins_mass[0], bins_mass[-1], unit="GeV/c^2")
        x_binning = ROOT.RooUniformBinning(bins_mass[0], bins_mass[-1], len(bins_mass)-1, "x_binning")
        ws.Import(x)

    for dataset_obj in import_datasets.values():
        dataset_obj["filenames"] = [ROOT.TFile(v, "READ") for v in dataset_obj["filenames"]]

    if lightMode_bkg is True:
        # Importing only the total background histogram
        bkg_merging_dict = {k : v for k, v in import_datasets.items() if k.startswith("bkg_")}
        import_datasets  = {k : v for k, v in import_datasets.items() if not k.startswith("bkg_")}

    
    # Loop over the pt-eta bins
    for bin_key in bins.keys():

        print(f"\n\n### Importing datasets in bin {bin_key} ###\n")
        gl_idx, bin_pt, bin_eta = bins[bin_key]

        flags = ["pass", "fail"] if type_analysis == "indep" else ["sim"]
        for flag in flags:

            if import_existing_ws is False:
                axis = ROOT.RooRealVar(f"x_{flag}_{bin_key}", "TP M_inv", bins_mass[0], bins_mass[-1], unit="GeV/c^2")
                axis.setBinning(x_binning)
            else:
                axis = ws.var(f"x_{flag}_{bin_key}")

            # Importing the datasets
            for dset_name, dset_obj in import_datasets.items():

                files, lumi_scale = dset_obj["filenames"], dset_obj["lumi_scale"]

                if "bkg" in dset_name:
                    dset_type = "bkg" if not "SameCharge" in dset_name else "data_SC"
                    dset_subs = dset_name.replace("bkg", "")
                    if altBinning_bkg is True: bin_pt, bin_eta = bin_idx_dict[str(gl_idx)]
                elif "SS" in dset_name:
                    dset_type = dset_name.replace("_SS", "")
                    dset_subs = "_SS"
                else:
                    dset_type = dset_name
                    dset_subs = ""

                print(dset_name, dset_type, dset_subs)  

                histo = get_roohist(files, dset_type, flag, axis, bin_key, bin_pt, bin_eta, global_scale=lumi_scale)
                histo.SetName(f"{histo.GetName()}{dset_subs}")
                ws.Import(histo)
            
            # Importing the total background histogram if light mode is called (only for OS)
            if lightMode_bkg is True:
                bkg_total_histo = ROOT.RooDataHist(f"Minv_bkg_{flag}_{bin_key}_total", f"Minv_bkg_{flag}_{bin_key}_total", 
                                                   ROOT.RooArgSet(axis), "x_binning")
                print(bkg_merging_dict.keys())
                for bkg_cat, bkg_obj in bkg_merging_dict.items():
                    print(bkg_cat)
                    if "SS" not in bkg_cat and bool(bkg_sum_selector[bkg_cat.replace("bkg_","")]) is True:
                        files = bkg_obj["filenames"]
                        print(files)
                        lumi_scale = bkg_obj["lumi_scale"]
                        if altBinning_bkg is True: bin_pt, bin_eta = bin_idx_dict[str(gl_idx)]
                        dset_type = "bkg" if not "SameCharge" in bkg_cat else "data_SC"
                        histo = get_roohist(files, dset_type, flag, axis, bin_key, bin_pt, bin_eta, global_scale=lumi_scale)
                        bkg_total_histo.add(histo)
                    else:
                        pass
                ws.Import(bkg_total_histo)

    return ws


      
###############################################################################     
###############################################################################


if __name__ == '__main__':

    type_eff = "tracking"

    base_folder = "/scratchnvme/rajarshi/Latest_3D_Steve_Histograms_22_Sep_2023"

    d = gen_import_dictionary(base_folder, type_eff, ["data", "mc", "bkg_WW", "bkg_WZ", "bkg_SameCharge"], 
                              ch_set=["plus", "minus"], scale_MC=True, add_SS_mc=False, add_SS_bkg=False)
    

    ws= ws_init(d, "indep", "pt_singlebin", "eta_singlebin", "mass_50_130", lightMode_bkg=True)

    ws.Print()


