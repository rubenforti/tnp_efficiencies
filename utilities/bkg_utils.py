"""
"""
import sys
from array import array
from copy import copy as cp
import ROOT
import utilities.base_lib as base_lib
from utilities.binning_utils import get_pt_binning_ref, bin_dictionary, bin_global_idx_dict
from utilities.results_utils import create_histograms, init_pass_fail_histos
from utilities.dataset_utils import get_totbkg_roohist
from utilities.plot_utils import plot_bkg, plot_2d_bkg_distrib, plot_projected_bkg, plot_bkg_comparison



def base_parser_bkg(parser):
    """
    """
    parser.add_argument("-b", "--bkg_categories", type=str, nargs="+", choices=base_lib.bkg_categories+["all"], default=["all"])
    parser.add_argument("--import_bkg_SS", action="store_true", 
                        help="Import the same-sign background")
    parser.add_argument("--lightMode_bkg", action="store_true", 
                        help="Import only the total background histograms in the workspace")
    parser.add_argument("--mergedbins_bkg", type=str, nargs=2, default=["", ""],
                        help="Use different binning (pt, eta) for background w.r.t. signal and data (that are evaluated in standard bins)")

    return parser

###############################################################################


def bkg_dict_for_plots(ws, import_categories, flag, bin_key, bin_pt, bin_eta,
                       import_fit_pars=False, import_fit_pdf_bkg=False, pdf_bkg_shape="expo"):
    """
    Creates a dictionary containing the bkg MC datasets and the useful information related to them.
    """
    axis = ws.var(f"x_{flag}_{bin_key}")

    binning = axis.getBinning("x_binning")
    axis.setBins(binning.numBins(), "plot_binning")

    datasets = {}
    datasets.update({"axis" : axis})

    for cat in import_categories:
        try:
            type_dataset, subs = cat.split("_", 1)
            subs = f"_{subs}"
        except:
            type_dataset, subs = cat, ""

        datasets[cat] = ws.data(f"Minv_{type_dataset}_{flag}_{bin_key}{subs}")

    import_categories_bkgOS = [cat for cat in import_categories if "SS" not in cat]
    import_categories_bkgSS = [cat for cat in import_categories if "SS" in cat]

    if len(import_categories_bkgOS) > 0 and import_categories_bkgOS!=["bkg_SameCharge"]:
        #print("import_categories_bkgOS")
        datasets["bkg_total"] = get_totbkg_roohist([ws, import_categories_bkgOS], flag, axis, bin_key, bin_pt, bin_eta)
    
    if len(import_categories_bkgSS) > 0:
        datasets["bkg_total_SS"] = get_totbkg_roohist([ws, import_categories_bkgSS], flag, axis, bin_key, bin_pt, bin_eta,
                                                      use_SS=True)

    if import_fit_pars:
        res_obj = ws.obj(f"results_{flag}_{bin_key}")
        pars = res_obj.floatParsFinal()
        datasets["fit_pars"] = {"nsig" : pars.find(f"nsig_{flag}_{bin_key}"),
                                "nbkg" : pars.find(f"nbkg_{flag}_{bin_key}")}

    '''
    if import_fit_pdf_bkg:
        pdf_bkg_fit = ws.pdf(f"{pdf_bkg_shape}_bkg_{flag}_{bin_key}")
        if type(pdf_bkg_fit) is ROOT.TObject:
            sys.exit(f"ERROR: bkg pdf not found at bin {bin_key}")

        res_obj = ws.obj(f"results_{flag}_{bin_key}")
        pars = res_obj.floatParsFinal()
        norm_bkg = pars.find(f"nbkg_{flag}_{bin_key}")

        if type(norm_bkg) is not ROOT.RooRealVar:
            norm_bkg = ROOT.RooRealVar(f"nbkg_{flag}_{bin_key}", f"nbkg_{flag}_{bin_key}", 0)

        datasets.update({"pdf_bkg_fit" : {
            "pdf" : pdf_bkg_fit,
            "norm" : norm_bkg}
            })
    '''

    return datasets

###############################################################################


def analyze_bkg(type_eff, ws_filename, binning_pt, binning_eta, bkg_categories, analysis_list, an_opt, output_path=""):
    """
    """
    file = ROOT.TFile(ws_filename, "READ")
    ws = file.Get("w")

    bins_pt, bins_eta = base_lib.binnings[binning_pt], base_lib.binnings[binning_eta]

    bin_dict = bin_dictionary(binning_pt, binning_eta, pt_binning_ref=get_pt_binning_ref(type_eff))

    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## Section with definition of histograms to be filled and other variables

    if "negweighted_bins" in analysis_list:
        negweight_histos = init_pass_fail_histos(["nw_single", "nw_totBkg", "nsum", "nsum_totBkg"], 
                                                 bins_pt=bins_pt, bins_eta=bins_eta)

    if "2D_distrib" in analysis_list:
        bkg_histos_2D = init_pass_fail_histos(bkg_categories, bins_pt=bins_pt, bins_eta=bins_eta)
        isNorm = an_opt["2D_distrib"]["cmp_data"] or an_opt["2D_distrib"]["cmp_signal"] or an_opt["2D_distrib"]["cmp_totBkg"]

    if "cross_cmp" in analysis_list:
        comp_1, comp_2 = an_opt["cross_cmp"]["cmp_cat"]
        cmp_name_1, cmp_name_2 = cp(comp_1), cp(comp_2)
        
        if "total" in comp_1:
            if "SS" in comp_1:
                cmp_name_1 = "mcSS"
            else:
                cmp_name_1 = "mcOS"
        else:
            cmp_name_1 = comp_1.replace("bkg_", "")

        if "total" in comp_2:
            if "SS" in comp_2:
                cmp_name_2 = "mcSS"
            else:
                cmp_name_2 = "mcOS"
        else:
            cmp_name_2 = comp_2.replace("bkg_", "")

        
        print(cmp_name_1, cmp_name_2)

        histos_crosscmp = init_pass_fail_histos([f"cross_cmp_{cmp_name_1}-{cmp_name_2}"], 
                                               bins_pt=bins_pt, bins_eta=bins_eta)
        print(histos_crosscmp.keys())
        import_categories = [in_cat for in_cat in [comp_1, comp_2] if in_cat in base_lib.bkg_categories]
        if "total" in comp_1 or "total" in comp_2:
            import_categories = base_lib.bkg_categories
        if "SS" in comp_1 or "SS" in comp_2:
            import_categories = import_categories + [cat+"_SS" for cat in import_categories if "SameCharge" not in cat]


    ## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ## Loop on the bins and on pass-fail
    
    for bin_key, [gl_idx, bin_pt, bin_eta] in bin_dict.items():

        for flag in ["pass", "fail"]:

            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            if "negweighted_bins" in analysis_list:

                cnt_weight_total, cnt_weight_single, cnt_sum_total, cnt_sum_single = 0, 0, 0, 0

                for bkg_cat in bkg_categories:
                    bkg_data = ws.data(f"Minv_bkg_{flag}_{bin_key}_{bkg_cat.replace('bkg_', '')}")
                    cnt_sum_single = cnt_sum_single+1 if bkg_data.sumEntries() < 0 else cnt_sum_single
                    for idx in range(bkg_data.numEntries()):
                        cnt_weight_single = cnt_weight_single+1 if bkg_data.weight(idx) < 0 else cnt_weight_single

                histo_totBkg = ws.data(f"Minv_bkg_{flag}_{bin_key}_total")
                if type(histo_totBkg) is ROOT.RooDataHist:
                    cnt_sum_total = 1 if histo_totBkg.sumEntries() < 0 else 0
                    for idx in range(histo_totBkg.numEntries()):
                        cnt_weight_total = cnt_weight_total+1 if histo_totBkg.weight(idx) < 0 else cnt_weight_total
                print(negweight_histos.keys())
                negweight_histos[f"nw_single_{flag}_2d"].SetBinContent(bin_pt, bin_eta, cnt_weight_single)
                negweight_histos[f"nw_totBkg_{flag}_2d"].SetBinContent(bin_pt, bin_eta, cnt_weight_total)
                negweight_histos[f"nsum_{flag}_2d"].SetBinContent(bin_pt, bin_eta, cnt_sum_single)
                negweight_histos[f"nsum_totBkg_{flag}_2d"].SetBinContent(bin_pt, bin_eta, cnt_sum_total)
            
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            if "minv_distrib" in analysis_list:

                setlog = an_opt["minv_distrib"]["logscale"]
                if an_opt["minv_distrib"]["logscale"]=="hybrid": 
                    setlog=True if flag=="pass" else False

                ch_sel = "SS" if an_opt["minv_distrib"]["is_SS_bkg"] else "OS"

                if an_opt["minv_distrib"]["cmp_data"] is True:
                    import_categories = ["data"] + bkg_categories
                    figpath_postfix = "_w_data"
                elif an_opt["minv_distrib"]["cmp_signal"] is True:
                    import_categories = ["mc" if not an_opt["minv_distrib"]["is_SS_bkg"] else "mc_SS"] + bkg_categories
                    figpath_postfix = "_w_signal"
                else:
                    import_categories = bkg_categories
                    figpath_postfix = ""


                datasets = bkg_dict_for_plots(ws, import_categories, flag, bin_key, bin_pt, bin_eta)

                plot_bkg(datasets, flag, bin_key, logscale=setlog, figpath=f"{output_path}/minv_plots{figpath_postfix}")

            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            if "2D_distrib" in analysis_list:

                if an_opt["2D_distrib"]["cmp_data"]: 
                    h_norm = ws.data(f"Minv_data_{flag}_{bin_key}")
                elif an_opt["2D_distrib"]["cmp_signal"]:
                    h_norm = ws.data(f"Minv_mc_{flag}_{bin_key}")
                elif an_opt["2D_distrib"]["cmp_totBkg"]:
                    if type(ws.data(f"Minv_bkg_{flag}_{bin_key}_total")) is ROOT.RooDataHist:
                        h_norm = ws.data(f"Minv_bkg_{flag}_{bin_key}_total")
                    else:
                        print("Warning! Total background not found in the workspace, normalization results impossible")
                        isNorm = False
                        '''
                        # Seems superflous
                        axis = ws.var(f"x_{flag}_{bin_key}")
                        if binning_pt not in ["pt", "pt_tracking"] or binning_eta!="eta":
                            extended_bin_dict = bin_global_idx_dict(binning_pt, binning_eta)
                            bin_pt_ext, bin_eta_ext = extended_bin_dict[str(gl_idx[0])]
                        else:
                            bin_pt_ext, bin_eta_ext = bin_pt, bin_eta
                        h_norm = get_totbkg_roohist([ws, bkg_categories], flag, axis, bin_key, bin_pt_ext, bin_eta_ext)
                        '''
                else:
                    pass

                for bkg_cat in bkg_categories:
                    type_dset, subs = ["bkg", bkg_cat] if bkg_cat != "mc_SS" else ["mc", "SS"]

                    print(type_dset, subs)
                    
                    bkg_hist = ws.data(f"Minv_{type_dset}_{flag}_{bin_key}_{subs.replace('bkg_', '')}")

                    n_events, err_events = bkg_hist.sumEntries(), base_lib.sumw2_error(bkg_hist)

                    if isNorm:
                        if h_norm.sumEntries() > 0 and n_events > 0:
                            if h_norm.sumEntries() < n_events: print(n_events, h_norm.sumEntries())
                            n_events = n_events/h_norm.sumEntries()
                            err_events = n_events*( (err_events/(n_events*h_norm.sumEntries()))**2 + 
                                                    (base_lib.sumw2_error(h_norm)/h_norm.sumEntries())**2 )**0.5
                        else:
                            n_events, err_events = 0, 0
                    else:
                        pass

                    bkg_histos_2D[f"{bkg_cat}_{flag}_2d"].SetBinContent(bin_pt, bin_eta,n_events)
                    bkg_histos_2D[f"{bkg_cat}_{flag}_2d"].SetBinError(bin_pt, bin_eta, err_events)
                    print(bkg_histos_2D[f"{bkg_cat}_{flag}_2d"].Integral())
            
            # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            
            if "cross_cmp" in analysis_list:
                
                plot_dict = bkg_dict_for_plots(ws, import_categories, flag, bin_key, bin_pt, bin_eta)
                plot_dict = {k : v for k, v in plot_dict.items() if k.replace(f"_{flag}_{bin_key}", "") in ["axis", comp_1, comp_2]}
                histos_crosscmp[
                    f"cross_cmp_{cmp_name_1}-{cmp_name_2}_{flag}_2d"].SetBinContent(
                        bin_pt, bin_eta, list(plot_dict.values())[1].sumEntries()/list(plot_dict.values())[2].sumEntries())
                
                plot_bkg_comparison(plot_dict, flag, bin_key, figpath=f"{output_path}/crossCmp_{cmp_name_1}-{cmp_name_2}")
                

    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Stuff to be done when the loop has terminated

    

    if "negweighted_bins" in analysis_list:
        file_out = ROOT.TFile(f"{output_path}/prob_bins_pt{len(bins_pt)-1}_eta{len(bins_eta)-1}.root", "RECREATE")
        file_out.cd()
        [histo.Write() for histo in negweight_histos.values()]
        file_out.Close()


    if "2D_distrib" in analysis_list:
        print(bkg_histos_2D.values())
        for bkg_cat in bkg_categories:
            print(bkg_cat, bkg_histos_2D[f"{bkg_cat}_pass_2d"])
            hist_dict_bkgcat = {"pass" : bkg_histos_2D[f"{bkg_cat}_pass_2d"], "fail" : bkg_histos_2D[f"{bkg_cat}_fail_2d"]}
            print(type(hist_dict_bkgcat["pass"]))
            plot_2d_bkg_distrib(hist_dict_bkgcat, bkg_cat, figpath=f"{output_path}/pt-eta_distrib")

        file_out = ROOT.TFile(f"{output_path}/bkg_2d_distrib.root", "RECREATE")
        file_out.cd()
        [histo.Write() for histo in bkg_histos_2D.values()]
        file_out.Close()

        if an_opt["2D_distrib"]["projected"]: 
            for flag in ["pass", "fail"]:
                plot_projected_bkg({k : v for k, v in bkg_histos_2D.items() if flag in k}, 
                                   "pt_singlebin", binning_eta, flag, figpath=f"{output_path}/pt-eta_distrib")
                plot_projected_bkg({k : v for k, v in bkg_histos_2D.items() if flag in k}, 
                                   binning_pt, "eta_singlebin", flag, figpath=f"{output_path}/pt-eta_distrib")
                
    if "cross_cmp" in analysis_list:
        file_out = ROOT.TFile(f"{output_path}/crossCmp_{cmp_name_1}-{cmp_name_2}.root", "RECREATE")
        file_out.cd()
        [histo.Write() for histo in histos_crosscmp.values()]
        file_out.Close()

    file.Close()