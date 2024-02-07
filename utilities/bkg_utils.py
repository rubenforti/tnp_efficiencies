"""
"""
import ROOT
import sys
from utilities.base_library import binning, bin_dictionary, bin_global_idx_dict, sumw2_error
from utilities.base_library import eval_efficiency as eval_bkgfrac  #La formula è la stessa
from utilities.results_utils import init_pass_fail_histos
from utilities.dataset_utils import get_totbkg_roohist
from utilities.plot_utils import plot_bkg, plot_2d_bkg_distrib, plot_projected_bkg
from array import array


# bkg_categories = ["WW", "WZ", "ZZ", "TTFullyleptonic", "Ztautau"]

###############################################################################

def show_negweighted_bins(ws_filename, bkg_categories, binning_pt, binning_eta, filepath="bkg_studies"):
    """
    """

    file = ROOT.TFile(ws_filename, "READ")
    ws = file.Get("w")

    bins_pt, bins_eta = binning(binning_pt), binning(binning_eta)
    nbins_pt, nbins_eta = len(bins_pt)-1, len(bins_eta)-1

    bin_dict = bin_dictionary(binning_pt, binning_eta)

    neg_weight_single = ROOT.TH2D("neg_weight_single", "neg_weight_single", 
                                  nbins_pt, bins_pt, nbins_eta, bins_eta)
    neg_weight_total = ROOT.TH2D("neg_weight_total", "neg_weight_total", 
                                 nbins_pt, bins_pt, nbins_eta, bins_eta)

    neg_sumEntries_single = ROOT.TH2D("neg_sumEntries_single", "neg_sumEntries_single", 
                                      nbins_pt, bins_pt, nbins_eta, bins_eta)
    neg_sumEntries_total = ROOT.TH2D("neg_sumEntries_total", "neg_sumEntries_total", 
                                     nbins_pt, bins_pt, nbins_eta, bins_eta)

    for bin_key, [gl_idx, bin_pt, bin_eta] in bin_dict.items():

        print(bin_key)

        cnt_sum_total, cnt_sum_single, cnt_weight_total, cnt_weight_single = 0, 0, 0, 0

        for flag in ["pass", "fail"]:
            axis = ws.var(f"x_{flag}_{bin_key}")
            tot_histo = ROOT.RooDataHist("tot_histo", "tot_histo", ROOT.RooArgSet(axis), "")

            for bkg_cat in bkg_categories:
                bkg_data = ws.data(f"Minv_bkg_{flag}_{bin_key}_{bkg_cat}")
                tot_histo.add(bkg_data) 
                cnt_sum_single = cnt_sum_single+1 if bkg_data.sumEntries() < 0 else cnt_sum_single
                    
            cnt_sum_total = cnt_sum_total+1 if tot_histo.sumEntries() < 0 else cnt_sum_total

            for idx in range(tot_histo.numEntries()):
                cnt_weight_total = cnt_weight_total+1 if tot_histo.weight(idx) < 0 else cnt_weight_total
                for bkg_cat in bkg_categories:
                    histo = ws.data(f"Minv_bkg_{flag}_{bin_key}_{bkg_cat}")
                    cnt_weight_single = cnt_weight_single+1 if histo.weight(idx) < 0 else cnt_weight_single

        neg_weight_single.SetBinContent(bin_pt, bin_eta, cnt_weight_single)
        neg_weight_total.SetBinContent(bin_pt, bin_eta, cnt_weight_total)
        neg_sumEntries_single.SetBinContent(bin_pt, bin_eta, cnt_sum_single)
        neg_sumEntries_total.SetBinContent(bin_pt, bin_eta, cnt_sum_total)    
            
    
    file = ROOT.TFile(f"{filepath}/prob_bins_pt{nbins_pt}_eta{nbins_eta}.root", "RECREATE")
    file.cd()
    neg_weight_single.Write()
    neg_weight_total.Write()
    neg_sumEntries_single.Write()
    neg_sumEntries_total.Write()
    file.Close()


###############################################################################

def make_bkg_dictionary(ws, import_categories, flag, bin_key, bin_pt, bin_eta,
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

    datasets["bkg_total"] = get_totbkg_roohist([ws, import_categories], flag, axis, bin_key, bin_pt, bin_eta)

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

def bkg_mass_distribution(type_eff, ws_filename, bkg_categories, binning_pt, binning_eta, 
                          study_SS_bkg=False,                            
                          plot_on_data=False, plot_fit_bkgpdf=False,
                          plot_on_signal=False, compare_bkgfrac=False,
                          logscale=True, figpath='figs/bkg_and_sig_mc'):
    """
    Plots the M_inv distribution of the total bkg and the various background samples, for each (pt,eta) bin.
    These bkg distributions can be compared to the MC signal or the data distributions.
    """

    '''
    if plot_on_data and plot_on_signal:
        print("ERROR: plotting on data AND signal is not available now")
        sys.exit()
    '''

    file = ROOT.TFile(ws_filename, "READ")
    ws = file.Get("w")

    bins_pt, bins_pt = binning(binning_pt), binning(binning_eta)
    bin_dict = bin_dictionary(binning_pt, binning_eta)
    
    bins_ratio = [round(0.0 + 0.2*i, 2) for i in range(6)] + \
                    [round(1.0 + i, 2) for i in range(1, 10)] + \
                    [round(10.0 + 20*i, 2) for i in range(1, 10)] + [200.0]

    bins_sigma = [round(-10.0 + i, 2) for i in range(51)]
    bins_pull = [round(-2 + 0.08*i, 2) for i in range(51)]

    histos = {}

    if plot_on_data and plot_fit_bkgpdf:
        histos.update(init_pass_fail_histos("h_nbkg_pull", "Nbkg MC vs fit [nsigma]", 
                                            array('d', bins_sigma), bins_pt, bins_pt))
    if plot_on_signal and compare_bkgfrac:
        histos.update(init_pass_fail_histos("h_bkgfrac_pull", "Bkg fraction mc vs fit pull", 
                                            array('d', bins_pull), bins_pt, bins_pt))
        histos.update(init_pass_fail_histos("h_bkgfrac_ratio", "Bkg fraction mc/fit ratio", 
                                            array('d', bins_ratio), bins_pt, bins_pt))
        
    # if study_SS_bkg: bkg_categories = [cat+"_SS" if cat != "bkg_SameCharge" else cat for cat in bkg_categories]


    for bin_key, [gl_idx, bin_pt, bin_eta] in bin_dict.items():

        for flag in ["pass", "fail"]:

            # if bin_key != "[24.0to35.0][-0.1to0.0]": continue

            if logscale=="hybrid": setlog=True if flag=="pass" else False
            
            if plot_on_data:
                
                import_categories = ["data"] + bkg_categories

                if study_SS_bkg: import_categories += ["mc_SS"]

                datasets = make_bkg_dictionary(ws, import_categories, flag, bin_key, bin_pt, bin_eta)

                
                plot_bkg(datasets, flag, bin_key, logscale=setlog, figpath=f"{figpath}/minv_plots_w_data")

                '''
                ## NOT USED STUFF
                if plot_fit_bkgpdf and 1==0:
                    integral_histo_bkg = datasets["total_bkg"].sumEntries()
                    norm_bkg = datasets["pdf_bkg_fit"]["norm"]
                    err_integral_histo_bkg = sumw2_error(datasets["bkg_total"])

                    nbkg_pull = (integral_histo_bkg - norm_bkg.getVal())/norm_bkg.getError()
                    histos[f"h_nbkg_pull_{flag}"].Fill(nbkg_pull)
                    histos[f"h_nbkg_pull_{flag}_2d"].SetBinContent(bin_pt, bin_eta, nbkg_pull)
                '''


            if plot_on_signal:

                import_categories = ["mc"] + bkg_categories

                if study_SS_bkg: import_categories += ["mc_SS"]
        
                datasets = make_bkg_dictionary(ws, import_categories, flag, bin_key, bin_pt, bin_eta)

                plot_bkg(datasets, flag, bin_key, logscale=setlog, figpath=f"{figpath}/minv_plots_w_sig")
                
                '''
                ## NOT USED STUFF
                nbkg = datasets["bkg_total"].sumEntries()
                err_nbkg = sumw2_error(datasets["bkg_total"])

                nsignal = datasets["mc"].sumEntries()
                err_nsignal = sumw2_error(datasets["mc"])

                bkgfrac_mc, err_bkgfrac_mc = eval_bkgfrac(nbkg, nsignal, err_nbkg, err_nsignal)
            
                if compare_bkgfrac:
                    if type(datasets["fit_pars"]["nsig"]) is ROOT.RooRealVar:
                        fit_par_nsig = datasets["fit_pars"]["nsig"]
                    else:
                        fit_par_nsig = ROOT.RooRealVar(fit_par_nsig.GetName(), fit_par_nsig.GetTitle(), 1)
                    
                    if type(datasets["fit_pars"]["nbkg"]) is ROOT.RooRealVar:
                        fit_par_nbkg = datasets["fit_pars"]["nbkg"]
                    else:
                        fit_par_nbkg = ROOT.RooRealVar(fit_par_nbkg.GetName(), fit_par_nbkg.GetTitle(), 0)

                    bkgfrac_fit, err_bkgfrac_fit = eval_bkgfrac(fit_par_nbkg.getVal(), 
                                                                fit_par_nsig.getVal(),
                                                                fit_par_nbkg.getError(), 
                                                                fit_par_nsig.getError())
                
                    bkg_frac_pull = (bkgfrac_mc - bkgfrac_fit)/((err_bkgfrac_mc**2 + (err_bkgfrac_fit**2))**0.5)

                    if bkgfrac_fit == 0: 
                        bkgfrac_ratio = 199.55
                        #err_bkgfrac_ratio = bkgfrac_mc
                    else:
                        bkgfrac_ratio = bkgfrac_mc/bkgfrac_fit
                        # err_bkgfrac_ratio = bkgfrac_ratio*((err_bkgfrac_mc/bkgfrac_mc)**2 + (err_bkgfrac_fit/bkgfrac_fit)**2)**0.5

                    histos[f"h_bkgfrac_ratio_{flag}"].Fill(bkgfrac_ratio)
                    histos[f"h_bkgfrac_ratio_{flag}_2d"].SetBinContent(bin_pt, bin_eta, bkgfrac_ratio)
                    histos[f"h_bkgfrac_pull_{flag}"].Fill(bkg_frac_pull)
                    histos[f"h_bkgfrac_pull_{flag}_2d"].SetBinContent(bin_pt, bin_eta, bkg_frac_pull)
                '''

            if plot_on_data is False and plot_on_signal is False:
                datasets = make_bkg_dictionary(ws, bkg_categories, flag, bin_key, bin_pt, bin_eta)
                plot_bkg(datasets, flag, bin_key, logscale=setlog, figpath=f"{figpath}/minv_plots")
            
            
    saveHists = False
    if plot_on_data and plot_fit_bkgpdf:
        filename = f"{figpath}/nbkg_mc_vs_datafit.root"
        saveHists = True
    elif plot_on_signal and compare_bkgfrac:
        filename = f"{figpath}/bkgfrac_mc_vs_datafit.root"
        saveHists = True

    if saveHists:
        file = ROOT.TFile(filename, "RECREATE")
        file.cd()
        [histo.Write() for histo in histos.values()]
        file.Close()



def bkg_2d_distrib(ws_filename, bkg_categories, binning_pt, binning_eta,
                   study_SS_bkg=False, norm_data=False, norm_sig=False, norm_tot_bkg=False, 
                   plot_projected=False, filepath="bkg_studies"):
    """
    Makes 2d histograms containing the differential distribution of the various
    background samples; these distributions can be normalized w.r.t. data,
    the total amount of bkg events or the MC signal.
    """

    file = ROOT.TFile(ws_filename, "READ")
    ws = file.Get("w")

    if norm_data and norm_tot_bkg:
        sys.exit("ERROR: can't normalize both on data and total bkg")
    elif norm_data:
        add_title = " norm on data"
    elif norm_tot_bkg:
        add_title = " norm on total bkg"
    else:
        add_title = ""

    bins_pt, bins_eta = binning(binning_pt), binning(binning_eta)
    bin_dict = bin_dictionary(binning_pt, binning_eta)    

    plot_categories = [cat.replace("bkg_", "") for cat in bkg_categories]

    # if study_SS_bkg: 
        # plot_categories = [cat+"_SS" if cat != "SameCharge" else cat for cat in plot_categories]
        #plot_categories += ["mc_SS"]
        
    
    histos = {}

    isNorm = norm_data or norm_sig or norm_tot_bkg


    for bkg_cat in plot_categories:


        histos.update(init_pass_fail_histos(bkg_cat, bkg_cat, bins_pt, bins_pt, bins_eta))

        del histos[f"{bkg_cat}_pass"]
        del histos[f"{bkg_cat}_fail"]

        for bin_key, [gl_idx, bin_pt, bin_eta] in bin_dict.items():

            # if bin_key != "[24.0to26.0][-2.4to-2.3]": sys.exit()

            num_events = {}

            if bkg_cat != "mc_SS":
                type_dset, subs = "bkg", bkg_cat
            else:
                type_dset, subs = "mc", "SS"
        

            for flag in ["pass", "fail"]:
                bkg_hist = ws.data(f"Minv_{type_dset}_{flag}_{bin_key}_{subs}")
                num_events.update({
                    flag : bkg_hist.sumEntries(),
                    f"error {flag}" : sumw2_error(bkg_hist)
                    })
                
                if norm_data: 
                    h_norm = ws.data(f"Minv_data_{flag}_{bin_key}")
                elif norm_sig:
                    h_norm = ws.data(f"Minv_mc_{flag}_{bin_key}")
                elif norm_tot_bkg:
                    if type(ws.data(f"Minv_bkg_{flag}_{bin_key}_total")) is ROOT.RooDataHist:
                        h_norm = ws.data(f"Minv_bkg_{flag}_{bin_key}_total")
                    else:
                        axis = ws.var(f"x_{flag}_{bin_key}")
                        if binning_pt not in ["pt", "pt_tracking"] or binning_eta!="eta":
                            extended_bin_dict = bin_global_idx_dict(binning_pt, binning_eta)
                            bin_pt_ext, bin_eta_ext = extended_bin_dict[str(gl_idx[0])]
                        else:
                            bin_pt_ext, bin_eta_ext = bin_pt, bin_eta
                        h_norm = get_totbkg_roohist([ws, bkg_categories], flag, axis, bin_key, bin_pt_ext, bin_eta_ext)
                else:
                    pass

                if isNorm:
                    if h_norm.sumEntries() > 0 and num_events[flag] > 0:
                        if h_norm.sumEntries() < num_events[flag]:
                            print(num_events[flag], h_norm.sumEntries())
                        num_events[flag] = num_events[flag]/h_norm.sumEntries()
                        num_events[f"error {flag}"] = num_events[flag]*(
                            (num_events[f"error {flag}"]/num_events[flag])**2 + (sumw2_error(h_norm)/h_norm.sumEntries())**2)**0.5
                    else:
                        num_events[flag], num_events[f"error {flag}"] = 0, 0
                else:
                    pass
            
                histos[f"{bkg_cat}_{flag}_2d"].SetBinContent(bin_pt, bin_eta, num_events[flag])
                histos[f"{bkg_cat}_{flag}_2d"].SetBinError(bin_pt, bin_eta, num_events[f"error {flag}"])


        hist_dict_bkgcat = {"pass" : histos[f"{bkg_cat}_pass_2d"], "fail" : histos[f"{bkg_cat}_fail_2d"]}
        plot_2d_bkg_distrib(hist_dict_bkgcat, bkg_cat, figpath=f"{filepath}/bkg_2d_distrib")


    rootfile_distrib = ROOT.TFile(f"{filepath}/bkg_2d_distrib.root", "RECREATE")
    rootfile_distrib.cd()
    [histo.Write() for histo in histos.values()]
    rootfile_distrib.Close()

    if plot_projected: 
        hist_pass_dict, hist_fail_dict = {}, {}
        for bkg_cat in bkg_categories: 
            hist_pass_dict.update({f"{bkg_cat}_bkg" : histos[f"{bkg_cat}_pass_2d"]})
            hist_fail_dict.update({f"{bkg_cat}_bkg" : histos[f"{bkg_cat}_fail_2d"]})

        plot_projected_bkg(hist_pass_dict, binning_pt, "eta_singlebin", "pass", figpath=filepath)
        plot_projected_bkg(hist_fail_dict, binning_pt, "eta_singlebin", "fail", figpath=filepath)
        plot_projected_bkg(hist_pass_dict, "pt_singlebin", binning_eta, "pass", figpath=filepath)
        plot_projected_bkg(hist_fail_dict, "pt_singlebin", binning_eta, "fail", figpath=filepath)


