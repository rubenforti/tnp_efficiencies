"""
"""
import ROOT
import sys
from utilities.base_library import lumi_factors, binning, bin_dictionary, sumw2_error
from utilities.base_library import eval_efficiency as eval_bkgfrac  #La formula è la stessa
from utilities.results_utils import init_pass_fail_histos
from utilities.plot_utils import plot_bkg
from array import array


bkg_categories = ["WW", "WZ", "ZZ", "TTSemileptonic", "Ztautau"]

colors = { 
    "WW" : ROOT.kOrange+1,
    "WZ" : ROOT.kYellow+3,
    "ZZ" : ROOT.kGreen+1,
    "TTSemileptonic" : ROOT.kCyan+1,
    "Ztautau" : ROOT.kMagenta+1,
    "SameCharge" : ROOT.kOrange+10,
    "total_bkg" : ROOT.kRed,
    "pdf_bkg_fit" : ROOT.kRed,
    "signal" : ROOT.kBlue
}

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
    
    cnt_mergedpt=0
    for bin_key in bin_dict.keys():

        print(bin_key)

        _, bin_pt, bin_eta = bin_dict[bin_key]

        # Bin transformation needed in case the bins are merged
        if type(bin_eta) is list:
            bin_eta = int(1+(nbins_eta*(bin_eta[0]-1)/48.))
        if type(bin_pt) is list:
            bin_pt_list = bin_pt
            bin_pt = int(bin_pt_list[0] - cnt_mergedpt)
            cnt_mergedpt += bin_pt_list[-1]-bin_pt_list[0] if bin_eta==nbins_eta else 0

        cnt_sum_total = 0
        cnt_sum_single = 0
        cnt_weight_total = 0
        cnt_weight_single = 0

        for flag in ["pass", "fail"]:
            axis = ws.var(f"x_{flag}_{bin_key}")
            tot_histo = ROOT.RooDataHist("tot_histo", "tot_histo", ROOT.RooArgSet(axis), "")

            for cat in bkg_categories:
                bkg_data = ws.data(f"Minv_bkg_{flag}_{bin_key}_{cat}")
                tot_histo.add(bkg_data) 
                cnt_sum_single = cnt_sum_single+1 if bkg_data.sumEntries() < 0 else cnt_sum_single
                    
            cnt_sum_total = cnt_sum_total+1 if tot_histo.sumEntries() < 0 else cnt_sum_total

            for idx in range(tot_histo.numEntries()):
                cnt_weight_total = cnt_weight_total+1 if tot_histo.weight(idx) < 0 else cnt_weight_total
                for cat in bkg_categories:
                    histo = ws.data(f"Minv_bkg_{flag}_{bin_key}_{cat}")
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

def make_bkg_dictionary(type_eff, ws, flag, bin_key, bkg_categories, 
                        import_data=False, import_mc_signal=False, 
                        import_fit_pars=False, import_fit_pdf_bkg=False, pdf_bkg_shape="expo"):
    """
    Creates a dictionary containing the bkg MC datasets and the useful information related to them.
    """
    
    lumi_scales = lumi_factors(type_eff, bkg_categories)

    axis = ws.var(f"x_{flag}_{bin_key}")
    axis.setBins(60, "plot_binning")

    datasets = {}
    datasets.update({"axis" : axis})

    for cat in bkg_categories:

        bkg_histo = ws.data(f"Minv_bkg_{flag}_{bin_key}_{cat}")

        datasets.update({f"{cat}_bkg" : {
            "roohisto" : bkg_histo,
            "histo_pdf" : ROOT.RooHistPdf(f"{cat}_bkg_pdf", f"{cat}_bkg_pdf", ROOT.RooArgSet(axis), bkg_histo, 0),
            "lumi_scale" : lumi_scales[cat],
            "integral" : bkg_histo.sumEntries(),
            "color" : colors[cat]
            }
        })

    bkg_total_histo = ROOT.RooDataHist(f"Minv_bkg_{flag}_{bin_key}_total", "bkg_total_histo", 
                                       ROOT.RooArgSet(axis), "plot_binning")

    for bkg_cat in bkg_categories:
        bkg_total_histo.add(datasets[f"{bkg_cat}_bkg"]["roohisto"])

    datasets.update({"total_bkg" : {
        "roohisto": bkg_total_histo,
        "lumi_scale" : 1,
        "integral" : bkg_total_histo.sumEntries(),
        "color" : colors["total_bkg"]
        }
    })

    if import_data:
        histo_data = ws.data(f"Minv_data_{flag}_{bin_key}")
        datasets.update({"data" : histo_data})
    
    if import_mc_signal:
        histo_mc = ws.data(f"Minv_mc_{flag}_{bin_key}")
        datasets.update({"MC_signal" : {
            "roohisto" : histo_mc,
            "histo_pdf" : ROOT.RooHistPdf(f"MC_signal_pdf", f"MC_signal_pdf", ROOT.RooArgSet(axis), histo_mc, 0),
            "lumi_scale" : lumi_scales["Zmumu"],
            "integral" : histo_mc.sumEntries(),
            "color" : colors["signal"]
            }
        })
    
    if import_fit_pars:
        res_obj = ws.obj(f"results_{flag}_{bin_key}")
        pars = res_obj.floatParsFinal()
        datasets.update({"fit_pars" : {
            "nsig" : pars.find(f"nsig_{flag}_{bin_key}"),
            "nbkg" : pars.find(f"nbkg_{flag}_{bin_key}"),
            }
        })

    if import_fit_pdf_bkg:
        pdf_bkg_fit = ws.pdf(f"{pdf_bkg_shape}_bkg_{flag}_{bin_key}")
        if type(pdf_bkg_fit) is ROOT.TObject:
            print(f"ERROR: bkg pdf not found at bin {bin_key}")
            sys.exit()

        res_obj = ws.obj(f"results_{flag}_{bin_key}")
        pars = res_obj.floatParsFinal()
        norm_bkg = pars.find(f"nbkg_{flag}_{bin_key}")

        if type(norm_bkg) is not ROOT.RooRealVar:
            norm_bkg = ROOT.RooRealVar(f"nbkg_{flag}_{bin_key}", f"nbkg_{flag}_{bin_key}", 0)

        datasets.update({"pdf_bkg_fit" : {
            "pdf" : pdf_bkg_fit,
            "norm" : norm_bkg,
            "color" : colors["pdf_bkg_fit"]
            }
        })

    return datasets

###############################################################################

def bkg_mass_distribution(type_eff, ws_filename, bkg_categories, binning_pt, binning_eta, 
                           plot_on_data=False, plot_fit_bkgpdf=False,
                           plot_on_signal=False, compare_bkgfrac=False,
                           logscale=True, figpath='figs/bkg_and_sig_mc'):
    """
    Plots the M_inv distribution of the total bkg and the various background samples, for each (pt,eta) bin.
    These bkg distributions can be compared to the MC signal or the data distributions.
    """

    if plot_on_data and plot_on_signal:
        print("ERROR: plotting on data AND signal is not available now")
        sys.exit()

    file = ROOT.TFile(ws_filename, "READ")
    ws = file.Get("w")

    bins_pt, bins_pt = binning(binning_pt), binning(binning_eta)
    nbins_pt, nbins_eta = len(bins_pt)-1, len(bins_pt)-1
    bin_dict = bin_dictionary(binning_pt, binning_eta)
    

    bins_ratio = [round(0.0 + 0.2*i, 2) for i in range(6)] + \
                    [round(1.0 + i, 2) for i in range(1, 10)] + \
                    [round(10.0 + 20*i, 2) for i in range(1, 10)] + [200.0]

    bins_sigma = [round(-10.0 + i, 2) for i in range(51)]
    bins_pull = [round(-2 + 0.08*i, 2) for i in range(51)]

    histos = {}

    if plot_on_data:
        histos.update(init_pass_fail_histos("h_nbkg_pull", "Nbkg MC vs fit [nsigma]", 
                                            array('d', bins_sigma), bins_pt, bins_pt))
    if plot_on_signal and compare_bkgfrac:
        histos.update(init_pass_fail_histos("h_bkgfrac_pull", "Bkg fraction mc vs fit pull", 
                                            array('d', bins_pull), bins_pt, bins_pt))
        histos.update(init_pass_fail_histos("h_bkgfrac_ratio", "Bkg fraction mc/fit ratio", 
                                            array('d', bins_ratio), bins_pt, bins_pt))

    for bin_key in bin_dict:

        _, bin_pt, bin_eta = bin_dict[bin_key]

        # Bin transformation needed in case the bins are merged
        if type(bin_pt) is list:
            bin_pt = int(1+(nbins_pt*(bin_pt[0]-1)/15.))
        if type(bin_eta) is list:
            bin_eta = int(1+(nbins_eta*(bin_eta[0]-1)/48.))
    

        for flag in ["pass", "fail"]:
            
            if plot_on_data:

                datasets = make_bkg_dictionary(type_eff, ws, flag, bin_key, bkg_categories, 
                                               import_data=True, import_fit_pdf_bkg=plot_fit_bkgpdf)
                
                integral_histo_bkg = datasets["total_bkg"]["integral"]
                
                norm_bkg = datasets["pdf_bkg_fit"]["norm"]
                err_integral_histo_bkg = sumw2_error(datasets["total_bkg"]["roohisto"])

                nbkg_pull = (integral_histo_bkg - norm_bkg.getVal())/norm_bkg.getError()
                histos[f"h_nbkg_pull_{flag}"].Fill(nbkg_pull)
                histos[f"h_nbkg_pull_{flag}_2d"].SetBinContent(bin_pt, bin_eta, nbkg_pull)
                
                plot_bkg(datasets, flag, bin_key, figpath=f"{figpath}/minv_plots_w_data") 


            if plot_on_signal:
        
                datasets = make_bkg_dictionary(type_eff, ws, flag, bin_key, bkg_categories, 
                                               import_mc_signal=True, import_fit_pars=compare_bkgfrac)
                
                plot_bkg(datasets, flag, bin_key, logscale=logscale, figpath=f"{figpath}/minv_plots_w_sig")
                
                nbkg = datasets["total_bkg"]["integral"]
                err_nbkg = sumw2_error(datasets["total_bkg"]["roohisto"])

                nsignal = datasets["MC_signal"]["integral"]
                err_nsignal = sumw2_error(datasets["MC_signal"]["roohisto"])

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

            else:
                datasets = make_bkg_dictionary(type_eff, ws, flag, bin_key, bkg_categories)
            
            
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



def gen_bkg_2d_distrib(ws_filename, bkg_categories, binning_pt, binning_eta, 
                       norm_data=False, norm_sig=False, norm_tot_bkg=False, filepath="bkg_studies"):
    """
    Makes 2d histograms containing the differential distribution of the various
    background samples; these distributions can be normalized w.r.t. data,
    the total amount of bkg events or the MC signal.
    """

    file = ROOT.TFile(ws_filename, "READ")
    ws = file.Get("w")

    if norm_data and norm_tot_bkg:
        print("ERROR: can't normalize both on data and total bkg")
        sys.exit()
    elif norm_data:
        add_title = " norm on data"
    elif norm_tot_bkg:
        add_title = " norm on total bkg"
    else:
        add_title = ""

    bins_pt, bins_eta = binning(binning_pt), binning(binning_eta)
    bin_dict = bin_dictionary(binning_pt, binning_eta)    

    nbins_pt = len(bins_pt) - 1
    nbins_eta = len(bins_eta) - 1

    rootfile_distrib = ROOT.TFile(f"{filepath}/bkg_2d_distrib.root", "RECREATE")
    
    histos = {}

    for cat in bkg_categories:

        histos.update(init_pass_fail_histos(cat, cat, bins_pt, bins_pt, bins_eta))

        cnt_mergedpt = 0

        for bin_key in bin_dict:

            _, bin_pt, bin_eta = bin_dict[bin_key]

            # Bin transformation needed in case the bins are merged
            if type(bin_eta) is list:
                bin_eta = int(1+(nbins_eta*(bin_eta[0]-1)/48.))
            if type(bin_pt) is list:
                bin_pt_list = bin_pt
                bin_pt = int(bin_pt_list[0] - cnt_mergedpt)
                cnt_mergedpt += bin_pt_list[-1]-bin_pt_list[0] if bin_eta==nbins_eta else 0
                    
            h_pass = ws.data(f"Minv_bkg_pass_{bin_key}_{cat}")
            h_fail = ws.data(f"Minv_bkg_fail_{bin_key}_{cat}")
            npass = h_pass.sumEntries()
            d_npass = sumw2_error(h_pass)
            nfail = h_fail.sumEntries()
            d_nfail = sumw2_error(h_fail)


            if norm_data:
                h_data_pass = ws.data(f"Minv_data_pass_{bin_key}")
                h_data_fail = ws.data(f"Minv_data_fail_{bin_key}")
                # print(n_pass, h_data_pass.sumEntries())
                npass, nfail = npass/h_data_pass.sumEntries(), nfail/h_data_fail.sumEntries()   
                d_npass = npass*((d_npass/npass)**2 + (1/h_data_pass.sumEntries())**2)**0.5
                d_nfail = nfail*((d_nfail/nfail)**2 + (1/h_data_fail.sumEntries())**2)**0.5

            if norm_sig:
                h_mc_pass = ws.data(f"Minv_mc_pass_{bin_key}")
                h_mc_fail = ws.data(f"Minv_mc_fail_{bin_key}")
                # print(n_pass, h_data_pass.sumEntries())
                d_npass = npass*((d_npass/npass)**2 + 
                                 (sumw2_error(h_mc_pass)/h_mc_pass.sumEntries())**2)**0.5
                d_nfail = nfail*((d_nfail/nfail)**2 + 
                                 (sumw2_error(h_mc_fail)/h_mc_fail.sumEntries())**2)**0.5

            if norm_tot_bkg:
                h_totbkg_pass = ws.data(f"Minv_bkg_pass_{bin_key}_total")
                h_totbkg_fail = ws.data(f"Minv_bkg_fail_{bin_key}_total")
                # print(n_pass, h_data_pass.sumEntries())
                npass, nfail = npass/h_data_pass.sumEntries(), nfail/h_data_fail.sumEntries()
                d_npass = npass*((d_npass/npass)**2 + 
                                 (sumw2_error(h_totbkg_pass)/h_totbkg_pass.sumEntries())**2)**0.5
                d_nfail = nfail*((d_nfail/nfail)**2 + 
                                 (sumw2_error(h_totbkg_fail)/h_totbkg_fail.sumEntries())**2)**0.5
            
        
            histos[f"{cat}_pass_2d"].SetBinContent(bin_pt, bin_eta, npass)
            histos[f"{cat}_pass_2d"].SetBinError(bin_pt, bin_eta, d_npass)
            histos[f"{cat}_fail_2d"].SetBinContent(bin_pt, bin_eta, nfail)
            histos[f"{cat}_fail_2d"].SetBinError(bin_pt, bin_eta, d_nfail)

        rootfile_distrib.cd()
        histos[f"{cat}_pass_2d"].Write()
        histos[f"{cat}_fail_2d"].Write()
    
    rootfile_distrib.Close()
