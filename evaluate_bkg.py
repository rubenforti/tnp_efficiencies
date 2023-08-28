"""
"""

import ROOT
import sys
from utilities.base_library import lumi_factors, binning, bin_dictionary, sumw2_error
from utilities.base_library import eval_efficiency as eval_bkgfrac  #La formula è la stessa
from utilities.results_utils import init_pass_fail_histos
from utilities.plot_utils import plot_bkg_on_histo
from utilities.fit_utils import check_existing_fit
from utilities.dataset_utils import show_negweighted_bins
from array import array
# from indep_eff_pseudodata import independent_efficiency


# bkg_categories = ["WW", "WZ", "ZZ", "TTSemileptonic", "Ztautau"]

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
            "histo_pdf" : ROOT.RooHistPdf(f"{cat}_bkg_pdf", f"{cat}_bkg_pdf", ROOT.RooArgSet(axis), bkg_histo),
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
            "histo_pdf" : ROOT.RooHistPdf(f"MC_signal_pdf", f"MC_signal_pdf", ROOT.RooArgSet(axis), histo_mc),
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

def fit_bkg(type_eff, type_analysis, ws, flag, bin_key, fit_shape="expo"):
    """
    """

    bkg_dict = make_bkg_dictionary(type_eff, ws, flag, bin_key, bkg_categories, import_data=True)


    axis = bkg_dict["axis"]

    histo_bkg = bkg_dict["total_bkg"]["roohisto"]

    
    if fit_shape == "expo":
        tau = ROOT.RooRealVar(f"tau_bkg_{flag}_{bin_key}", "tau", 0.0, -5., 5.)
        pdf_bkg = ROOT.RooExponential(f"pdf_bkg_{flag}_{bin_key}", "pdf_bkg", axis, tau)
    elif fit_shape == "cmsshape":
        pass

    nbkg = ROOT.RooRealVar(f"nbkg_{flag}_{bin_key}", "nbkg", histo_bkg.sumEntries(), 0.5, 1.5*histo_bkg.sumEntries())

    model = ROOT.RooExtendPdf(f"model_bkg_{flag}_{bin_key}", "model_bkg", pdf_bkg, nbkg)

    res_bkg = model.fitTo(histo_bkg,
                          ROOT.RooFit.Save(True),
                          ROOT.RooFit.Extended(True),
                          ROOT.RooFit.Strategy(2),
                          ROOT.RooFit.Minimizer("Minuit2", "migrad"),
                          ROOT.RooFit.SumW2Error(True), 
                          ROOT.RooFit.PrintLevel(-1))

    res_bkg.SetName(f"results_{flag}_{bin_key}")
    
    ws.Import(histo_bkg)
    ws.Import(model)
    ws.Import(res_bkg)



###############################################################################

def bkg_mass_distribution(type_eff, bkg_categories, ws, binning_pt, binning_eta, 
                           plot_on_data=False, plot_on_signal=False, logscale=True,
                           figpath='figs/bkg_and_sig_mc'):
    """
    Plots the M_inv distribution of the total bkg and the various background samples, for each (pt,eta) bin.
    These bkg distributions can be compared to the MC signal or the data distributions.
    """

    if plot_on_data and plot_on_signal:
        print("ERROR: plotting on data AND signal is not available now")
        sys.exit()


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
    if plot_on_signal:
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
                                               import_data=True, import_fit_pdf_bkg=True)
                
                integral_histo_bkg = datasets["total_bkg"]["integral"]
                
                norm_bkg = datasets["pdf_bkg_fit"]["norm"]
                err_integral_histo_bkg = sumw2_error(datasets["total_bkg"]["roohisto"])

                nbkg_pull = (integral_histo_bkg - norm_bkg.getVal())/norm_bkg.getError()
                histos[f"h_nbkg_pull_{flag}"].Fill(nbkg_pull)
                histos[f"h_nbkg_pull_{flag}_2d"].SetBinContent(bin_pt, bin_eta, nbkg_pull)
                
                plot_bkg_on_histo(datasets, flag, bin_key, figpath='figs/bkg_mc_on_data') 


            if plot_on_signal:
        
                datasets = make_bkg_dictionary(type_eff, ws, flag, bin_key, bkg_categories, 
                                               import_mc_signal=True, import_fit_pars=True)

                nbkg = datasets["total_bkg"]["integral"]
                err_nbkg = sumw2_error(datasets["total_bkg"]["roohisto"])

                nsignal = datasets["MC_signal"]["integral"]
                err_nsignal = sumw2_error(datasets["MC_signal"]["roohisto"])

                print(nsignal, err_nsignal)
                
                if type(datasets["fit_pars"]["nsig"]) is ROOT.RooRealVar:
                    fit_par_nsig = datasets["fit_pars"]["nsig"]
                else:
                    fit_par_nsig = ROOT.RooRealVar(fit_par_nsig.GetName(), fit_par_nsig.GetTitle(), 1)
                
                if type(datasets["fit_pars"]["nbkg"]) is ROOT.RooRealVar:
                    fit_par_nbkg = datasets["fit_pars"]["nbkg"]
                else:
                    fit_par_nbkg = ROOT.RooRealVar(fit_par_nbkg.GetName(), fit_par_nbkg.GetTitle(), 0)

                bkgfrac_mc, err_bkgfrac_mc = eval_bkgfrac(nbkg, nsignal, err_nbkg, err_nsignal)
                
                bkgfrac_fit, err_bkgfrac_fit = eval_bkgfrac(fit_par_nbkg.getVal(), fit_par_nsig.getVal(),
                                                                fit_par_nbkg.getError(), fit_par_nsig.getError())
                
                bkg_frac_pull = (bkgfrac_mc - bkgfrac_fit)/((err_bkgfrac_mc**2 + (err_bkgfrac_fit**2))**0.5)

                if bkgfrac_fit == 0:
                    bkgfrac_ratio = 199.55
                    err_bkgfrac_ratio = bkgfrac_mc
                else:
                    bkgfrac_ratio = bkgfrac_mc/bkgfrac_fit
                    err_bkgfrac_ratio = bkgfrac_ratio*((err_bkgfrac_mc/bkgfrac_mc)**2 + (err_bkgfrac_fit/bkgfrac_fit)**2)**0.5

                histos[f"h_bkgfrac_ratio_{flag}"].Fill(bkgfrac_ratio)
                histos[f"h_bkgfrac_ratio_{flag}_2d"].SetBinContent(bin_pt, bin_eta, bkgfrac_ratio)
                histos[f"h_bkgfrac_pull_{flag}"].Fill(bkg_frac_pull)
                histos[f"h_bkgfrac_pull_{flag}_2d"].SetBinContent(bin_pt, bin_eta, bkg_frac_pull)

                plot_bkg_on_histo(datasets, flag, bin_key, figpath=figpath)

            else:
                datasets = make_bkg_dictionary(type_eff, ws, flag, bin_key, bkg_categories)
                plot_bkg_on_histo(datasets, flag, bin_key, logscale=logscale, figpath=figpath)
                



    '''
    if plot_on_data:
        filename = f"bkg_results/nbkg_mc_vs_datafit.root"
    elif plot_on_signal:
        filename = f"bkg_results/bkgfrac_mc_vs_datafit.root"
        #filename = "prova.root"

    file = ROOT.TFile(filename, "RECREATE")
    file.cd()
    for histo in histos.values():
        histo.Write()
    file.Close()
    '''


###############################################################################

def fit_on_pseudodata(ws, binning_pt, binning_eta, bkg_shape, bkg_categories,
                      refit_numbkg=False, test_bkg=False, verb=-1, figs=False):
    """

    """

    bin_dict = bin_dictionary(binning_pt, binning_eta)

    prob_bins = []

    bins_pt = binning(binning_pt)
    bins_eta = binning(binning_eta)

    nbins_pt, nbins_eta = len(bins_pt)-1, len(bins_eta)-1

    
    Nproblems = 0

    # key_prova = "[24.0to26.0][0.6to1.2]"
    # bin_dict = {key_prova : bin_dict[key_prova]}


    for bin_key in bin_dict:

        _, bin_pt, bin_eta = bin_dict[bin_key]

        existingRes = check_existing_fit("indep", ws, bin_key)

        print(existingRes)

        if existingRes == 0:
            res_pass, res_fail, status = independent_efficiency(
                ws, bin_key, bkg_shape, bkg_categories,
                refit_numbkg=refit_numbkg, test_bkg=test_bkg, verb=verb, figs=figs)
        else:
            res_pass, res_fail = existingRes
            status = True
            

        if status is False:
            print(f"\nBin {bin_key} ({bin_pt}|{bin_eta}) has problems!\n")
            Nproblems += 1
            prob_bins.append(f"{bin_key}")
            print("****")
            res_pass.Print()
            res_pass.correlationMatrix().Print()
            '''
            pars_pass = res_pass.floatParsFinal()
            nsig_pass = pars_pass.find(f"nsig_pass_{bin_key}")
            print((nsig_pass.getVal()**0.5, nsig_pass.getError()))
            '''
            print("****")
            res_fail.Print()
            res_fail.correlationMatrix().Print()
            '''

            pars_fail = res_fail.floatParsFinal()
            nsig_fail = pars_fail.find(f"nsig_fail_{bin_key}")
            print((nsig_fail.getVal()**0.5, nsig_fail.getError()))
            '''
            print("****")
            print(res_pass.status(), res_fail.status())
            print(res_pass.covQual(), res_fail.covQual())
            print(res_pass.edm(), res_fail.edm())
            print('\n')
        else:
            pass
            #res_object.add_result(ws, bin_pt, bin_eta)


    print(len(prob_bins))
    print(prob_bins)

###############################################################################

def plot_backgrounds_2d(ws, binning_pt, binning_eta, bkg_categories, 
                        norm_data=False, norm_tot_bkg=False, filename="root_files/backgrounds/bkg_2d_distributions.root"):
    """
    Makes 2d histograms containing the differential distribution of the various background samples; these
    distributions can be normalized w.r.t. data or the total amount of bkg events.
    """

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

    rootfile_distrib = ROOT.TFile(filename, "RECREATE")
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
                d_npass, d_nfail = d_npass/h_data_pass.sumEntries(), d_nfail/h_data_fail.sumEntries()

            if norm_tot_bkg:
                h_totbkg_pass = ws.data(f"Minv_bkg_pass_{bin_key}_total")
                h_totbkg_fail = ws.data(f"Minv_bkg_fail_{bin_key}_total")
                # print(n_pass, h_data_pass.sumEntries())
                npass, nfail = npass/h_data_pass.sumEntries(), nfail/h_data_fail.sumEntries()
                d_npass, d_nfail = d_npass/h_data_pass.sumEntries(), d_nfail/h_data_fail.sumEntries()
            
        
            histos[f"{cat}_pass_2d"].SetBinContent(bin_pt, bin_eta, npass)
            histos[f"{cat}_fail_2d"].SetBinContent(bin_pt, bin_eta, nfail)

        '''
        c = ROOT.TCanvas("c", "c", 1200, 1200)
        c.Divide(1,2)

        c.cd(1)
        ROOT.gStyle.SetOptStat("")
        histos_pass[cat].Draw("COLZ")
        c.cd(2)
        ROOT.gStyle.SetOptStat("")
        histos_fail[cat].Draw("COLZ")

        if norm_data:
            c.SaveAs(f"figs/backgrounds/{cat}_bkg_h2d_norm_data.pdf")
        elif norm_tot_bkg:
            c.SaveAs(f"figs/backgrounds/{cat}_bkg_h2d_norm_totbkg.pdf") 
        else:
            c.SaveAs(f"figs/backgrounds/{cat}_bkg_h2d.pdf")
        '''

        rootfile_distrib.cd()
        histos[f"{cat}_pass_2d"].Write()
        histos[f"{cat}_fail_2d"].Write()
    
    rootfile_distrib.Close()

       
###############################################################################
###############################################################################


if __name__ == "__main__":

    ROOT.gROOT.SetBatch(True)
    ROOT.PyConfig.IgnoreCommandLineOptions = True

    bkg_categories = ["WW", "WZ", "ZZ", "TTSemileptonic", "Ztautau"]
    # bkg_categories = ["Ztautau"]

    type_eff = "triggerminus"

    type_analysis = "indep"

    lumi_scales = lumi_factors(type_eff, bkg_categories)

    sig_lumi_scale = lumi_scales.pop("Zmumu")
    # lumi_scales = {"WW":1, "WZ":1, "ZZ":1, "TTSemileptonic":1, "Ztautau":1}

    '''
    filename_data = "/scratchnvme/wmass/Steve_root_files/Standard_SF_files/tnp_iso_data_vertexWeights1_oscharge1.root"
    filename_mc = "/scratchnvme/wmass/Steve_root_files/Standard_SF_files/tnp_iso_mc_vertexWeights1_oscharge1.root"
    filename_mc_weightsum = "/scratchnvme/rajarshi/Signal_TNP_3D_Histograms/OS/tnp_iso_mc_vertexWeights1_oscharge1.root"
    dirname_bkg = "/scratchnvme/rajarshi/Bkg_TNP_3D_Histograms/OS"
    '''

    filename_data = "root_files/datasets/tnp_triggerminus_data_vertexWeights1_oscharge1.root"
    filename_mc = "root_files/datasets/tnp_triggerminus_mc_vertexWeights1_oscharge1.root"
    dirname_bkg = "root_files/datasets"

    bkg_filenames = {}
    [bkg_filenames.update({cat : 
        f"{dirname_bkg}/tnp_{type_eff}_{cat}_vertexWeights1_oscharge1.root"}) for cat in bkg_categories]


    import_dictionary = {
        #"data" : filename_data,
        "mc" : {
            "filename" : filename_mc,
            "lumi_scale" : sig_lumi_scale
        },
        "bkg" : {
            "filenames" : bkg_filenames,
            "lumi_scales" : lumi_scales
        }
    }
    
    binning_pt_key = "pt_12bins"
    binning_eta_key = "eta_16bins"

    binning_mass = "mass_60_120"


    # bin_set = bin_dict()
    bin_set = bin_dictionary(binning_pt_key, binning_eta_key)

    ws_filename = "root_files/ws_triggerminus_pseudodata.root"

    # ws = ws_init(import_dictionary, an, binning_pt_key, binning_eta_key, binning_mass)
    # ws.writeToFile(ws_filename)

    # file = ROOT.TFile("root_files/ws_bkg_studies.root", "READ")
    # file = ROOT.TFile("root_files/ws_bkg_pseudodata.root", "READ")

    file = ROOT.TFile(ws_filename, "READ")
    ws = file.Get("w")

    binning_mass = binning("mass_60_120")

    '''
    binnings_list = [["pt", "eta"], ["pt", "eta_24bins"], ["pt", "eta_16bins"], 
                     ["pt", "eta_8bins"], ["pt_12bins", "eta"], ["pt_9bins", "eta"],
                     ["pt_6bins", "eta"], ["pt_12bins", "eta_24bins"], ["pt_9bins", "eta_16bins"],
                     ["pt_12bins", "eta_16bins"], ["pt_9bins", "eta_24bins"]]
    
    for binning_pt, binning_eta in binnings_list:
        bin_set = bin_dictionary(binning_pt, binning_eta)
        w = ws_init(import_dictionary, an, bin_set, binning_mass)
        show_negweighted_bins(w, bkg_categories, binning_pt, binning_eta)
    '''
    show_negweighted_bins(type_eff, ws, bkg_categories, "pt", "eta")


    '''
    plot_backgrounds_2d(ws, binning_pt_key, binning_eta_key, bkg_categories,
                        norm_data=False, norm_tot_bkg=False, 
                        filename="bkg_2d_distr_pt12_eta16.root")
    '''
                        
    '''
    bkg_mass_distribution("iso", bkg_categories, ws, binning_pt_key, binning_eta_key,
                          plot_on_data=False, plot_on_signal=False, logscale=False,
                          figpath="figs/bkg_pt12_eta16")
    '''


    # fit_on_pseudodata(ws, "pt", "eta_8bins", "expo", bkg_categories, refit_numbkg=True, verb=-1, figs=True)
    # ws.writeToFile(ws_filename)

    # ws.Print()

    # plot_pseudodata_eff_comparison(ws, "pt", "eta_8bins", "bkg_results/pseudodata_eff_comparison.root")


    #ws.writeToFile("root_files/ws/ws_data_bkg.root")

    '''
    for bin_key in bin_set:
        fit_bkg("iso", an, ws, "pass", bin_key)

    ws.writeToFile("root_files/ws/ws_bkg_fit.root")
    '''



    '''
    h_pass_gen = ROOT.TH2D(h_pass[0])
    h_pass_gen.Sumw2()
    print(h_pass[0].GetBinContent(1,1), h_pass[0].GetBinError(1,1))

    h_fail_gen = ROOT.TH2D(h_fail[0])
    h_fail_gen.Sumw2()

    for i in range(1, 5):
        print(h_pass[i].GetBinContent(1,1), h_pass[i].GetBinError(1,1))
        h_pass_gen.Add(h_pass[i])
        h_fail_gen.Add(h_fail[i])
    

    c = ROOT.TCanvas()
    c.Divide(1,2)
    c.cd(1)
    ROOT.gStyle.SetOptStat("en")
    h_pass_gen.Draw("COLZ")
    c.cd(2)
    ROOT.gStyle.SetOptStat("en")
    h_fail_gen.Draw("COLZ")
    '''

    '''
    res = results_manager("indep")
    res.Open("results/benchmark_iso/results_iso_indep_benchmark.pkl")
    res_dict = res.dictionary()

    print("*********")
    # print(h_pass_gen.GetBinContent(1,1), h_pass_gen.GetBinError(1,1))
    print(res_dict["1,1"]["pars_pass"].find("nbkg_pass_(1|1)"))
    print(res_dict["1,1"]["pars_fail"].find("nbkg_fail_(1|1)"))
    '''













