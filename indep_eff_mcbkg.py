"""
"""
import ROOT
import time
import os
import sys
from multiprocessing import Pool
from utilities.results_utils import results_manager, efficiency_from_res
from utilities.plot_utils import plot_fitted_pass_fail
from utilities.dataset_utils import ws_init
from utilities.base_library import eval_efficiency, lumi_factors, bin_dictionary, binning, sumw2_error
from utilities.fit_utils import fit_quality, check_chi2, check_existing_fit


def independent_efficiency(ws, bin_key, bkg_categories, refit_numbkg=False, verb=-1, figs=False):
  
    path = os.path.dirname(__file__)
    ROOT.gSystem.cd(path)

    # ----------------------------------
    #  Initial parameters - MODIFY HERE
    # ----------------------------------
    NBINS = [2000, 2000]
    bufFractions = [0.5, 0.5]

    fit_strategy = [2, 2]

    # Lower limit for NSIG and upper limit for NBKG, in terms of total events in the histogram, for pass and fail respectively
    lims_num = [[0.5, 0.2], [0.5, 0.2]]

    mean_p = ROOT.RooRealVar(
        f"mean_pass_{bin_key}", "mean pass", 0, -5.0, 5.0)
    sigma_p = ROOT.RooRealVar(
        f"sigma_pass_{bin_key}", "sigma pass", 0.5, 0.1, 5.0)
    mean_f = ROOT.RooRealVar(
        f"mean_fail_{bin_key}", "mean fail", 0, -5.0, 5.0)
    sigma_f = ROOT.RooRealVar(
        f"sigma_fail_{bin_key}", "sigma fail",0.5, 0.1, 5.0)

    # ------------------------------
    #  Histograms and PDFs building
    # ------------------------------
    mean, sigma = [mean_f, mean_p], [sigma_f, sigma_p]
    axis, histo_data, pdf_mc, bkg_hist = [0, 0], [0, 0], [0, 0], [0, 0]
    smearing, conv_pdf, background = [0, 0], [0, 0], [0, 0]

    nexp, Nsig, Nbkg = [0, 0], [0, 0], [0, 0]
    sum_func, model, results, fit_status = [0, 0], [0, 0], [0, 0], [0, 0]

    for cond in ["pass", "fail"]:  # PASS has index 1, FAIL has index 0

        idx = 1 if cond == "pass" else 0

        axis[idx] = ws.var(f"x_{cond}_{bin_key}")
        axis[idx].setRange("fitRange", 60.0, 120.0)

        # binning = ROOT.RooUniformBinning(
        #     axis.getRange("fitRange")[0], axis.getRange("fitRange")[1], NBINS[idx])
        axis[idx].setBinning(
            ROOT.RooUniformBinning(axis[idx].getRange("fitRange")[0],
                                   axis[idx].getRange("fitRange")[1], NBINS[idx]), "cache")

        histo_data[idx] = ws.data(f"Minv_data_{cond}_{bin_key}")

        pdf_mc[idx] = ROOT.RooHistPdf(
            f"pdf_mc_{cond}_{bin_key}", "pdf MC", axis[idx],
            ws.data(f"Minv_mc_{cond}_{bin_key}"), 3)

        smearing[idx] = ROOT.RooGaussian(f"smearing_{cond}_{bin_key}",
                                         "Gaussian smearing", axis[idx], mean[idx], sigma[idx])

        conv_pdf[idx] = ROOT.RooFFTConvPdf(f"conv_{cond}_{bin_key}", f"Convolution pdf",
                                           axis[idx], pdf_mc[idx], smearing[idx])
        conv_pdf[idx].setBufferFraction(bufFractions[idx])
        # conv_pdf[idx].setBufferStrategy(0)
        # conv_pdf[idx].setOperMode(3)

        axis[idx].setBins(60, "plot_binning")
        
        bkg_hist[idx] = ROOT.RooDataHist(f"Minv_bkg_{cond}_{bin_key}_total", "bkg_total_histo", 
                                    ROOT.RooArgSet(axis[idx]), "plot_binning")
        for cat in bkg_categories:
            bkg_hist[idx].add(ws.data(f"Minv_bkg_{cond}_{bin_key}_{cat}"))

        background[idx] = ROOT.RooHistPdf(f"bkg_{cond}_{bin_key}", "bkg pdf", ROOT.RooArgSet(axis[idx]), bkg_hist[idx])
        

        # -----------------------
        #  Final models and fits
        # -----------------------

        n_events = histo_data[idx].sumEntries()

        nexp[idx] = ROOT.RooRealVar(
            f"nexp_{cond}_{bin_key}", f"nexp {cond}",
            n_events, n_events-5*(n_events**0.5), n_events+5*(n_events**0.5))

        '''
        Nsig[idx] = ROOT.RooRealVar(
            f"nsig_{cond}_{bin_key}", "#signal events",
            nexp[idx].getVal(), lims_num[idx][0]*nexp[idx].getVal(), nexp[idx].getMax())
        Nbkg[idx] = ROOT.RooRealVar(
            f"nbkg_{cond}_{bin_key}", "#background events",
            lims_num[idx][1]*nexp[idx].getVal()/2., 0.0, lims_num[idx][1]*nexp[idx].getVal())
        '''
 
        Nsig[idx] = ROOT.RooRealVar(f"nsig_{cond}_{bin_key}", f"nsig {cond}",
                                    0.9*nexp[idx].getVal(), 0.5, 1.5*nexp[idx].getVal())
        Nbkg[idx] = ROOT.RooRealVar(f"nbkg_{cond}_{bin_key}", f"nbkg {cond}",
                                    0.1*nexp[idx].getVal(), 0.5, 1.5*nexp[idx].getVal())           

        sum_func[idx] = ROOT.RooAddPdf(
                            f"sum_{cond}_{bin_key}", "Signal+Bkg",
                            ROOT.RooArgList(conv_pdf[idx], background[idx]),
                            ROOT.RooArgList(Nsig[idx], Nbkg[idx]))

        # sum_func[idx].setNormRange("fitRange")

        model[idx] = ROOT.RooAddPdf(
            sum_func[idx], f'PDF_{cond}_{bin_key}')

        model[idx].setNormRange("fitRange")

        results[idx] = model[idx].fitTo(histo_data[idx],
                                        # ROOT.RooFit.Extended(1),
                                        ROOT.RooFit.Range("fitRange"),
                                        ROOT.RooFit.Minimizer("Minuit2"),
                                        ROOT.RooFit.Strategy(fit_strategy[idx]),
                                        # ROOT.RooFit.MaxCalls(100000),
                                        ROOT.RooFit.Save(1),
                                        ROOT.RooFit.PrintLevel(verb)
                                        )

        status_chi2 = check_chi2(histo_data[idx], model[idx], results[idx])
        status = bool(status_chi2*fit_quality(results[idx], type_checks="benchmark"))

        print(status_chi2)
        print(status)

        low_nbkg = (Nbkg[idx].getVal() < 0.005*Nsig[idx].getVal())
        print(Nbkg[idx].getVal())
        print(Nsig[idx].getVal())
        print(low_nbkg)

        if (status is False) and refit_numbkg and low_nbkg:
            results[idx].Print()
            Nbkg[idx].setVal(0)
            Nbkg[idx].setConstant()
            print("REFITTING WITHOUT BACKGROUND")
            print("\n\n\n\n")
            results[idx] = model[idx].fitTo(histo_data[idx],
                                            # ROOT.RooFit.Extended(1),
                                            ROOT.RooFit.Range("fitRange"),
                                            ROOT.RooFit.Minimizer("Minuit2"),
                                            ROOT.RooFit.Strategy(fit_strategy[idx]),
                                            # ROOT.RooFit.MaxCalls(100000),
                                            ROOT.RooFit.Save(1),
                                            ROOT.RooFit.PrintLevel(verb)
                                            )
    
        if check_chi2(histo_data[idx], model[idx], results[idx]) is False:
            results[idx].SetTitle("Chi2_not_passed")

        fit_status[idx] = fit_quality(results[idx], type_checks="benchmark")

        results[idx].SetName(f"results_{cond}_{bin_key}")

    
    res_fail, res_pass = results

    eff, d_eff = efficiency_from_res(res_pass, res_fail)

    eff_mc, deff_mc = eval_efficiency(ws.data(f"Minv_mc_pass_{bin_key}").sumEntries(), 
                                      ws.data(f"Minv_mc_fail_{bin_key}").sumEntries(),
                                      sumw2_error(ws.data(f"Minv_mc_pass_{bin_key}")),
                                      sumw2_error(ws.data(f"Minv_mc_fail_{bin_key}")))
                
    scale_factor = eff/eff_mc
    d_scale_factor = scale_factor*((d_eff/eff)**2 + (deff_mc/eff_mc)**2)**0.5
    nsigma = (eff-eff_mc)/(d_eff**2 + deff_mc**2)**0.5 

    print(f"Efficiency: {eff} +- {d_eff}")
    print(f"Efficiency MC: {eff_mc} +- {deff_mc}")
    print(f"Scale factor: {scale_factor} +- {d_scale_factor}")


    
    # ----------------------------------
    #  Quality checks (and plots) - NEW
    # ----------------------------------

    status = bool(fit_status[0]*fit_status[1])

    figpath = "figs/fit_iso_indep_mcbkg_mergedbins" if status is True else "figs/check_fits/mcbkg_mergedbins"

    print(f"Num Bkg_histos = {bkg_hist[1].sumEntries()}, {bkg_hist[0].sumEntries()}")

    if figs:
        plot_objects = {
                "pass" : {
                    "axis" : axis[1],
                    "data" : histo_data[1],
                    "model" : model[1],
                    "res" : results[1],
                },
                "fail" : {
                    "axis": axis[0],
                    "data": histo_data[0],
                    "model": model[0],
                    "res": results[0],
                },
                "efficiency" : [eff, d_eff],
                "efficiency_mc" : [eff_mc, deff_mc],
                "scale_factor" : [scale_factor, d_scale_factor]
            }
        plot_fitted_pass_fail("indep", plot_objects, bin_key, pull=False, figpath=figpath)
        

    if status:
        ws.Import(model[0]), ws.Import(model[1])
        ws.Import(results[0]), ws.Import(results[1])                 
    else:
        res_pass.Print()
        res_pass.correlationMatrix().Print()
        print("****")
        res_fail.Print()
        res_fail.correlationMatrix().Print()
        print("****")
        print(res_pass.status(), res_fail.status())
        print(res_pass.covQual(), res_fail.covQual())
        print(res_pass.edm(), res_fail.edm())
        print('\n')


    return res_pass, res_fail, status



if __name__ == '__main__':

    
    ROOT.gROOT.SetBatch(True)
    ROOT.PyConfig.IgnoreCommandLineOptions = True

    t0 = time.time()

    type_eff = "iso"
    type_analysis = "indep"

    filename_data = "/scratchnvme/wmass/Steve_root_files/Standard_SF_files/tnp_iso_data_vertexWeights1_oscharge1.root"
    filename_mc = "/scratchnvme/wmass/Steve_root_files/Standard_SF_files/tnp_iso_mc_vertexWeights1_oscharge1.root"
    dirname_bkg = "/scratchnvme/rajarshi/Bkg_TNP_3D_Histograms/OS"
    
    bkg_categories= ["WW", "WZ", "ZZ", "TTSemileptonic", "Ztautau"]
    
    bkg_filenames = {}
    [bkg_filenames.update({cat : 
        f"{dirname_bkg}/tnp_{type_eff}_{cat}_vertexWeights1_oscharge1.root"}) for cat in bkg_categories]    

    lumi_scales = lumi_factors(type_eff, bkg_categories)

    lumi_scale_signal = lumi_scales.pop("Zmumu")

    bins = bin_dictionary("pt", "eta")

    import_dictionary = {
        "data" : filename_data,
        "mc" : {
            "filename": filename_mc,
            "lumi_scale" : lumi_scale_signal
        },
        "bkg" : {
            "filenames" : bkg_filenames,
            "lumi_scales" : lumi_scales
        }
    }

    workspace_name = "root_files/ws_indep_mcbkg.root"
    
    # ws = ws_init(import_dictionary, type_analysis, bins, binning("mass_60_120"))
    # ws.writeToFile(workspace_name)
    
    file = ROOT.TFile(workspace_name, "READ")
    ws = file.Get("w")

    Nproblems = 0

    for bin_key in bins.keys():



        _, bin_pt, bin_eta = bins[bin_key]

        existingFit = check_existing_fit(type_analysis, ws, bin_key)    

        if existingFit == 0:
            res_pass, res_fail, status = independent_efficiency(ws, bin_key, bkg_categories, refit_numbkg=True, figs=True)
        else :
            res_pass, res_fail = existingFit
            status = True

        if status is False:
            print(f"\nBin {bin_key} ({bin_pt}|{bin_eta}) has problems!\n")
            Nproblems += 1
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
    

    print(f"NPROBLEMS ={Nproblems}")


    t1 = time.time()

    print(f"TIME ELAPSED = {t1-t0}")