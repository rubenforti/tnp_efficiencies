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
from utilities.fit_utils import fit_quality, check_chi2



def independent_efficiency(ws, bin_key, bkg_shape, bkg_categories,
                           refit_numbkg=False, test_bkg=False, verb=-1, figs=False):
  
    path = os.path.dirname(__file__)
    ROOT.gSystem.cd(path)

    # ----------------------------------
    #  Initial parameters - MODIFY HERE
    # ----------------------------------
    NBINS = [5000, 5000]
    bufFractions = [0.2, 0.2]

    sumw2_error_option = False

    fit_strategy = [2, 2]

    mean_p = ROOT.RooRealVar(
        f"mean_pass_{bin_key}", "mean pass", 0, -5.0, 5.0)
    sigma_p = ROOT.RooRealVar(
        f"sigma_pass_{bin_key}", "sigma pass", 1, 0.1, 5.0)
    mean_f = ROOT.RooRealVar(
        f"mean_fail_{bin_key}", "mean fail", 0, -5.0, 5.0)
    sigma_f = ROOT.RooRealVar(
        f"sigma_fail_{bin_key}", "sigma fail", 0.5, 0.1, 5.0)

    if bkg_shape == "expo":
        tau_p = ROOT.RooRealVar(
            f"tau_pass_{bin_key}", "tau pass", 0.0, -5, 5)
        tau_f = ROOT.RooRealVar(
            f"tau_fail_{bin_key}", "tau fail", 0.0, -5, 5)
        tau = [tau_f, tau_p]
    elif bkg_shape == "mixed":
        pass
    elif bkg_shape == "cmsshape":
        pass
    else:
        print("REQUESTED BACKGROUND SHAPE IS NOT SUPPORTED")
        sys.exit()

    # ------------------------------
    #  Histograms and PDFs building
    # ------------------------------
    mean, sigma = [mean_f, mean_p], [sigma_f, sigma_p]
    axis, histo_pseudodata, histo_total_bkg, pdf_mc = [0, 0], [0, 0], [0, 0], [0, 0]
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
    
        axis[idx].setBins(60, "plot_binning")
        
        histo_total_bkg[idx] = ROOT.RooDataHist(f"Minv_bkg_{cond}_{bin_key}_total", "bkg_total_histo", 
                                                ROOT.RooArgSet(axis[idx]), "plot_binning")
        histo_pseudodata[idx] = ROOT.RooDataHist(f"Minv_pseudodata_{cond}_{bin_key}", "pseudodata_histo",
                                                 ROOT.RooArgSet(axis[idx]), "plot_binning")

        
        [histo_total_bkg[idx].add(ws.data(f"Minv_bkg_{cond}_{bin_key}_{cat}")) for cat in bkg_categories]

        histo_pseudodata[idx].add(ws.data(f"Minv_mc_{cond}_{bin_key}"))
        histo_pseudodata[idx].add(histo_total_bkg[idx])


        pdf_mc[idx] = ROOT.RooHistPdf(f"pdf_mc_{cond}_{bin_key}", "pdf MC", axis[idx],
                                      ws.data(f"Minv_mc_{cond}_{bin_key}"), 3)

        smearing[idx] = ROOT.RooGaussian(f"smearing_{cond}_{bin_key}",
                                         "Gaussian smearing", axis[idx], mean[idx], sigma[idx])

        conv_pdf[idx] = ROOT.RooFFTConvPdf(f"conv_{cond}_{bin_key}", "Convolution pdf",
                                           axis[idx], pdf_mc[idx], smearing[idx])
        conv_pdf[idx].setBufferFraction(bufFractions[idx])
        # conv_pdf[idx].setBufferStrategy(0)
        # conv_pdf[idx].setOperMode(3)

        if bkg_shape == "expo":
            background[idx] = ROOT.RooExponential(f"expo_bkg_{cond}_{bin_key}",
                                                  "Exponential background", axis[idx], tau[idx])
        elif bkg_shape == "mixed":
            pass
        elif bkg_shape == "cmsshape":
            pass
        else:
            print("REQUESTED BACKGROUND SHAPE IS NOT SUPPORTED")
            sys.exit()

        # -----------------------
        #  Final models and fits
        # -----------------------

        n_events = histo_pseudodata[idx].sumEntries()

        nexp[idx] = ROOT.RooRealVar(
            f"nexp_{cond}_{bin_key}", f"nexp {cond}",
            n_events, n_events-5*(n_events**0.5), n_events+5*(n_events**0.5))


        Nsig[idx] = ROOT.RooRealVar(f"nsig_{cond}_{bin_key}", f"nsig {cond}",
                                    nexp[idx].getVal(), 0.5, 5*nexp[idx].getVal())
        Nbkg[idx] = ROOT.RooRealVar(f"nbkg_{cond}_{bin_key}", f"nbkg {cond}",
                                    0.01*nexp[idx].getVal(), 0.5, 0.5*nexp[idx].getVal())           

        sum_func[idx] = ROOT.RooAddPdf(
                            f"sum_{cond}_{bin_key}", "Signal+Bkg",
                            ROOT.RooArgList(conv_pdf[idx], background[idx]),
                            ROOT.RooArgList(Nsig[idx], Nbkg[idx]))

        # sum_func[idx].setNormRange("fitRange")

        model[idx] = ROOT.RooAddPdf(sum_func[idx], f'PDF_{cond}_{bin_key}')

        model[idx].setNormRange("fitRange")

        results[idx] = model[idx].fitTo(histo_pseudodata[idx],
                                        # ROOT.RooFit.Extended(1),
                                        ROOT.RooFit.Range("fitRange"),
                                        ROOT.RooFit.Minimizer("Minuit2", "migrad"),
                                        ROOT.RooFit.SumW2Error(sumw2_error_option),
                                        ROOT.RooFit.Strategy(fit_strategy[idx]),
                                        # ROOT.RooFit.AsymptoticError(True),
                                        # ROOT.RooFit.MaxCalls(100000),
                                        ROOT.RooFit.Save(1),
                                        ROOT.RooFit.PrintLevel(verb)
                                        )

        status_chi2 = check_chi2(histo_pseudodata[idx], model[idx], results[idx], type="pearson")
        status = bool(status_chi2*fit_quality(results[idx], type_checks="benchmark"))

        print(status_chi2)
        print(status)

        low_nbkg = (Nbkg[idx].getVal() < 0.005*Nsig[idx].getVal())
        print(Nbkg[idx].getVal())
        print(Nsig[idx].getVal())
        print(low_nbkg)

        if (status is False) and refit_numbkg and low_nbkg and bkg_shape=='expo':
            results[idx].Print()
            tau[idx].setVal(1)
            tau[idx].setConstant()
            Nbkg[idx].setVal(0)
            Nbkg[idx].setConstant()
            print("REFITTING WITHOUT BACKGROUND")
            print("\n\n\n\n")
            results[idx] = model[idx].fitTo(histo_pseudodata[idx],
                                            # ROOT.RooFit.Extended(1),
                                            ROOT.RooFit.Range("fitRange"),
                                            ROOT.RooFit.Minimizer("Minuit2", "migrad"),
                                            ROOT.RooFit.SumW2Error(sumw2_error_option),
                                            # ROOT.RooFit.AsymptoticError(True),
                                            ROOT.RooFit.Strategy(fit_strategy[idx]),
                                            # ROOT.RooFit.MaxCalls(100000),
                                            ROOT.RooFit.Save(1),
                                            ROOT.RooFit.PrintLevel(verb)
                                            )
    
        if check_chi2(histo_pseudodata[idx], model[idx], results[idx], type="pearson") is False:
            results[idx].SetTitle("Chi2_not_passed")

        fit_status[idx] = fit_quality(results[idx], type_checks="benchmark")

        results[idx].SetName(f"results_{cond}_{bin_key}")


    eff, d_eff = efficiency_from_res(results[1], results[0])

    # ----------------------------------
    #  Quality checks (and plots) - NEW
    # ----------------------------------

    res_fail, res_pass = results


    eff_mc, deff_mc = eval_efficiency(ws.data(f"Minv_mc_pass_{bin_key}").sumEntries(), 
                                      ws.data(f"Minv_mc_fail_{bin_key}").sumEntries(),
                                      sumw2_error(ws.data(f"Minv_mc_pass_{bin_key}")),
                                      sumw2_error(ws.data(f"Minv_mc_fail_{bin_key}")))
                
    scale_factor = eff/eff_mc
    d_scale_factor = scale_factor*((d_eff/eff)**2 + (deff_mc/eff_mc)**2)**0.5
    nsigma = (eff-eff_mc)/(d_eff**2 + deff_mc**2)**0.5 


    status = bool(fit_status[0]*fit_status[1])

    figpath = "figs/fit_iso_indep_mergedbins_pseudodata" if status is True else "figs/check_fits/pseudodata"

    if figs:
        plot_objects = {
                "pass" : {
                    "axis" : axis[1],
                    "data" : histo_pseudodata[1],
                    "model" : model[1],
                    "res" : results[1],
                },
                "fail" : {
                    "axis": axis[0],
                    "data": histo_pseudodata[0],
                    "model": model[0],
                    "res": results[0],
                },
                "efficiency" : [eff, d_eff],
                "efficiency_mc" : [eff_mc, deff_mc],
                "scale_factor" : [scale_factor, d_scale_factor]
            }
        plot_fitted_pass_fail("indep", plot_objects, bin_key, pull=False, figpath=figpath)

    if status is True:
        ws.Import(histo_pseudodata[0]), ws.Import(histo_pseudodata[1])
        ws.Import(model[0]), ws.Import(model[1])
        ws.Import(results[0]), ws.Import(results[1])


    return res_pass, res_fail, status


'''


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

    bins = bin_dictionary("pt", "eta_8bins")

    import_dictionary = {
        # "data" : filename_data,
        "mc" : {
            "filename": filename_mc,
            "lumi_scale" : lumi_scale_signal
        },
        "bkg" : {
            "filenames" : bkg_filenames,
            "lumi_scales" : lumi_scales
        }
    }

    workspace_name = "root_files/ws_bkg_pseudodata.root"
    
    ws = ws_init(import_dictionary, type_analysis, bins, binning("mass_60_120"))
    ws.writeToFile(workspace_name)
    
    file = ROOT.TFile(workspace_name, "READ")
    ws = file.Get("w")

    resfile = ROOT.TFile("root_files/pseudodata_eff_comparison.root", "RECREATE")
    resfile.cd()


    print(len(bins.keys()))
 
    prob_bins = []

    # pool = Pool()
    # pool.map(independent_efficiency(), bins.keys())

    
    for bin_key in bins.keys():
        res_pass, res_fail = independent_efficiency(ws, bin_key, "expo", bkg_categories,
                                                    refit_numbkg=True, figs=True)



    t1 = time.time()

    print(f"TIME ELAPSED = {t1-t0}")
    
    '''
