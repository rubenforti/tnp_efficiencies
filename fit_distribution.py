"""
"""

import ROOT
from utilities import import_Steve_histos, makeAndSavePlot, pearson_chi2_eval, llr_test_bkg


def fit_without_bkg(axis, t, histo_data, histo_mc, flags, events_data, saveplot=False):
    """
    """

    if 'pass' in histo_data.GetTitle():
        id_flag = 'pass'
    elif 'fail' in histo_data.GetTitle():
        id_flag = 'fail'
    else:
        pass   


    pdf_mc = ROOT.RooHistPdf("pdf_mc", "pdf_mc", axis, histo_mc)
    mean = ROOT.RooRealVar("mean", "mean", 0, -2, 2)
    sigma = ROOT.RooRealVar("sigma", "sigma", 0.5, 0.001, 2)
    smearing = ROOT.RooGaussian("smearing", "smearing", axis, mean, sigma)

    axis.setBins(1000, "cache")
    conv_func = ROOT.RooFFTConvPdf("conv", "conv", axis, pdf_mc, smearing, 3)
    conv_func.setBufferFraction(0.1)

    res = conv_func.fitTo(histo_data, Extended=True, Save=True)

    pearson_chi2_eval(histo_data, conv_func, histo_data.numEntries(), res)

    if saveplot is True:
        makeAndSavePlot(axis, histo_data, conv_func, pull=False,
                        name=f"figs/fit_{t}/{id_flag}_{flags[0]}_{flags[1]}.pdf")

    return res


def fit_with_bkg(axis, t, histo_data, histo_mc, pdf_bkg, flags, events_data, saveplot=False):
    """
    """

    if 'pass' in histo_data.GetTitle():
        id_flag = 'pass'
    elif 'fail' in histo_data.GetTitle():
        id_flag = 'fail'
    else:
        pass

    pdf_mc = ROOT.RooHistPdf("pdf_mc", "pdf_mc", axis, histo_mc)
    mean = ROOT.RooRealVar("mean", "mean", 0, -2, 2)
    sigma = ROOT.RooRealVar("sigma", "sigma", 0.5, 0.001, 2)
    smearing = ROOT.RooGaussian("smearing", "smearing", axis, mean, sigma)

    '''
    tau = ROOT.RooRealVar("tau", "tau", -10, 0)
    expo = ROOT.RooExponential("expo", "expo", axis, tau)
    '''

    axis.setBins(1000, "cache")
    conv_func = ROOT.RooFFTConvPdf("conv", "conv", axis, pdf_mc, smearing, 3)
    conv_func.setBufferFraction(0.1)

    Nsig = ROOT.RooRealVar("nsig", "#signal events", 0, events_data)
    Nbkg = ROOT.RooRealVar("nbkg", "#background events", 0, events_data)

    sum_func = ROOT.RooAddPdf("sum", "sum", [conv_func, pdf_bkg], [Nsig, Nbkg])

    model = ROOT.RooAddPdf(sum_func)

    res = model.fitTo(histo_data, Extended=True, Save=True)

    pearson_chi2_eval(histo_data, model, histo_data.numEntries(), res)

    if saveplot is True:
        makeAndSavePlot(axis, histo_data, model, pull=False,
                        name=f"figs/fit_{t}/{id_flag}_{flags[0]}_{flags[1]}.pdf")

    null_bkg = llr_test_bkg(histo_data, model)
    present_bkg = not null_bkg
    print(f"Background is accepted? {present_bkg}")

    return res


if __name__ == '__main__':

    custom_pdfs = ['RooCBExGaussShape',
                   'RooDoubleCBFast', 'RooCMSShape', 'my_double_CB']
    # import_pdf_library(custom_pdfs[2])

    type_eff = ("sa", "global", "ID", "iso", "trigger", "veto")
    t = type_eff[3]

    idx_cond = 1  # One for "pass", zero for fail
    id_flag = "fail" if idx_cond == 0 else "pass"

    NBINS_MASS = 80

    h_data, h_mc, n_events, x = import_Steve_histos(t, [1], [1])

    print(
        f"Num events in data and mc = {n_events[0][idx_cond]}, {n_events[1][idx_cond]}")

    tau = ROOT.RooRealVar("tau", "tau", -10, 0)
    expo = ROOT.RooExponential("expo", "expo", x, tau)

    res_pass = fit_with_bkg(
        x, t, h_data[idx_cond], h_mc[idx_cond], expo, idx_cond, n_events[0][idx_cond])



    '''
    pdf_mc = ROOT.RooHistPdf("pdf_mc", "pdf_mc", x, h_mc[idx_cond])

    mean = ROOT.RooRealVar("mean", "mean", 0, -2, 2)
    sigma = ROOT.RooRealVar("sigma", "sigma", 0.5, 0.001, 2)
    smearing = ROOT.RooGaussian("smearing", "smearing", x, mean, sigma)

    tau = ROOT.RooRealVar("tau", "tau", -10, 0)
    expo = ROOT.RooExponential("expo", "expo", x, tau)

    x.setBins(1000, "cache")
    conv_func = ROOT.RooFFTConvPdf("conv", "conv", x, pdf_mc, smearing, 3)
    conv_func.setBufferFraction(0.1)

    Nsig = ROOT.RooRealVar("nsig", "#signal events", 0, n_events[0][idx_cond])
    Nbkg = ROOT.RooRealVar("nbkg", "#background events",
                           0, n_events[0][idx_cond])

    sum_func = ROOT.RooAddPdf("sum", "sum", [conv_func, expo], [Nsig, Nbkg])

    model = ROOT.RooAddPdf(sum_func)

    #pdf_sig = make_convolution(x, h_data[idx_cond], pdf_mc, smearing)
    #model = add_pdfs(x, h_data[idx_cond], pdf_sig, expo)

    res = model.fitTo(h_data[idx_cond], Extended=True, Save=True, Hesse=False)
    '''

    '''
    npars = NBINS_MASS - res.floatParsFinal().getSize()
    chi2_sqrtvar = (2*npars)**(1/2.)
    print(f"Expected chi2 pars: mu={npars}, sqrt(var)={chi2_sqrtvar}")
    chi2_obj = ROOT.RooChi2Var(
        "chi2", "chi2", model, h_data[idx_cond], Verbose=True)
    print(f"Measured chi2 = {chi2_obj.getVal()}")
    print(f"Distance in sigma = {(chi2_obj.getVal()-npars)/chi2_sqrtvar}")
    '''
