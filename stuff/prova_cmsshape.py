

import ROOT
import time
from utilities.dataset_utils import import_pdf_library


def enableBinIntegrator(func, num_bins):
    """
    Force numeric integration and do this numeric integration with the
    RooBinIntegrator, which sums the function values at the bin centers.
    """
    custom_config = ROOT.RooNumIntConfig(func.getIntegratorConfig())
    custom_config.method1D().setLabel("RooBinIntegrator")
    custom_config.getConfigSection(
        "RooBinIntegrator").setRealValue("numBins", num_bins)
    func.setIntegratorConfig(custom_config)
    func.forceNumInt(True)


def enableTrapezoidIntegrator(func, max_steps):
    """
    Force numeric integration and do this numeric integration with the
    RooBinIntegrator, which sums the function values at the bin centers.
    """
    custom_config = ROOT.RooNumIntConfig(func.getIntegratorConfig())
    custom_config.method1D().setLabel("RooIntegrator1D")
    custom_config.getConfigSection(
        "RooIntegrator1D").setRealValue("maxSteps", int(max_steps))
    func.setIntegratorConfig(custom_config)
    func.forceNumInt(True)


alpha = ROOT.RooRealVar("alpha", "alpha", 60.0, 40.0, 130.0)
beta = ROOT.RooRealVar("beta", "beta", 5, 0.1, 20)
gamma = ROOT.RooRealVar("gamma", "gamma", 0.07, 0.0001, 0.2)
peak = ROOT.RooRealVar("peak", "peak", 91.0)  # 88.0, 92.0)


def prova_cmsshape_default():
    """
    Attempt to fit an (extended) cmsshape on a dataset generated from the
    cmsshape itself. The fit is run with MINUIT/MIGRAD
    """

    import_pdf_library('RooCMSShape')

    axis = ROOT.RooRealVar("x", "x", 50, 130)
    axis.setBins(1000)

    background = ROOT.RooCMSShape("cmsshape_bkg", "CMSShape background",
                                  axis, alpha, beta, gamma, peak)

    axis.setRange("fitRange", 50, 130)

    mu = ROOT.RooRealVar("mu", "mu", 91, 85, 97)
    sigma = ROOT.RooRealVar("sigma", "sigma", 4, 0.5, 10)
    # background = ROOT.RooGaussian("gaus", "gaus", axis, mu, sigma)

    expected_num = ROOT.RooRealVar("nexp", "nexp", 10000, 9000, 11000)
    extend = ROOT.RooExtendPdf(
        "Extended", "extend", background, expected_num, rangeName="fitRange")

    data = background.generate({axis}, 10000)

    model = ROOT.RooExtendPdf(extend)
    '''
    enableBinIntegrator(background, 1000)
    enableBinIntegrator(model, 1000)
    '''

    res = model.fitTo(data,
                      # NumCPU=(10, 1),
                      Extended=True,
                      Range="fitRange",
                      #Constrain=ROOT.RooArgSet("nexp"),
                      Save=True,
                      # IntegrateBins=1e-6,
                      # Minimizer=("Minuit2", "migrad"),
                      PrintLevel=0)

    res.Print()
    res.correlationMatrix().Print()

    print(res.status())
    print(res.covQual())


def prova_cmsshape_minuit2(int_strategy, num):
    """
    Attempt to fit an (extended) cmsshape on a dataset generated from the
    cmsshape itself. The fit is run with MINUIT2, that has been made the
    default algorithm. A convergence is found (requiring less than 1h of
    operation) by using either the bin_integrator or the trapezoid_integrator
    with at most 12 steps.
    """
    t0 = time.time()

    import_pdf_library('RooCMSShape')

    axis = ROOT.RooRealVar("x", "x", 50, 130)
    axis.setBins(1000)
    axis.setRange("fitRange", 50, 130)

    background = ROOT.RooCMSShape("cmsshape_bkg", "CMSShape background",
                                  axis, alpha, beta, gamma, peak)

    if int_strategy == 'bin_center':
        enableBinIntegrator(background, num)
    elif int_strategy == 'trapezoid':
        enableTrapezoidIntegrator(background, num)
    else:
        pass

    expected_num = ROOT.RooRealVar("nexp", "nexp", 10000, 6000, 14000)
    model = ROOT.RooExtendPdf(
        "Extended", "extend", background, expected_num, rangeName="fitRange")

    data = background.generate({axis}, 10000)

    ROOT.Math.MinimizerOptions.SetDefaultMinimizer("Minuit2")
    ROOT.Math.MinimizerOptions.SetDefaultTolerance(2e-2)
    ROOT.Math.MinimizerOptions.SetDefaultMaxIterations(100000)
    ROOT.Math.MinimizerOptions.SetDefaultErrorDef(0.5)
    ROOT.Math.MinimizerOptions.SetDefaultStrategy(2)
    '''
    res1 = model.fitTo(data,
                       #NumCPU=(10, 1),
                       Extended=True,
                       Range="fitRange",
                       ExternalConstraints=ROOT.RooArgSet("nexp"),
                       Save=True,
                       Minimizer=("Minuit", "migrad"),
                       PrintLevel=-1)
    '''
    res2 = model.fitTo(data,
                       # NumCPU=(10, 1),
                       Extended=True,
                       Range="fitRange",
                       ExternalConstraints=ROOT.RooArgSet("nexp"),
                       MaxCalls=1000000,
                       Save=True,
                       PrintLevel=0)
    '''
    res1.Print()
    res1.correlationMatrix().Print()
    print(res1.status())
    print(res1.covQual())
    '''
    print("***********************************")
    res2.Print()
    res2.correlationMatrix().Print()
    print(res2.status())
    print(res2.covQual())

    t1 = time.time()
    print(f"TIME ELAPSED = {t1-t0}")


if __name__ == '__main__':

    # prova_cmsshape_default()

    prova_cmsshape_minuit2("analytical", 10000)
