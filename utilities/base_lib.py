"""
"""
import os, sys
import ROOT
import argparse
from array import array


binnings = {
    "pt": array('d', [24., 26., 28., 30., 32., 34., 36., 38., 40., 42., 44., 47., 50., 55., 60., 65.]),
    "eta" : array('d', [round(-2.4 + i*0.1, 2) for i in range(49)]),
    "mass_60_120" : array('d', [60 + i for i in range(61)]),
    "mass_50_130" : array('d', [50 + i for i in range(81)]),
    "pt_reco" : array('d', [24., 26., 30., 34., 38., 42., 46., 50., 55., 60., 65.]),
    "pt_tracking" : array('d', [24., 35., 45., 55., 65.]),

    "pt_singlebin" : array('d', [24., 65.]),
    "eta_singlebin" : array('d', [round(-2.4 + i*4.8, 2) for i in range(2)]),

    "mass_2GeV" : array('d', [60 + 2*i for i in range(31)]),
    "mass_3GeV" : array('d', [60 + 3*i for i in range(21)]),
    "mass_4GeV" : array('d', [60 + 4*i for i in range(16)]),
    "pt_12bins" : array('d', [24., 28., 30., 32., 34., 36., 38., 40., 44., 50., 55., 60., 65.]),
    "pt_9bins" : array('d', [24., 28., 32., 36., 40., 44., 50., 55., 60., 65.]),
    "pt_6bins" : array('d', [24., 30., 36., 42., 50., 55., 65.]),

    "eta_24bins" : array('d', [round(-2.4 + i*0.2, 2) for i in range(25)]),
    "eta_16bins" : array('d', [round(-2.4 + i*0.3, 2) for i in range(17)]),
    "eta_8bins" : array('d', [round(-2.4 + i*0.6, 2) for i in range(9)]),
    "eta_4bins" : array('d', [round(-2.4 + i*1.2, 2) for i in range(5)]),
    }

default_binnings_by_eff = {
    "reco" : ["pt_reco", "eta", "mass_60_120"],
    "tracking" : ["pt_tracking", "eta", "mass_50_130"],
    "idip" : ["pt", "eta", "mass_60_120"],
    "trigger" : ["pt", "eta", "mass_60_120"],
    "iso" : ["pt", "eta", "mass_60_120"]
}

bkg_categories = ["bkg_SameCharge",
                  "bkg_QCD",
                  "bkg_WW", "bkg_WZ", "bkg_ZZ", 
                  "bkg_TTSemileptonic", "bkg_TTFullyleptonic", "bkg_Ztautau",
                  "bkg_WplusJets", "bkg_WminusJets", "bkg_Zjets"
                  ]

default_steve_hists_folder = "/scratch/rforti/steve_histograms_2016/"


lumi_data = 16.8  # fb^-1, for 2016postVFP


BR_TAUToMU = 0.1739
BR_TAUToE = 0.1782
Z_TAU_TO_LEP_RATIO = (1.-(1. - BR_TAUToMU - BR_TAUToE)**2)
xsec_ZmmPostVFP = 2001.9

cross_sections_bkg = {
    # Unit = pb
    "QCD" : 238800,
    "WW" : 12.6,
    "WZ" : 27.59,  #5.4341,
    "ZZ" : 0.60,
    "TTFullyleptonic" : 88.29,
    "TTSemileptonic" : 366.34,
    "WplusJets" : 11765.9,
    "WminusJets" : 8703.87,
    "Ztautau" : xsec_ZmmPostVFP*Z_TAU_TO_LEP_RATIO,
    "Zjets" : xsec_ZmmPostVFP  # obtained from the signal MC file with the "reverse" gen-matching option
}

###############################################################################


def base_parser():
    """
    """
    parser = argparse.ArgumentParser(description='Command line basic parser of tnp_efficiencies')
    
    parser.add_argument('-e', '--eff', type=str, choices=['reco', 'tracking', 'idip', 'trigger', 'iso'], default='reco',
                        help='Type of efficiency to be studied')
    parser.add_argument('-g', '--generate_ws', action='store_true',
                        help="Generate the workspace")
    parser.add_argument('-i', '--input', type=str, default='',
                        help='Folder containing the input histograms (not used when the ws generation is off)')
    parser.add_argument('-o', '--output', type=str, default='./',
                        help='Folder where to store the results')
    parser.add_argument('-ws', '--ws_filename', type=str, default='',
                        help='Name of the root file that contains the workspace (or will contain, when generating it)')
    parser.add_argument('-y', '--year', type=str, choices=['2016', '2017', '2018'], default='2016')
    parser.add_argument('-ch', '--charge_selection', type=str, choices=['plus', 'minus', 'all', ''], default='all')
    parser.add_argument('-an', '--type_analysis', type=str, choices=['indep', 'sim', 'sim_sf'], default='indep',
                        help='Type of fitting analysis to be performed')
    parser.add_argument('--binning_pt', type=str, choices=[b_key for b_key in binnings.keys() if 'pt' in b_key])
    parser.add_argument('--binning_eta', type=str, choices=[b_key for b_key in binnings.keys() if 'eta' in b_key])
    parser.add_argument('--binning_mass', type=str, choices=[b_key for b_key in binnings.keys() if 'mass' in b_key])

    parser.add_argument("--ws_postfix", type=str, default="", 
                        help="Postfix to be added to the workspace name (when generating it)")
    parser.add_argument('--OS_tracking', action='store_true',
                        help='Study the tracking efficiency only in the OS region')
    parser.add_argument('--all_SA_fail', action='store_true',
                        help='Build the signal template for the failing probes by using all the probes with their SA inv. mass [tested for tracking step, but introduces a bias]')

    return parser

###############################################################################


def control_parsing(args):
    """
    Apply control conditions on the arguments returned by the parser
    """

    # Setting the default binning
    if args.binning_pt is None:   args.binning_pt   = default_binnings_by_eff[args.eff][0]
    if args.binning_eta is None:  args.binning_eta  = default_binnings_by_eff[args.eff][1]
    if args.binning_mass is None: args.binning_mass = default_binnings_by_eff[args.eff][2]

    # Control on output folder
    if not os.path.exists(args.output): os.makedirs(args.output)
    
    # Control on workspace filename. The resulting 'args.ws_filename' object is
    # the full path of the workspace file, ready to be used in the routine
    if args.ws_filename == "":
        eff_postfix = "" if args.charge_selection in ["all", ""] else args.charge_selection
        args.ws_filename = args.output+f"/ws_{args.eff}{eff_postfix}.root"
        if args.ws_postfix != "": 
            args.ws_filename = args.ws_filename.replace(".root", f"_{args.ws_postfix}.root")
    else:
        if not args.ws_filename.endswith(".root"): 
            args.ws_filename += ".root"
        if "/" not in args.ws_filename:
            args.ws_filename = args.output_folder+"/"+args.ws_filename
    
    # Charge selection
    args.charge_selection = ["plus", "minus"] if args.charge_selection=="all" else [args.charge_selection]

    # Control on the workspace existence
    if args.generate_ws and os.path.exists(args.ws_filename): 
        confirm_gen = input("Workspace already exists, do you REALLY want to generate a new one from scratch? (y/n) ")
        if confirm_gen != "y": args.generate_ws = False

    return args

###############################################################################


def lumi_factor(filepath, process):
    """
    Returns the lumi factor for the process in the given file
    """
    file = ROOT.TFile(filepath)
    wsum_histo = file.Get("weightSum")
    num_init = wsum_histo.Integral()

    if "mc" in filepath:
        xsection = xsec_ZmmPostVFP*1000 # has to be put in fb
    else:
        xsection = cross_sections_bkg[process]*1000
        
    lumi_process = num_init/xsection

    print(process, num_init, lumi_process)

    scale = lumi_data/lumi_process

    return scale

###############################################################################


def sumw2_error(histo):
    """
    Calculates the error on the integral of a RooDataHist as the (square root 
    of the) sum of the errors associated to each bin. The errors considered are
    the "SumW2", already stored in the RooDataHist.
    """
    variance = 0
    for i in range(0, histo.numEntries()):
        histo.get(i)
        variance += histo.weightError(ROOT.RooAbsData.SumW2)**2
    sum_error = variance**0.5

    return sum_error

###############################################################################


def eval_efficiency(npass, nfail, sigma_npass, sigma_nfail):
    """
    """
    eff = npass/(npass+nfail)
    var1 = (nfail**2)*(sigma_npass**2)
    var2 = (npass**2)*(sigma_nfail**2)
    sigma_eff = ROOT.TMath.Sqrt(var1+var2)/((npass+nfail)**2)

    return eff, sigma_eff

###############################################################################


def efficiency_from_res(res_pass, res_fail):
    """
    """
    pars_pass, pars_fail = res_pass.floatParsFinal(), res_fail.floatParsFinal()

    for par in pars_pass:
        if "nsig" in par.GetName():
            npass = par
    for par in pars_fail:
        if "nsig" in par.GetName():
            nfail = par

    eff, d_eff = eval_efficiency(npass.getVal(), nfail.getVal(), npass.getError(), nfail.getError())

    return eff, d_eff

###############################################################################


def eval_norm_corrected(Ndata, Nbkg_raw, f, df):
    """
    """
    Nsig_corr = (Ndata - Nbkg_raw)/(1-f)
    sigma_Nsig_corr = ROOT.TMath.Sqrt(Ndata + Nbkg_raw + (Nsig_corr*df)**2)/(1-f)
    
    Nbkg_corr = Nbkg_raw - (f*Nsig_corr)
    sigma_Nbkg_corr = ROOT.TMath.Sqrt(Nbkg_raw + (Nsig_corr*df)**2 + (f*sigma_Nsig_corr)**2)

    return Nsig_corr, sigma_Nsig_corr, Nbkg_corr, sigma_Nbkg_corr

###############################################################################


def import_pdf_library(*functions):
    """
    Imports the C++ defined functions into the ROOT interpreter. 
    """
    current_path = os.path.dirname(__file__)
    import_path = os.path.join(current_path, '..', 'libCpp')
    
    for function in functions:
        ctrl_head = ROOT.gInterpreter.Declare(f' #include "{import_path}/{function}.h"')
        ctrl_source = ROOT.gSystem.CompileMacro(f"{import_path}/{function}.cc", opt="k")

        if ctrl_head is not True: 
            sys.exit("ERROR in header loading")
        if ctrl_source != 1: 
            sys.exit("ERROR in sourcefile compiling and loading")

###############################################################################


def safeSystem(cmd, dryRun=False, quitOnFail=True):
    print(cmd)
    if dryRun is False:
        res = os.system(cmd)
        if res:
            print('-'*30)
            print("safeSystem(): error occurred when executing the following command. Aborting")
            print(cmd)
            print('-'*30)
            if quitOnFail:
                quit()
        return res
    else:
        return 0

