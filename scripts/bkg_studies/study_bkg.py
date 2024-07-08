"""
"""
import ROOT
import utilities.base_lib as base_lib
from utilities.dataset_utils import ws_init
from utilities.bkg_utils import base_parser_bkg, analyze_bkg 


an_options = {
    "negweighted_bins" : ["is_SS_bkg"],
    "minv_distrib" : ["is_SS_bkg", "cmp_data", "cmp_signal", "cmp_bkgfit", "logscale"],
    "2D_distrib" : ["cmp_data", "cmp_signal", "cmp_totBkg", "projected"],
    "cross_cmp" : ["cmp_cat"]
}


ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

parser = base_lib.base_parser()
parser = base_parser_bkg(parser)

parser.add_argument("-s", "--bkg_to_study", type=str, nargs="+", choices=an_options.keys(), 
                    default="minv_distrib", help="Type of study to be performed")

parser.add_argument("--is_SS_bkg", action="store_true", help="Study the SS background")

parser.add_argument("--aux_ws", type=str, default="", help="Auxiliary workspace")

parser.add_argument("--cmp_data", action="store_true")
parser.add_argument("--cmp_signal", action="store_true")
parser.add_argument("--cmp_totBkg", action="store_true")
parser.add_argument("--cmp_bkgfit", action="store_true",
                    help="Compare the backgrounds with the expected one from the fit (has to be present in the auxiliary ws)")

parser.add_argument("--logscale", type=str, choices=["pass", "fail", "all", ""], default="all",
                    help="Set the log scale for the Minv plots of bkg")

parser.add_argument("--projected", action="store_true", help="Plot the 2D bkg distributions projected on the axes")

parser.add_argument("--cmp_cat", type=str, nargs=2, default=["bkg_SameCharge", "bkg_total"],
                    help="Background categories to be compared")


args = parser.parse_args()

base_lib.control_parsing(args) 

an_opt = {}
for an_type in args.bkg_study: 
    an_opt[an_type] = {opt : vars(args)[opt] for opt in an_options[an_type]}


# -----------------------------------------------------------------------------------------------------------
#  DATASET GENERATION
# --------------------

if args.generate_ws is True:
    
    import_categories = args.bkg_categories
    if "minv_distr" in args.study and args.cmp_data:
        import_categories += ["data"]
    if "minv_distr" in args.study and args.cmp_signal:
        import_categories += ["mc"]

    ws = ws_init(args.input, import_categories, args.eff, args.binning_pt, args.binning_eta, args.binning_mass, 
                 ch_set=args.charge_selection, 
                 do_OS_tracking=args.OS_tracking, 
                 add_SS_mc=args.study_SS, 
                 add_SS_bkg=(args.study_SS or args.add_SS_bkg), 
                 lightMode_bkg=args.lightMode_bkg)
                 
    ws.writeToFile(args.ws_filename)


analyze_bkg(args.eff, args.ws_filename, args.binning_pt, args.binning_eta, args.bkg_categories, args.bkg_study, an_opt,
            output_path=args.output_folder)
