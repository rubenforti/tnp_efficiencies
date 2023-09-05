"""
"""
import ROOT
import time
import json
import sys
from utilities.base_library import bin_dictionary, lumi_factors
from utilities.dataset_utils import ws_init
from make_fits import runFits

t0 = time.time()

ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

# -----------------------------------------------------------------------------
#  GENERAL SETTINGS
# ------------------
type_eff = "iso"
type_analysis = "indep"

localDatasets = False
load_McBkg = False

generateDatasets = True
default_fit_settings = True

binning_pt, binning_eta, binning_mass = "pt", "eta", "mass_60_120"
mergedbins_bkg = False
binning_pt_bkg, binning_eta_bkg = "pt_12bins", "eta_16bins"

bkg_types = ["WW", "WZ", "ZZ", "TTSemileptonic", "Ztautau"]

workspace_name = f"root_files/ws_{type_eff}_{type_analysis}_prova.root"
import_pdfs = True

savefigs = False
figpath = {"good": "figs/pseudodata_trigminus", "check": "figs/pseudodata_trigminus/check"} 


if mergedbins_bkg and (binning_pt != "pt" or binning_eta != "eta"):
    print("ERROR: Evaluation of background in merged bins for its comparison on data is allowed only wrt reco-bins of pt and eta for data")
    sys.exit()

# -----------------------------------------------------------------------------------------------------------
#  DATASET GENERATION
# --------------------
if generateDatasets:

    if load_McBkg:
        lumi_scales = lumi_factors(type_eff, bkg_types)
        lumi_scale_signal = lumi_scales.pop("Zmumu")
    else: 
        lumi_scale_signal = 1

    if localDatasets:
        filename_data = f"root_files/datasets/tnp_{type_eff}_data_vertexWeights1_oscharge1.root"
        filename_mc = f"root_files/datasets/tnp_{type_eff}_mc_vertexWeights1_oscharge1.root"
        dirname_bkg = "root_files/datasets"
    else:
        filename_data = f"/scratchnvme/wmass/Steve_root_files/Standard_SF_files/tnp_{type_eff}_data_vertexWeights1_oscharge1.root"
        filename_mc = f"/scratchnvme/wmass/Steve_root_files/Standard_SF_files/tnp_{type_eff}_mc_vertexWeights1_oscharge1.root"
        dirname_bkg = "/scratchnvme/rajarshi/Bkg_TNP_3D_Histograms/OS"

    import_dictionary = {
        "data" : filename_data, 
        "mc" : {"filename": filename_mc, "lumi_scale" : lumi_scale_signal}}
    ws = ws_init(import_dictionary, type_analysis, binning_pt, binning_eta, binning_mass)
    ws.writeToFile(workspace_name)

    if load_McBkg is True:
        bkg_filenames = {}
        [bkg_filenames.update({cat : 
            f"{dirname_bkg}/tnp_{type_eff}_{cat}_vertexWeights1_oscharge1.root"}) for cat in bkg_types]

        import_dict_bkg = {"bkg" : {"filenames" : bkg_filenames, "lumi_scales" : lumi_scales}}

        if mergedbins_bkg:
            ws_bkg = ws_init(import_dict_bkg, type_analysis, binning_pt_bkg, binning_eta_bkg, 
                             binning_mass, import_existing_ws=True, existing_ws_filename=workspace_name, 
                             altBinning_bkg=True)
        else:
            ws_bkg = ws_init(import_dict_bkg, type_analysis, binning_pt, binning_eta, 
                             binning_mass, import_existing_ws=True, existing_ws_filename=workspace_name, 
                             altBinning_bkg=False)
    
        ws_bkg.writeToFile(workspace_name)


# -----------------------------------------------------------------------------------------------------------
# FIT SETTINGS
# -------------
fit_settings_filename = "default_fit_settings.json" if default_fit_settings else "custom_fit_settings.json"

with open(fit_settings_filename) as file:
    fit_settings = json.load(file)

if ("idip" in type_eff) or ("trigger" in type_eff) or ("iso" in type_eff):
    fit_settings = fit_settings["idip_trig_iso"]
else:
    fit_settings = fit_settings[type_eff]

fit_settings.update({"type_analysis" : type_analysis})

if load_McBkg:
    fit_settings.update({"bkg_categories" : bkg_types})

# -----------------------------------------------------------------------------------------------------------
#  RUNNING FITS
# --------------
ws_fit = runFits(workspace_name, type_eff, bin_dictionary(binning_pt, binning_eta), fit_settings, 
                 fit_verb=-1, refit_numbkg=True, import_pdfs=import_pdfs, savefigs=savefigs, figpath=figpath)
ws_fit.writeToFile(workspace_name)

# -----------------------------------------------------------------------------------------------------------

t1 = time.time()

print(f"TEMPO = {t1-t0}")
