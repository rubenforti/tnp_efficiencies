import ROOT
import argparse
from utilities.base_lib import binnings
from utilities.results_manager import results_manager
from utilities.binning_utils import bin_dictionary

ROOT.gROOT.SetBatch(True)
ROOT.PyConfig.IgnoreCommandLineOptions = True

parser = argparse.ArgumentParser()
parser.add_argument("-nomi", "--nominal", type=str, help="File with the nominal results")
parser.add_argument("-alt", "--alternative", type=str, help="File with the alternative results")
parser.add_argument("-eff", "--efficiency", type=str, help="Efficiency type")
parser.add_argument("-o", "--output", type=str, help="Output file")
parser.add_argument("--altModel_txt", type=str, help="Alternative model")
#parser.add_argument("--name_nomi", type=str, help="Name of the nominal results")
#parser.add_argument("--name_alt", type=str, help="Name of the alternative results")

args = parser.parse_args()

binning_pt = f"pt_{args.efficiency}" if args.efficiency in ["reco", "tracking"] else "pt"
binning_eta = "eta"

# Mi faceva fatica cambiare a mano tutto
def binning(arg): return binnings[arg]

res_nominal = results_manager(args.nominal, "indep", binning_pt, binning_eta, altModel_check=args.altModel_txt)
res_alt = results_manager(args.alternative, "indep", binning_pt, binning_eta, altModel_check=args.altModel_txt)

h2d_eff_mc = ROOT.TH2D("EffMC2D", "EffMC2D", len(binning(binning_eta))-1, binning(binning_eta), len(binning(binning_pt))-1, binning(binning_pt))
h2d_eff_nomi = ROOT.TH2D("EffData2D", "EffData2D", len(binning(binning_eta))-1, binning(binning_eta), len(binning(binning_pt))-1, binning(binning_pt))
h2d_eff_alt = ROOT.TH2D("EffDataAltBkg2D", "EffDataAltBkg2D", len(binning(binning_eta))-1, binning(binning_eta), len(binning(binning_pt))-1, binning(binning_pt))
h2d_sf_nomi = ROOT.TH2D("SF2D_nominal", "SF2D_nominal", len(binning(binning_eta))-1, binning(binning_eta), len(binning(binning_pt))-1, binning(binning_pt))
h2d_sf_alt = ROOT.TH2D("SF2D_dataAltBkg", "SF2D_dataAltBkg", len(binning(binning_eta))-1, binning(binning_eta), len(binning(binning_pt))-1, binning(binning_pt))


for bin_key, [gl_idx, pt_idx, eta_idx] in bin_dictionary(binning_pt, binning_eta).items():

    eff_nomi, d_eff_nomi = res_nominal.getEff(bin_key)
    eff_alt, d_eff_alt = res_alt.getEff(bin_key)

    eff_mc, d_eff_mc = res_nominal.getEffMC(bin_key)

    sf, sf_alt = eff_nomi/eff_mc, eff_alt/eff_mc

    d_sf = sf * ((d_eff_nomi/eff_nomi)**2 + (d_eff_mc/eff_mc)**2)**0.5
    d_sf_alt = sf_alt * ((d_eff_alt/eff_alt)**2 + (d_eff_mc/eff_mc)**2)**0.5

    h2d_eff_mc.SetBinContent(eta_idx, pt_idx, eff_mc)
    h2d_eff_mc.SetBinError(eta_idx, pt_idx, d_eff_mc)

    h2d_eff_nomi.SetBinContent(eta_idx, pt_idx, eff_nomi)
    h2d_eff_nomi.SetBinError(eta_idx, pt_idx, d_eff_nomi)

    h2d_eff_alt.SetBinContent(eta_idx, pt_idx, eff_alt)
    h2d_eff_alt.SetBinError(eta_idx, pt_idx, d_eff_alt)

    h2d_sf_nomi.SetBinContent(eta_idx, pt_idx, sf)
    h2d_sf_nomi.SetBinError(eta_idx, pt_idx, d_sf)

    h2d_sf_alt.SetBinContent(eta_idx, pt_idx, sf_alt)
    h2d_sf_alt.SetBinError(eta_idx, pt_idx, d_sf_alt)


file_out = ROOT.TFile(args.output, "RECREATE")

file_out.cd()
h2d_eff_mc.Write()
h2d_eff_nomi.Write()
h2d_eff_alt.Write()
h2d_sf_nomi.Write()
h2d_sf_alt.Write()
file_out.Close()


print("\n\nSaved file for smoothing:", file_out.GetName())


                                



