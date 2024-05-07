import ROOT
import sys
from utilities.base_lib import bin_dictionary, lumi_factor, binning
from utilities.dataset_utils import ws_init, get_roohist


base_folder = "/scratch/rforti/steve_histograms_2016/tracking"


filepaths = [f"{base_folder}/tnp_trackingplus_mc_vertexWeights1_genMatching1_oscharge0.root",
             f"{base_folder}/tnp_trackingminus_mc_vertexWeights1_genMatching1_oscharge0.root"]


def get_roohist_passSA(files, axis, bin_key, bin_pt, bin_eta,
                       global_scale=-1.0):

    if type(bin_pt) is int:  bin_pt = [bin_pt, bin_pt]
    if type(bin_eta) is int: bin_eta = [bin_eta, bin_eta]    

    numBins = axis.getBinning("x_binning").numBins()

    roohisto = ROOT.RooDataHist(f"Minv_mc_passSA_{bin_key}", f"Minv_mc_passSA_{bin_key}",
                                ROOT.RooArgSet(axis), "x_binning")
    
    tmp_histos, tmp_roohistos = [], []



    print(f"\nImporting histograms in bin {bin_key}")

    idx_file = 0

    histo_names = []

    for file in files:

        histo_name = f"pass_mu_DY_postVFP_alt"
        

        histo_names.append(histo_name)
        histo3d = file.Get(histo_name)

        # In the projection, option "e" is specified tofail_template_with_all_SA calculate the bin content errors in the new histogram for 
        # generic selection of bin_pt and bin_eta. Without it, all works as expected ONLY IF the projection is done 
        # on a single bin of pt-eta
        th1_histo = histo3d.ProjectionX(f"Histo_mc_passSA_{idx_file}", bin_pt[0], bin_pt[1], bin_eta[0], bin_eta[1], "e")
        print("Under/overflow events:", th1_histo.GetBinContent(0), th1_histo.GetBinContent(numBins+1))

        if global_scale > 0: th1_histo.Scale(global_scale)

        if th1_histo.Integral() < 0:
            print(f"ERROR: negative entries in TH1 in bin {bin_pt}-{bin_eta}")
            #sys.exit()

        th1_histo.Rebin(int(th1_histo.GetNbinsX()/numBins))

        tmp_histos.append(th1_histo)

        tmp_roohistos.append(ROOT.RooDataHist(f"Minv_mc_passSA_{bin_key}_{idx_file}", f"Minv_mc_passSA_{bin_key}",
                                              ROOT.RooArgList(axis), th1_histo))

        roohisto.add(tmp_roohistos[-1])

        idx_file += 1

    print("Beginning consistency check...")
    th1_integral = 0
    for tmp_histo in tmp_histos: th1_integral += tmp_histo.Integral()
    if round(th1_integral, 5) != round(roohisto.sumEntries(), 5):
        sys.exit(f"ERROR: TH1 and RooDataHist have different integrals in bin {bin_key}")
    for i in range(numBins):
        roohisto.get(i)
        th1_bincontent = 0
        for it_file in range(len(files)): 
            th1_bincontent += tmp_histos[it_file].GetBinContent(i+1)
        if round(roohisto.weight(i), 5) != round(th1_bincontent, 5):
            # print(roohisto.weight(i), th1_bincontent)
            sys.exit(f"ERROR: bin {i} has different content in TH1 and RooDataHist")
    print("Consistency check passed, RooDataHist imported correctly\n")
            
    return roohisto




if __name__ == "__main__":

    nbins_pt = 4
    nbins_eta = 48



    files = [ROOT.TFile(v, "READ") for v in filepaths]


    bins_mass = binning("mass_50_130")
    x_binning = ROOT.RooUniformBinning(bins_mass[0], bins_mass[-1], len(bins_mass)-1, "x_binning")

    axis = ROOT.RooRealVar(f"x", "TP M_inv", bins_mass[0], bins_mass[-1], unit="GeV/c^2")
    axis.setBinning(x_binning)

    bin_dict = bin_dictionary("pt_tracking", "eta")

    for b_key, [gl_index, pt_idx, eta_idx] in bin_dict.items():
        print(b_key, gl_index, pt_idx, eta_idx)
        roohist_passSA = get_roohist_passSA(files, axis, b_key, pt_idx, eta_idx, 
                                     global_scale=lumi_factor(filepaths[-1], "mc"))
        roohist_fail = get_roohist(files, "mc", "fail", axis, b_key, pt_idx, eta_idx,
                                   
                                   global_scale=lumi_factor(filepaths[-1], "mc"))
        roohist_allSA = get_roohist(files, "mc", "fail", axis, b_key, pt_idx, eta_idx,
                                    fail_template_with_all_SA=True,
                                    global_scale=lumi_factor(filepaths[-1], "mc"))
        
        roopdf_passSA = ROOT.RooHistPdf(f"pdf_passSA_{b_key}", f"pdf_passSA_{b_key}", ROOT.RooArgSet(axis), roohist_passSA)
        roopdf_fail = ROOT.RooHistPdf(f"pdf_fail_{b_key}", f"pdf_fail_{b_key}", ROOT.RooArgSet(axis), roohist_fail)
        roopdf_allSA = ROOT.RooHistPdf(f"pdf_allSA_{b_key}", f"pdf_allSA_{b_key}", ROOT.RooArgSet(axis), roohist_allSA)

        print(roohist_passSA.sumEntries(), roohist_fail.sumEntries(), roohist_allSA.sumEntries())
        

        c = ROOT.TCanvas()
        c.cd()
        frame = axis.frame()
        roopdf_passSA.plotOn(frame,
                             ROOT.RooFit.LineColor(ROOT.kCyan),
                             ROOT.RooFit.Normalization(roohist_fail.sumEntries(), ROOT.RooAbsReal.NumEvent)
                             )
        roohist_fail.plotOn(frame,
                            ROOT.RooFit.LineColor(ROOT.kBlue),
                            ROOT.RooFit.MarkerColor(ROOT.kBlue)
                            # ROOT.RooFit.Normalization(roohist_fail.sumEntries(), ROOT.RooAbsReal.NumEvent)
                            )
        roopdf_allSA.plotOn(frame,
                            ROOT.RooFit.LineColor(ROOT.kRed),
                            ROOT.RooFit.Normalization(roohist_fail.sumEntries(), ROOT.RooAbsReal.NumEvent)
                            )
        
        frame.Draw()

        
        pad_legend = ROOT.TPad("pad_legend", "pad_legend", 0.75, 0.75, 0.95, 0.9)
        pad_legend.SetMargin(0.1, 0.1, 0.1, 0.1)
        pad_legend.Draw()

        legend = ROOT.TLegend(0.1, 0.1, 0.9, 0.9)
        legend.SetFillColor(0)
        legend.SetTextSize(0.15)
        legend.SetTextAlign(12)
        legend.SetBorderSize(0)

        pad_legend.cd()

        legend.AddEntry(roopdf_passSA, "SA pass", "l")
        legend_obj = legend.GetListOfPrimitives().Last()
        legend_obj.SetLineColor(ROOT.kCyan)
        legend_obj.SetLineWidth(3)

        legend.AddEntry(roohist_fail, "SA fail", "lep")
        legend_obj = legend.GetListOfPrimitives().Last()
        legend_obj.SetLineColor(ROOT.kBlue)
        legend_obj.SetLineWidth(3)

        legend.AddEntry(roopdf_allSA, "SA pass+fail", "l")
        legend_obj = legend.GetListOfPrimitives().Last()
        legend_obj.SetLineColor(ROOT.kRed)
        legend_obj.SetLineWidth(3)

        legend.Draw()

        #c.SetLogy()

        
        c.SaveAs(f"/home/users/ruben/tnp_efficiencies/utilities_true/mc_templates_norm/mc_{b_key}.png")
        