import ROOT
import sys
from utilities.base_lib import bin_dictionary
from utilities.dataset_utils import import_pdf_library
from utilities.res_tools import results_manager


gen_folder = "/scratch/rforti/tnp_efficiencies_results/tracking"

run_altSig = True

import_pdf_library("RooCMSShape")
bin_dict = bin_dictionary("pt_tracking", "eta")


filen = ROOT.TFile(gen_folder+"/benchmark/ws_tracking.root")
ws = filen.Get("w")

if run_altSig is False:
    
    res_nomi = results_manager("indep", "pt_tracking", "eta", import_ws=ws)
    file_sys = ROOT.TFile(gen_folder+"/BBlight_legacySettings/ws_tracking_BBlight.root")
    ws_syst = file_sys.Get("w")
    res_syst = results_manager("indep", "pt_tracking", "eta", import_ws=ws_syst)

else:
    #file_data = ROOT.TFile(gen_folder+"/../egm_tnp_results/tracking/mu_RunGtoH_mu_tracking_both.root")
    filen = ROOT.TFile(gen_folder+"/../egm_tnp_results/tracking/mu_RunGtoH_mu_tracking_both_nominalFit.root")
    file_sys = ROOT.TFile(gen_folder+"/../egm_tnp_results/tracking/mu_RunGtoH_mu_tracking_both_altSigFit.root")
    with open(gen_folder+"/../egm_tnp_results/tracking/allEfficiencies.txt", "r") as file_bmark:
        row_list = file_bmark.readlines()
    res_nomi = results_manager("indep", "pt_tracking", "eta", import_txt=row_list)
    res_syst = results_manager("indep", "pt_tracking", "eta", import_txt=row_list, altSig_check=True)

    peak = ROOT.RooRealVar("peakF", "peakF", 90.0)
    peak.setConstant() 


for bin_key, [gl_idx, pt_idx, eta_idx] in bin_dict.items():

    
    # if bin_key != "[24.0to35.0][-2.4to-2.3]": continue

    eff_nomi, d_eff_nomi = res_nomi.getEff(bin_key)
    eff_syst, d_eff_syst = res_syst.getEff(bin_key)

    syst, d_syst = eff_syst/eff_nomi, (eff_syst/eff_nomi)*((d_eff_syst/eff_syst)**2 + (d_eff_nomi/eff_nomi)**2)**0.5

    histo = ws.data(f"Minv_data_fail_{bin_key}")

    if run_altSig is False:
        
        axis = ws.var(f"x_fail_{bin_key}")

        pdf_syst = ws_syst.pdf(f"fitPDF_fail_{bin_key}")
        pars_syst = ws_syst.obj(f"results_fail_{bin_key}").floatParsFinal()
        nbkg_syst, nbkg_syst_err = pars_syst.find(f"nbkg_fail_{bin_key}").getVal(), pars_syst.find(f"nbkg_fail_{bin_key}").getError()
        
        #norm_tot_syst = nbkg_syst + pars_syst.find(f"nsig_fail_{bin_key}").getVal()

        pdf_nomi = ws.pdf(f"fitPDF_fail_{bin_key}")
        pars_nomi = ws.obj(f"results_fail_{bin_key}").floatParsFinal()
        nbkg_nomi, nbkg_nomi_err = pars_nomi.find(f"nbkg_fail_{bin_key}").getVal(), pars_nomi.find(f"nbkg_fail_{bin_key}").getError()
    
    else:

        axis = ROOT.RooRealVar("x", "x", 50, 130)
        axis.setBins(80)

        new_histo = ROOT.RooDataHist(f"Minv_data_fail_{bin_key}", f"Minv_data_fail_{bin_key}",
                                    ROOT.RooArgList(axis), "")
        
        for i in range(histo.numEntries()):
            bin_cont = histo.weight(i)
            new_histo.set(i, bin_cont, abs(bin_cont)**0.5)
        
        histo = new_histo

        gl_idx_corr = gl_idx+48-1
        res_key = ""
        
        if (gl_idx_corr < 100): 
            idx_string = "0"+str(gl_idx_corr)
        else:
            idx_string = str(gl_idx_corr)

        for key in filen.GetListOfKeys():
            if f"bin"+idx_string in key.GetName() and "_resF" in key.GetName():
                res_key = key.GetName()
                break
        if res_key == "":
            print(f"Error: no key found for bin {bin_key}")
            sys.exit()

        print(res_key, bin_key)


        ## Nominal fit
        pars_nomi = filen.Get(res_key).floatParsFinal()
        nbkg_nomi, nbkg_nomi_err = pars_nomi.find(f"nBkgF").getVal(), pars_nomi.find(f"nBkgF").getError()

        gaus_nomi = ROOT.RooGaussian("gaus_nomi", "gaus_nomi", axis, pars_nomi.find("meanF"), pars_nomi.find("sigmaF"))

        axis.setBinning(ROOT.RooUniformBinning(50, 130, 2000), "cache")

        mc_hist = ws.data(f"Minv_mc_fail_{bin_key}")
        mc_new_histo = ROOT.RooDataHist(f"Minv_mc_fail_{bin_key}", f"Minv_mc_fail_{bin_key}", ROOT.RooArgList(axis), "")
        for i in range(mc_hist.numEntries()):
            bin_cont = mc_hist.weight(i)
            mc_new_histo.set(i, bin_cont, abs(bin_cont)**0.5)

        mc_template = ROOT.RooHistPdf(f"mc_template_{bin_key}", f"mc_template_{bin_key}", 
                                      ROOT.RooArgSet(axis), mc_new_histo, 3)
        
        sigPdf_nomi = ROOT.RooFFTConvPdf("sigPdf_nomi", "sigPdf_nomi", axis, mc_template, gaus_nomi)
        sigPdf_nomi.setBufferFraction(0.1)

        if pars_nomi.getSize() == 7:
            alpha = pars_nomi.find(f"acmsF")
            beta = pars_nomi.find(f"betaF")
            gamma = pars_nomi.find(f"gammaF")
            bkgPdf_nomi = ROOT.RooCMSShape("cmsshape_bkg", "cmsshape_bkg", axis, alpha, beta, gamma, peak)
        else:
            c1 = pars_nomi.find(f"c1F")
            c2 = pars_nomi.find(f"c2F")
            c3 = pars_nomi.find(f"c3F")
            c4 = pars_nomi.find(f"c4F")
            bkgPdf_nomi = ROOT.RooChebychev("chebyshev_bkg", "chebyshev_bkg", axis, ROOT.RooArgList(c1, c2, c3, c4))

        pdf_nomi = ROOT.RooAddPdf("pdf_nomi", "pdf_nomi", 
                                  ROOT.RooArgList(sigPdf_nomi, bkgPdf_nomi), 
                                  ROOT.RooArgList(pars_nomi.find(f"nSigF"), pars_nomi.find(f"nBkgF")))
        
        for par in pdf_nomi.getParameters(histo):
            par.setConstant()
            par.Print()

        ## AltSig fit
        pars_syst = file_sys.Get(res_key).floatParsFinal()
        nbkg_syst, nbkg_syst_err = pars_syst.find(f"nBkgF").getVal(), pars_syst.find(f"nBkgF").getError()

        if pars_syst.getSize() == 7:
            alpha = pars_syst.find(f"acmsF")
            beta = pars_syst.find(f"betaF")
            gamma = pars_syst.find(f"gammaF")
            bkgPdf_alt = ROOT.RooCMSShape("cmsshape_bkg", "cmsshape_bkg", axis, alpha, beta, gamma, peak)
        else:
            c1 = pars_syst.find(f"c1F")
            c2 = pars_syst.find(f"c2F")
            c3 = pars_syst.find(f"c3F")
            c4 = pars_syst.find(f"c4F")
            bkgPdf_alt = ROOT.RooChebychev("chebyshev_bkg", "chebyshev_bkg", axis, ROOT.RooArgList(c1, c2, c3, c4))

        # dummy, just to make the following code easier
        pdf_syst = ROOT.RooAddPdf("pdf_syst", "pdf_syst", 
                                  ROOT.RooArgList(sigPdf_nomi, bkgPdf_alt), 
                                  ROOT.RooArgList(pars_syst.find(f"nSigF"), pars_syst.find(f"nBkgF")))

    c = ROOT.TCanvas(bin_key, bin_key, 1200, 900)
    c.cd()

    frame = axis.frame()



    pdf_nomi.plotOn(frame, Name=f"pdf_nomi",
                    LineColor=ROOT.kRed, 
                    LineStyle=ROOT.kDashed,
                    Normalization=histo.sumEntries())

    for comp in pdf_nomi.getComponents():
            if "bkg" in comp.GetName():
                bkg_set = ROOT.RooArgSet(comp)
                pdf_nomi.plotOn(frame, Name="pdf_bkg_nomi", Components=bkg_set,
                                LineColor=ROOT.kGreen, LineWidth=4,
                                Normalization=histo.sumEntries())
                
    if not run_altSig:
        pdf_syst.plotOn(frame, Name=f"pdf_syst",
                        LineColor=ROOT.kOrange, 
                        LineStyle=ROOT.kDashed,
                        Normalization=histo.sumEntries())
    
    for comp in pdf_syst.getComponents():
        if "bkg" in comp.GetName():
            bkg_set = ROOT.RooArgSet(comp)
            pdf_syst.plotOn(frame, Name="pdf_bkg_syst", Components=bkg_set,
                            LineColor=ROOT.kBlue, LineWidth=4,
                            Normalization=histo.sumEntries())
            
    histo.plotOn(frame)

    frame.SetTitle(f"Nomi and alt. fits - failing probes, bin {bin_key}")
    #frame.SetTitleOffset(0.8)
    frame.GetYaxis().SetTitle("Events/(1 GeV)")
    frame.GetXaxis().SetTitle("TP M_{inv} (GeV/c^{2})")
    
    frame.Draw()

    legend = ROOT.TLegend(0.65, 0.65, 0.9, 0.9)
    legend.AddEntry("", f"Syst. effect = {syst:.4f}", "")
    legend.AddEntry(histo, f"Data ({histo.sumEntries():.0f} entries)", "lep")
    legend.AddEntry(frame.findObject("pdf_nomi"), "Total PDF (nominal)", "l")
    if not run_altSig:
        legend.AddEntry(frame.findObject("pdf_syst"), "Total PDF (syst)", "l")
    else:
        legend.AddEntry("", "", "")
    legend.AddEntry(frame.findObject("pdf_bkg_nomi"), f"Ref. bkg ({nbkg_nomi:.0f} #pm {nbkg_nomi_err:.0f} evts)", "l")
    legend.AddEntry(frame.findObject("pdf_bkg_syst"), f"Alt. bkg ({nbkg_syst:.0f} #pm {nbkg_syst_err:.0f} evts)", "l")
    
    legend.Draw()

    folder = "BBlight_legacySettings" if run_altSig is False else "../egm_tnp_results/tracking"
    c.SaveAs(f"/scratch/rforti/tnp_efficiencies_results/tracking/{folder}/pdf_comparison/syst_{bin_key}.png")
    

    

