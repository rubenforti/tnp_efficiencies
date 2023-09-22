import ROOT
from utilities.base_library import bin_dictionary

file = ROOT.TFile("results/benchmark_iso/ws_iso_indep_benchmark.root")

ws = file.Get("w")


bin_dict = bin_dictionary("pt", "eta")


for bin_key in bin_dict.keys():

    for flag in ["pass", "fail"]:

        if bin_key != "[60.0to65.0][2.0to2.1]": continue
 
        res  = ws.obj(f"results_{flag}_{bin_key}")
        # res.Print("v")
        pars = res.floatParsFinal()

        Nsig = pars.find(f"nsig_{flag}_{bin_key}")
        Nbkg = pars.find(f"nbkg_{flag}_{bin_key}")

        # pars.Print("v")
        # res.correlationMatrix().Print("v")
        # res.globalCorr().Print("v")


        if round((Nsig.getError() - Nsig.getVal()**0.5)*100/Nsig.getError(), 0)  < 0 : 
            print(bin_key, "nsig", Nsig.getError(), Nsig.getVal()**0.5, res.globalCorr(f"nsig_{flag}_{bin_key}"))
        
        if type(Nbkg) is ROOT.RooRealVar:
            if round((Nbkg.getError() - Nbkg.getVal()**0.5)*100/Nbkg.getError(), 0)  < 0 :
                print(bin_key, "nbkg", Nbkg.getError(), Nbkg.getVal()**0.5, res.globalCorr(f"nbkg_{flag}_{bin_key}"))
