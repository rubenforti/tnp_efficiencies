
from multiprocessing import Pool
import ROOT
from itertools import repeat


x = ROOT.RooRealVar("x", "x", 0, 10)
mu = ROOT.RooRealVar("mu", "mu", 5)
sigma = ROOT.RooRealVar("sigma", "sigma", 1, 0.1, 5)
    


def doFit(ws, val):

    var = ROOT.RooRealVar(f"mu_{val}", f"mu_{val}", val, 0, 10)
    gauss_fit = ROOT.RooGaussian(f"gauss_{val}", f"gauss_{val}", x, var, sigma)
    data = ws.data("data")
    res = gauss_fit.fitTo(data, ROOT.RooFit.Save(), ROOT.RooFit.PrintLevel(-1))
    print(res.status(), res.covQual())
    return gauss_fit



if __name__ == '__main__':
    ROOT.gROOT.SetBatch(True)

    ws = ROOT.RooWorkspace("w", "w")
    gauss_gen = ROOT.RooGaussian("gauss", "gauss", x, mu, sigma)
    data = gauss_gen.generate(ROOT.RooArgSet(x), 1000)
    data.SetName("data")
    ws.Import(data)


    '''
    for i in range(10):
        var = ROOT.RooRealVar(f"mu_{i}", f"mu_{i}", i, 0, 10)
        gauss_fit = ROOT.RooGaussian(f"gauss_{i}", f"gauss_{i}", x, var, sigma)
        ws.Import(gauss_fit)

    '''


    pool = Pool(processes=6)
    map_res = pool.starmap(doFit, zip(repeat(ws), range(10)))

    '''
    print(type(map_res))
    print(map_res)


    frame = x.frame()
    data.plotOn(frame)

    map_res[0].plotOn(frame)
    frame.Draw()
    frame.SaveAs("prova.png")
    '''    

