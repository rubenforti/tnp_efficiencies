"""
"""

import ROOT
import pickle


def fit_quality(res):
    """
    """
    check_migrad = (res.status() == 0)
    check_covm = (res.covQual() == 3)
    check_edm = (res.edm() < 1e-4)

    return bool(check_migrad*check_covm*check_edm)


def eval_efficiency(npass, nfail, sigma_npass, sigma_nfail):
    """
    """
    eff = npass/(npass+nfail)
    var1 = ((1-npass)**2)*(sigma_npass**2)
    var2 = (npass**2)*(sigma_nfail**2)
    sigma_eff = ROOT.TMath.Sqrt(var1+var2)/((npass+nfail)**2)

    return eff, sigma_eff


class res_manager_indep:
    """
    """

    def __init__(self):
        """
        """
        self._dict_results = {}

    def add_result(self, res_pass, res_fail, bin_pt, bin_eta, check=False):
        """
        """

        if check is True:
            goodfit = bool(fit_quality(res_pass)*fit_quality(res_fail))
            print(goodfit)
        else:
            goodfit = True

        if (goodfit is True) and (f"{bin_pt},{bin_eta}" not in self._dict_results):
            for par in res_pass.floatParsFinal():
                if par.GetName() == f'nsig_pass_({bin_pt},{bin_eta})':
                    Npass = par.getVal()
                    sigma_Npass = par.getError()
            for par in res_fail.floatParsFinal():
                if par.GetName() == f'nsig_fail_({bin_pt},{bin_eta})':
                    Nfail = par.getVal()
                    sigma_Nfail = par.getError()
            eff, d_eff = eval_efficiency(
                Npass, Nfail, sigma_Npass, sigma_Nfail)

            print(f'Measured efficiency is: {eff} +- {d_eff}')

            new_res = {
                f"{bin_pt},{bin_eta}": {
                    "efficiency": (eff, d_eff),
                    "fit_stat_pass": {
                        "migrad_status": res_pass.status(),
                        "parameters": res_pass.floatParsFinal(),
                        "cov_matrix": res_pass.covarianceMatrix(),
                        "cov_matrix_quality": res_pass.covQual(),
                        "global_correlation": res_pass.globalCorr(),
                        "EDM": res_pass.edm()
                        },
                    "fit_stat_fail": {
                        "migrad_status": res_fail.status(),
                        "parameters": res_fail.floatParsFinal(),
                        "cov_matrix": res_fail.covarianceMatrix(),
                        "cov_matrix_quality": res_fail.covQual(),
                        "global_correlation": res_fail.globalCorr(),
                        "EDM": res_fail.edm()
                        }
                    }
                }
            self._dict_results.update(new_res)
        else:
            pass

    def open(self, filename):
        """
        """
        with open(filename, "rb") as file:
            self._dict_results = pickle.load(file)

    def write(self, filename):
        """
        """
        with open(filename, "wb") as file:
            pickle.dump(self._dict_results, file)
            file.close()

    def dictionary(self):
        return self._dict_results

    def view_efficiencies(self, bin_pt=0, bin_eta=0):
        """
        """
        res = self._dict_results

        print(" Bins | Efficiency")
        if bin_pt == 0 and bin_eta == 0:
            for key in res.keys():
                eff = res[key]["efficiency"][0]
                d_eff = res[key]["efficiency"][1]
                if eff <= 1 and eff >= 0:
                    print(f'  {key} | {eff} +- {d_eff} ')
                else:
                    print(f'  {key} | {eff} +- {d_eff} !!!!!')
        else:
            eff = res[f"{bin_pt},{bin_eta}"]["efficiency"][0]
            d_eff = res[f"{bin_pt},{bin_eta}"]["efficiency"][1]
            if eff <= 1 and eff >= 0:
                print(f'  {key} | {eff} +- {d_eff} ')
            else:
                print(f'  {key} | {eff} +- {d_eff} !!!!!')

    def view_fits_statuses(self):
        """
        """
        res = self._dict_results
        print(" Bins | mig  covQual  EDM")
        print("     ")
        for key in res.keys():
            migr_p = res[key]["fit_stat_pass"]["migrad_status"]
            covq_p = res[key]["fit_stat_pass"]["cov_matrix_quality"]
            edm = res[key]["fit_stat_pass"]["EDM"]
            print(f"  {key} |  {migr_p}      {covq_p}      {edm}")
            migr_p = res[key]["fit_stat_fail"]["migrad_status"]
            covq_f = res[key]["fit_stat_fail"]["cov_matrix_quality"]
            edm = res[key]["fit_stat_fail"]["EDM"]
            print(f"  {key} |  {migr_p}      {covq_f}      {edm}")
            print("    ")

    def check_fit_status(self, bin_pt, bin_eta, conditions='all'):
        """
        """
        res = self._dict_results
        key = f"{bin_pt},{bin_eta}"

        eff = res[key]["efficiency"]
        migr_p = res[key]["fit_stat_pass"]["migrad_status"]
        covq_p = res[key]["fit_stat_pass"]["cov_matrix_quality"]
        edm_p = res[key]["fit_stat_pass"]["EDM"]
        migr_f = res[key]["fit_stat_fail"]["migrad_status"]
        covq_f = res[key]["fit_stat_fail"]["cov_matrix_quality"]
        edm_f = res[key]["fit_stat_fail"]["EDM"]

        if conditions == 'all':
            conditions = ['eff', 'migrad', 'covqual', 'edm']

        cond_list = []

        if 'eff' in conditions:
            cond_list.append(eff[0] > 0 and eff[0] < 1)
        if 'migrad' in conditions:
            cond_list.append(migr_p == 0 and migr_f == 0)
        if 'covqual' in conditions:
            cond_list.append(covq_p == 3 and covq_f == 3)
        if 'edm' in conditions:
            cond_list.append(edm_p < 1e-4 and edm_f < 1e-4)

        nfails = 0
        for cond in cond_list:
            if cond is False:
                nfails = nfails + 1

        return nfails

    def get_problematic_bins(self):
        """
        """
        res = self._dict_results
        bins = []
        for key in res.keys():
            bin_pt, bin_eta = key.split(',')
            status = self.check_fit_status(int(bin_pt), int(bin_eta))
            if status != 0:
                print(f'Bin {key} has {status} problems')
                bins.append(key)
        return bins


class res_manager_sim:
    """
    """

    def __init__(self):
        """
        """
        self._dict_results = {}

    def add_result(self, res, bin_pt, bin_eta):
        """
        """
        eff = 0
        for par in res.floatParsFinal():
            if par.GetName() == 'efficiency':
                eff = par.getVal()
                d_eff = par.getError()

        new_res = {f"{bin_pt},{bin_eta}": {
            "efficiency": (eff, d_eff),
            "migrad_status": res.status(),
            "parameters": res.floatParsFinal(),
            "cov_matrix": res.covarianceMatrix(),
            "cov_matrix_quality": res.covQual(),
            "global_correlation": res.globalCorr(),
            "EDM": res.edm()
            }
        }
        self._dict_results.update(new_res)

    def open(self, filename):
        """
        """
        with open(filename, "rb") as file:
            self._dict_results = pickle.load(file)

    def write(self, filename):
        """
        """
        with open(filename, "wb") as file:
            pickle.dump(self._dict_results, file)
            file.close()

    def view_efficiencies(self, bin_pt=0, bin_eta=0):
        """
        """
        print(" Bins | Efficiency")
        for key in res.keys():
            eff = res[key]["efficiency"]
            if eff > 1:
                print(f'  {key} | {eff}  !!!')
            else:
                print(f'  {key} | {eff}  ')

    def dictionary(self):
        return self._dict_results


if __name__ == '__main__':

    filename = 'indep_eff_results.pkl'

    res = res_manager_indep()

    res.open(filename)

    probs = res.problematic_bins('migrad')

    print(len(probs))
