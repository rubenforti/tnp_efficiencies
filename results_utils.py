"""
"""

import ROOT
import pickle
from utilities import fit_quality, eval_efficiency


class results_manager:
    """
    """

    def __init__(self, type_estimate):
        """
        """
        self._dict_results = {}
        self._analysis = type_estimate

    def open(self, filename):
        with open(filename, "rb") as file:
            self._dict_results = pickle.load(file)

    def write(self, filename):
        with open(filename, "wb") as file:
            pickle.dump(self._dict_results, file)
            file.close()

    def dictionary(self):
        return self._dict_results

    def add_result(self, bin_pt, bin_eta, check=False, *res):
        """
        """
        if check is True:
            quality = 1
            for result in res:
                quality *= result
            goodfit = bool(quality)
            print(goodfit)
        else:
            goodfit = True

        if goodfit and (f"{bin_pt},{bin_eta}" not in self._dict_results) and self._analysis == 'indep':
            res_pass, res_fail = res
            for par in res_pass.floatParsFinal():
                if par.GetName() == f'nsig_pass_({bin_pt}|{bin_eta})':
                    Npass = par.getVal()
                    sigma_Npass = par.getError()
            for par in res_fail.floatParsFinal():
                if par.GetName() == f'nsig_fail_({bin_pt}|{bin_eta})':
                    Nfail = par.getVal()
                    sigma_Nfail = par.getError()
            eff, d_eff = eval_efficiency(
                Npass, Nfail, sigma_Npass, sigma_Nfail)

            print(f'Measured efficiency is: {eff} +- {d_eff}')

            new_res = {
                f"{bin_pt},{bin_eta}": {
                    "efficiency": (eff, d_eff),
                    "fit_stat_pass": {
                        "parameters": res_pass.floatParsFinal(),
                        "cov_matrix": res_pass.covarianceMatrix(),
                        "global_correlation": res_pass.globalCorr(),
                        "migrad_status": res_pass.status(),
                        "cov_matrix_quality": res_pass.covQual(),
                        "EDM": res_pass.edm()
                        },
                    "fit_stat_fail": {
                        "parameters": res_fail.floatParsFinal(),
                        "cov_matrix": res_fail.covarianceMatrix(),
                        "global_correlation": res_fail.globalCorr(),
                        "migrad_status": res_fail.status(),
                        "cov_matrix_quality": res_fail.covQual(),
                        "EDM": res_fail.edm()
                        }
                    }
                }
            self._dict_results.update(new_res)

        elif goodfit and (f"{bin_pt},{bin_eta}" not in self._dict_results) and self._analysis == 'sim':
            new_res = {f"{bin_pt},{bin_eta}": {
                "parameters": res.floatParsFinal(),
                "cov_matrix": res.covarianceMatrix(),
                "global_correlation": res.globalCorr(),
                "migrad_status": res.status(),
                "cov_matrix_quality": res.covQual(),
                "EDM": res.edm()
                }
            }
            self._dict_results.update(new_res)
        else:
            pass

    def view_efficiency(self, bin_pt=0, bin_eta=0):
        """
        HAS TO BE EXTENDED FOR SIMULTANEOUS FITS
        """
        res = self._dict_results

        print(" Bins   Efficiency ")
        print("-------------------")
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
        HAS TO BE EXTENDED FOR SIMULTANEOUS FITS
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
        HAS TO BE EXTENDED FOR SIMULTANEOUS FITS
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
        HAS TO BE EXTENDED FOR SIMULTANEOUS FITS
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


if __name__ == '__main__':

    filename = 'indep_eff_results.pkl'

    res = results_manager()

    res.open(filename)

    probs = res.problematic_bins('migrad')

    print(len(probs))
