"""
"""

import pickle


class results_manager:
    """
    """

    def __init__(self):
        """
        """
        self._dict_results = {}

    def add_result(self, res_pass, res_fail, eff, bin_pt, bin_eta):
        """
        """
        new_res = {f"{bin_pt},{bin_eta}": {
            "efficiency": eff,
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
        for key in self._dict_results.keys():
            eff = self._dict_results[key]["efficiency"][0]
            if eff > 1:
                print(f'  {key} | {eff}  !!!')
            else:
                print(f'  {key} | {eff}  ')

    def view_fits_status(self):
        """
        """
        res = self._dict_results
        print(" Bins | mig  covQual  EDM")
        for key in res.keys():
            migr_p = res[key]["fit_stat_pass"]["migrad_status"]
            covq_p = res[key]["fit_stat_pass"]["cov_matrix_quality"]
            edm = res[key]["fit_stat_pass"]["EDM"]
            print(f"  {key} |  {migr_p}      {covq_p}      {edm}")
            migr_p = res[key]["fit_stat_pass"]["migrad_status"]
            covq_f = res[key]["fit_stat_pass"]["cov_matrix_quality"]
            edm = res[key]["fit_stat_pass"]["EDM"]
            print(f"  {key} |  {migr_p}      {covq_f}      {edm}")
            print("    ")

    def problematic_bins(self, conditions=''):
        """
        """
        res = self._dict_results
        bins = []
        for key in res.keys():
            eff = res[key]["efficiency"]
            migr_p = res[key]["fit_stat_pass"]["migrad_status"]
            covq_p = res[key]["fit_stat_pass"]["cov_matrix_quality"]
            edm_p = res[key]["fit_stat_pass"]["EDM"]
            migr_p = res[key]["fit_stat_fail"]["migrad_status"]
            covq_f = res[key]["fit_stat_fail"]["cov_matrix_quality"]
            edm_f = res[key]["fit_stat_fail"]["EDM"]

            cond_list = []

            if 'eff' in conditions:
                cond_list.append((eff[0] < 0) or (
                    eff[0] > 1) or (eff[1] > 1e-2))
            if 'migrad' in conditions:
                cond_migrad = False
            if 'covqual' in conditions:
                cond_covm = False
            if 'edm' in conditions:
                cond_list.append((edm_p > 1e-2) or (edm_f > 1e-2))

            for cond in cond_list:
                if cond:
                    bins.append(key)
                    print(f"Problematic bin:  {key}")

        return bins
