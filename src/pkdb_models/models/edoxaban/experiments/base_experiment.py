"""
Reusable functionality for multiple simulation experiments.
"""
import pandas as pd
from collections import namedtuple
from typing import Dict
from pkdb_models.models.edoxaban import MODEL_PATH
from sbmlsim.experiment import SimulationExperiment
from sbmlsim.model import AbstractModel
from sbmlsim.task import Task

from pkdb_models.models.edoxaban.edoxaban_pk import calculate_edoxaban_pk, calculate_edoxaban_pd

# Constants for conversion
MolecularWeights = namedtuple("MolecularWeights", "edo m4 m6 mx")


class EdoxabanSimulationExperiment(SimulationExperiment):
    """Base class for all SimulationExperiments."""

    font = {"weight": "bold", "size": 22}
    scan_font = {"weight": "bold", "size": 15}
    tick_font_size = 15
    legend_font_size = 9
    suptitle_font_size = 25

    # labels
    label_time = "time"
    label_edo = "edoxaban"
    label_edo_total = "EDO+MET"
    label_m4 = "M4"
    label_m6 = "M6"
    label_mx = "MX"

    label_edo_plasma = label_edo + " plasma"
    label_m4_plasma = label_m4 + " plasma"
    label_m6_plasma = label_m6 + " plasma"
    label_mx_plasma = label_mx + " plasma"
    label_edo_total_plasma = label_edo_total + " plasma"

    label_edo_urine = label_edo + " urine"
    label_m4_urine = label_m4 + " urine"
    label_m6_urine = label_m6 + " urine"
    label_mx_urine = label_mx + " urine"
    label_edo_total_urine = label_edo_total + " urine"

    label_edo_feces = label_edo + " feces"
    label_m4_feces = label_m4 + " feces"
    label_m6_feces = label_m6 + " feces"
    label_mx_feces = label_mx + " feces"
    label_edo_total_feces = label_edo_total + " feces"


    labels: Dict[str, str] = {
        "time": "time",
        "[Cve_edo]": label_edo_plasma,
        "[Cve_m4]": label_m4_plasma,
        "[Cve_m6]": label_m6_plasma,
        "[Cve_mx]": label_mx_plasma,
        "[Cve_edo_total]": label_edo_total_plasma,

        "Aurine_edo": label_edo_urine,
        "Aurine_m4": label_m4_urine,
        "Aurine_m6": label_m6_urine,
        "Aurine_mx": label_mx_urine,
        "Aurine_edo_total": label_edo_total_urine,

        "Afeces_edo": label_edo_feces,
        "Afeces_m4": label_m4_feces,
        "Afeces_m6": label_m6_feces,
        "Afeces_mx": label_mx_feces,
        "Afeces_edo_total": label_edo_total_feces,

        "PT": "PT",  # "prothrombin time",
        "PT_change": "PT change", # "prothrombin time change",
        "PT_ratio": "PT ratio",  # "prothrombin time ratio",
        "aPTT": "aPTT",
        "aPTT_change": "aPTT change",
        "aPTT_ratio": "aPTT ratio",
        "Xa_inhibition": "Xa inhibition",
    }

    # units
    unit_time = "hr"
    unit_metabolite = "nM"
    unit_metabolite_urine = "µmole"
    unit_metabolite_feces = "µmole"

    unit_edo = unit_metabolite
    unit_m4 = unit_metabolite
    unit_m6 = unit_metabolite
    unit_mx = unit_metabolite

    unit_edo_urine = unit_metabolite_urine
    unit_m4_urine = unit_metabolite_urine
    unit_m6_urine = unit_metabolite_urine
    unit_mx_urine = unit_metabolite_urine

    unit_edo_feces = unit_metabolite_feces
    unit_m4_feces = unit_metabolite_feces
    unit_m6_feces = unit_metabolite_feces
    unit_mx_feces = unit_metabolite_feces

    units: Dict[str, str] = {
        "time": unit_time,
        "[Cve_edo]": unit_edo,
        "[Cve_m4]": unit_m4,
        "[Cve_m6]": unit_m6,
        "[Cve_mx]": unit_mx,
        "[Cve_edo_total]": unit_edo,

        "Aurine_edo": unit_edo_urine,
        "Aurine_m4": unit_m4_urine,
        "Aurine_m6": unit_m6_urine,
        "Aurine_mx": unit_mx_urine,
        "Aurine_edo_total": unit_edo_urine,

        "Afeces_edo": unit_edo_feces,
        "Afeces_m4": unit_m4_feces,
        "Afeces_m6": unit_m6_feces,
        "Afeces_mx": unit_mx_feces,
        "Afeces_edo_total": unit_edo_feces,

        "PT": "s",
        "PT_change": "s",
        "PT_ratio": "dimensionless",
        "aPTT": "s",
        "aPTT_change": "s",
        "aPTT_ratio": "dimensionless",
        "Xa_inhibition": "dimensionless",
    }
    dose_colors = {
        10: "grey",
        30: "tab:blue",
        60: "tab:orange",
        90: "tab:red",
        120: "tab:green",
        150: "tab:purple",
        180: "tab:brown"
    }


    # ----------- Fasting/food -----
    # food changes the fraction absorbed
    fasting_map = {  # GU__F_edo_abs
        "not reported": 0.82,  # assuming fasted state if nothing is reported
        "fasted": 0.82,
        "fed": 1.0,
    }
    fasting_colors = {
        "fasted": "black",
        "fed": "tab:red",
    }

    # ----------- Renal map --------------
    renal_map = {
        "Normal renal function": 101.0 / 101.0,  # 1.0,
        "Mild renal impairment": 50.0 / 101.0,  # 0.5
        "Moderate renal impairment": 35.0 / 101.0,  # 0.35
        "Severe renal impairment": 20.0 / 101.0,  # 0.20
        # "End stage renal disease": 10.5 / 101.0,  # 0.1
    }
    renal_colors = {
        "Normal renal function": "black",
        "Mild renal impairment": "#66c2a4",
        "Moderate renal impairment": "#2ca25f",
        "Severe renal impairment": "#006d2c",
        # "End stage renal disease": "#006d5e"
    }

    # ----------- Cirrhosis map --------------
    cirrhosis_map = {
        "Control": 0,
        "Mild cirrhosis": 0.3994897959183674,  # CPT A
        "Moderate cirrhosis": 0.6979591836734694,  # CPT B
        "Severe cirrhosis": 0.8127551020408164,  # CPT C
    }
    cirrhosis_colors = {
        "Control": "black",
        "Mild cirrhosis": "#74a9cf",  # CPT A
        "Moderate cirrhosis": "#2b8cbe",  # CPT B
        "Severe cirrhosis": "#045a8d",  # CPT C
    }

    def models(self) -> Dict[str, AbstractModel]:
        Q_ = self.Q_
        return {
            "model": AbstractModel(
                source=MODEL_PATH,
                language_type=AbstractModel.LanguageType.SBML,
                changes={},
            )
        }

    @staticmethod
    def _default_changes(Q_):
        """Default changes to simulations."""

        changes = {
            # 	>>> !Optimal parameter 'LI__EDO2M6_f' within 5% of lower bound! <<<
            # 	>>> !Optimal parameter 'LI__EDO2MX_f' within 5% of upper bound! <<<
            # 	>>> !Optimal parameter 'LI__MXEXBI_k' within 5% of upper bound! <<<
            'GU__Ka_dis_edo': Q_(0.3559517961123874, '1/hr'),  # [0.001 - 100]
            'GU__EDOABS_k': Q_(0.021312888266795352, '1/min'),  # [0.001 - 10]
            'LI__EDO2M4_Vmax': Q_(0.02467742534326666, '1/min'),  # [0.0001 - 10]
            'LI__EDO2M6_f': Q_(0.20000000005912993, 'dimensionless'),  # [0.2 - 0.5]
            'LI__EDO2MX_f': Q_(9.999946086529825, 'dimensionless'),  # [2 - 10]
            'LI__MXEXBI_k': Q_(9.999998121216769e-05, '1/min'),  # [1e-06 - 0.0001]
            'KI__EDOEX_k': Q_(1.5766419182037374, '1/min'),  # [0.0001 - 10]
            'KI__M4EX_k': Q_(1.7146667975857133, '1/min'),  # [0.0001 - 10]
            'KI__M6EX_k': Q_(0.5765035347258279, '1/min'),  # [0.0001 - 10]
            'KI__MXEX_k': Q_(2.8819072042706586, '1/min'),  # [0.0001 - 10]

            'Emax_PT': Q_(3.6762994520462455, 'dimensionless'),  # [0.1 - 10]
            'EC50_edo_PT': Q_(0.0027062726583800974, 'mM'),  # [1e-07 - 0.01]
            'Emax_aPTT': Q_(1.9397792689345885, 'dimensionless'),  # [0.1 - 10]
            'EC50_edo_aPTT': Q_(0.0008459875594125807, 'mM'),  # [1e-07 - 0.01]
        }

        return changes

    def default_changes(self: SimulationExperiment) -> Dict:
        """Default changes to simulations."""
        return EdoxabanSimulationExperiment._default_changes(Q_=self.Q_)

    def tasks(self) -> Dict[str, Task]:
        if self.simulations():
            return {
                f"task_{key}": Task(model="model", simulation=key)
                for key in self.simulations()
            }
        return {}

    def data(self) -> Dict:
        self.add_selections_data(
            selections=[
                "time",
                "[Cve_edo]",
                "[Cve_m4]",
                "[Cve_m6]",
                "[Cve_mx]",
                "[Cve_edo_total]",

                "Aurine_edo",
                "Aurine_m4",
                "Aurine_m6",
                "Aurine_mx",
                "Aurine_edo_total",

                "Afeces_edo",
                "Afeces_m4",
                "Afeces_m6",
                "Afeces_mx",
                "Afeces_edo_total",

                # cases
                'KI__f_renal_function',
                'f_cirrhosis',

                "PT",
                "PT_change",
                "PT_ratio",
                "aPTT",
                "aPTT_change",
                "aPTT_ratio",
                "Xa_inhibition",

                "PODOSE_edo",
                "GU__F_edo_abs",
                "BW",
            ]
        )
        return {}

    @property
    def Mr(self):
        return MolecularWeights(
            edo=self.Q_(548.058, "g/mole"),
            m4=self.Q_(521.0, "g/mole"),
            m6=self.Q_(534.0, "g/mole"),
            mx=self.Q_(521.0, "g/mole"),
        )

    # --- Pharmacokinetic parameters ---
    pk_labels = {
        "auc": "AUCend",
        "aucinf": "AUC",
        "cl": "Total clearance",
        "cl_renal": "Renal clearance",
        "cl_hepatic": "Hepatic clearance",
        "cmax": "Cmax",
        "thalf": "Half-life",
        "kel": "kel",
        "vd": "vd",
    }

    pk_units = {
        "auc": "µmole/l*hr",
        "aucinf": "µmole/l*hr",
        "cl": "ml/min",
        "cl_renal": "ml/min",
        "cl_hepatic": "ml/min",
        "cmax": "µmole/l",
        "thalf": "hr",
        "kel": "1/hr",
        "vd": "l",
    }

    def calculate_edoxaban_pk(self, scans: list = []) -> Dict[str, pd.DataFrame]:
       """Calculate pk parameters for simulations (scans)"""
       pk_dfs = {}
       if scans:
           for sim_key in scans:
               xres = self.results[f"task_{sim_key}"]
               df = calculate_edoxaban_pk(experiment=self, xres=xres)
               pk_dfs[sim_key] = df
       else:
           for sim_key in self._simulations.keys():
               xres = self.results[f"task_{sim_key}"]
               df = calculate_edoxaban_pk(experiment=self, xres=xres)
               pk_dfs[sim_key] = df
       return pk_dfs

    def calculate_edoxaban_pd(self, scans: list = []) -> Dict[str, pd.DataFrame]:
       """Calculate pd parameters for simulations (scans)"""
       pd_dfs = {}
       if scans:
           for sim_key in scans:
               xres = self.results[f"task_{sim_key}"]
               df = calculate_edoxaban_pd(experiment=self, xres=xres)
               pd_dfs[sim_key] = df
       else:
           for sim_key in self._simulations.keys():
               xres = self.results[f"task_{sim_key}"]
               df = calculate_edoxaban_pd(experiment=self, xres=xres)
               pd_dfs[sim_key] = df
       return pd_dfs