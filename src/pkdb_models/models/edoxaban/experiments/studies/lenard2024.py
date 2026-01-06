from typing import Dict

from sbmlsim.data import DataSet, load_pkdb_dataframe
from sbmlsim.fit import FitMapping, FitData
from sbmlutils.console import console

from pkdb_models.models.edoxaban.experiments.base_experiment import (
    EdoxabanSimulationExperiment,
)
from pkdb_models.models.edoxaban.experiments.metadata import Tissue, Route, Dosing, ApplicationForm, Health, Health, \
    Fasting, EdoxabanMappingMetaData, Coadministration

from sbmlsim.plot import Axis, Figure
from sbmlsim.simulation import Timecourse, TimecourseSim

from pkdb_models.models.edoxaban.helpers import run_experiments


class Lenard2024(EdoxabanSimulationExperiment):
    """Simulation experiment of Lenard2024."""

    colors = {
        "API25, RIV25, EDO50, CLAR": "tab:blue",
        "API25, RIV25, EDO50": "black",
        "EDO60, CLAR": "tab:blue",
        "EDO60": "black",
    }

    interventions = list(colors.keys())

    doses =  {
        "API25, RIV25, EDO50": 0.050,  # [μg] (2 times dosing, first on day1, second on day 5)
        "API25, RIV25, EDO50, CLAR": 0.050,  # 50 µg
        "EDO60": 60, # [mg]
        "EDO60, CLAR": 60
    }

    labels = {
        "API25, RIV25, EDO50, CLAR": "EDO50, CLAR",
        "API25, RIV25, EDO50": "EDO50",
    }

    infos_pk = {
        "[Cve_edo]": "edoxaban"
    }


    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig2"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                if label.startswith("edoxaban"):
                        dset.unit_conversion("mean", 1 / self.Mr.edo)
                dsets[label] = dset
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        for intervention in self.interventions:
                dose = self.doses[intervention]
                tcsims[f"{intervention}"] = TimecourseSim(
                        [Timecourse(
                            start=0,
                            end=50 * 60,  # [min]
                            steps=500,
                            changes={
                                **self.default_changes(),
                                "PODOSE_edo": Q_(dose, "mg")
                            },
                        )])
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for intervention in self.interventions:
                # PK
                for ks, sid in enumerate(self.infos_pk):
                    name = self.infos_pk[sid]
                    mappings[f"fm_{name}_{intervention}"] = FitMapping(
                        self,
                        reference=FitData(
                            self,
                            dataset=f"{name}_{intervention}",
                            xid="time",
                            yid="mean",
                            yid_sd="mean_sd",
                            count="count",
                        ),
                        observable=FitData(
                            self, task=f"task_{intervention}", xid="time", yid=sid,
                        ),
                        metadata=EdoxabanMappingMetaData(
                            tissue=Tissue.PLASMA,
                            route=Route.PO,
                            application_form=ApplicationForm.TABLET,
                            dosing=Dosing.SINGLE,
                            health=Health.HEALTHY,
                            fasting=Fasting.NR,
                            coadministration=Coadministration.CLARITHROMYCIN if "CLAR" in intervention else Coadministration.NONE
                        ),
                    )

        return mappings

    def figures(self) -> Dict[str, Figure]:
        return {
            **self.figure_pk(),
        }

    def figure_pk(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig2",
            num_cols=2,
            name=f"{self.__class__.__name__} (healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_edo_plasma, unit=self.unit_edo)
        plots[1].set_yaxis(self.label_edo_plasma, unit=self.unit_edo)

        for intervention in self.interventions:
            name = "edoxaban"
            kp = 0 if "60" in intervention else 1
            # simulation
            plots[kp].add_data(
                task=f"task_{intervention}",
                xid="time",
                yid="[Cve_edo]",
                label=self.labels[intervention] if "50" in intervention else intervention,
                color=self.colors[f"{intervention}"],
            )
            # data
            plots[kp].add_data(
                dataset=f"{name}_{intervention}",
                xid="time",
                yid="mean",
                yid_sd="mean_sd",
                count="count",
                label=self.labels[intervention] if "50" in intervention else intervention,
                color=self.colors[f"{intervention}"],
            )

        return {
            fig.sid: fig
        }

if __name__ == "__main__":
    run_experiments(Lenard2024, output_dir=Lenard2024.__name__)
