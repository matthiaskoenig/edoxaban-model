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


class Chen2017a(EdoxabanSimulationExperiment):
    """Simulation experiment of Chen2017a."""


    colors = {
      # "EDO60": "black",
        "mED60": "black"
    }
    dose = 60

    interventions =  list(colors.keys())

    infos_pk = {
        "[Cve_edo]": "edoxaban",
        "[Cve_m4]": "M4"
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                if label.startswith("edoxaban"):
                        dset.unit_conversion("mean", 1 / self.Mr.edo)
                if label.startswith("M4_"):
                    dset.unit_conversion("mean", 1 / self.Mr.m4)
                dsets[label] = dset
        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        for intervention in self.interventions:
            dose = self.dose
            tc0 = Timecourse(
                start=0,
                end=48 * 60,  # [min]
                steps=500,
                changes={
                    **self.default_changes(),
                    "BW": Q_(61.8, "kg"),
                    "PODOSE_edo": Q_(dose, "mg")
                },
            )
            tc1 = Timecourse(
                start=0,
                end=24 * 60,  # [min]
                steps=500,
                changes={
                    "PODOSE_edo": Q_(dose, "mg")
                },
            )
            tc2 = Timecourse(
                start=0,
                end=50 * 60,  # [min]
                steps=500,
                changes={
                    "PODOSE_edo": Q_(dose, "mg")
                },
            )
            tcsims[f"{intervention}"] = TimecourseSim(
                [tc0] + [tc1 for _ in range(6)] + [tc2]
            )
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
                        dosing=Dosing.MULTIPLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.NR,
                        coadministration= Coadministration.NONE
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
            sid="Fig1",
            num_cols=2,
            #num_cols=2,
            name=f"{self.__class__.__name__} (healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_edo_plasma, unit=self.unit_edo)
        plots[1].set_yaxis(self.label_m4_plasma, unit=self.unit_m4)

        for intervention in self.interventions:
            for ks, sid in enumerate(self.infos_pk):
                    name = self.infos_pk[sid]

                    # simulation
                    plots[ks].add_data(
                        task=f"task_{intervention}",
                        xid="time",
                        yid=sid,
                        label=intervention,
                        color=self.colors[f"{intervention}"],
                    )
                    # data
                    plots[ks].add_data(
                        dataset=f"{name}_{intervention}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd",
                        count="count",
                        label=intervention,
                        color=self.colors[intervention],
                    )

        return {
            fig.sid: fig
        }


if __name__ == "__main__":
    run_experiments(Chen2017a, output_dir=Chen2017a.__name__)