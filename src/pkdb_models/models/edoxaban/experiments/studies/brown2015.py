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


class Brown2015(EdoxabanSimulationExperiment):
    """Simulation experiment of Brown2015.

    TODO: ETP data available, but not in model
    """

    colors = {
        "ED60": EdoxabanSimulationExperiment.dose_colors[60],
        "ED180": EdoxabanSimulationExperiment.dose_colors[180],
    }
    doses = {
        "ED60": 60,
        "ED180": 180
    }

    interventions = list(colors.keys())


    infos_pk = {
        "[Cve_edo]": "edoxaban"
    }
    infos_pd = {
        #"ETP (change relative)": "ETP",
        "PT_change": "PT"
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig2", "Fig3", "Fig4"]:
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
            tcsims[intervention] = TimecourseSim([
                Timecourse(
                    start=0,
                    end=73 * 60,
                    steps=800,
                    changes={
                        **self.default_changes(),
                        "PODOSE_edo": Q_(dose, "mg"),
                    },
                )
            ])

        return tcsims


    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for intervention in self.interventions:
            # PK
            for k, sid in enumerate(self.infos_pk.keys()):
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
                        self,
                        task=f"task_{intervention}",
                        xid="time",
                        yid=sid,
                    ),
                    metadata=EdoxabanMappingMetaData(
                        tissue=Tissue.URINE,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                    ),
                )
            # PD
            for ks, sid in enumerate(self.infos_pd):
                name = self.infos_pd[sid]
                mappings[f"fm_{intervention}_{name}"] = FitMapping(
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
                        fasting=Fasting.FASTED,
                        coadministration=Coadministration.NONE #COADMINISTRATION WITH RIF
                    ),
                )

        return mappings

    def figures(self) -> Dict[str, Figure]:
        return {
            **self.figure_pk(),
            **self.figure_pd(),
        }

    def figure_pk(self) -> Dict[str, Figure]:

        fig = Figure(
            experiment=self,
            sid="Fig1",
            name=f"{self.__class__.__name__} (healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_edo_plasma, unit=self.unit_edo)

        for intervention in self.interventions:
            for ks, sid in enumerate(self.infos_pk):
                name = self.infos_pk[sid]
                # simulation
                plots[ks].add_data(
                    task=f"task_{intervention}",
                    xid="time",
                    yid=sid,
                    label=intervention,
                    color=self.colors[intervention],
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
            fig.sid: fig,
        }

    def figure_pd(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig_PD",
            name=f"{self.__class__.__name__} (healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time),
            legend=True,
        )
        plots[0].set_yaxis(self.labels["PT_change"], unit=self.units["PT_change"])

        for intervention in self.interventions:
            for ks, sid in enumerate(self.infos_pd):
                name = self.infos_pd[sid]
                plots[ks].add_data(
                    task=f"task_{intervention}",
                    xid="time",
                    yid=sid,
                    label=f"{intervention}",
                    color=self.colors[intervention],
                )
                plots[ks].add_data(
                    dataset=f"{name}_{intervention}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=f"{intervention}",
                    color=self.colors[intervention],
                )

        return {
            fig.sid: fig,
        }

if __name__ == "__main__":
    run_experiments(Brown2015, output_dir=Brown2015.__name__)
