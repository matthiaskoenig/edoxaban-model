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

class Liu2022(EdoxabanSimulationExperiment):
    """Simulation experiment of Liu2022."""

    fasting_states = ["fasted", "fed"]
    fraction_absorbed = {
        "Fasted": EdoxabanSimulationExperiment.fasting_map["fasted"],
        "Fed": EdoxabanSimulationExperiment.fasting_map["fed"]
    }

    #fasting_states = ["Fasting", "Fed"]
    colors = {
        "fastedR": "black",
        "fastedT": "black",
        "fedR": "tab:blue",
        "fedT": "tab:blue",
    }
    interventions = list(colors.keys())

    bodyweights = {
        "fastedR": 69.2,
        "fastedT": 69.2,
        "fedR": 69.4,
        "fedT": 69.4,
    }

    infos_pk = {
        "[Cve_edo]": "edoxaban"
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Fig2"]:
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
            bw = self.bodyweights[intervention]
            fed = "fed" in intervention
            fasting_key = "fed" if fed else "fasted"
            fraction_absorbed = self.fasting_map[fasting_key]

            tcsims[intervention] = TimecourseSim([
                Timecourse(
                    start=0,
                    end=50 * 60,  # 30 hours in minutes
                    steps=500,
                    changes={
                        **self.default_changes(),
                        "BW": Q_(bw, "kg"),
                        "PODOSE_edo": Q_(60, "mg"),
                        "GU__F_edo_abs": Q_(fraction_absorbed, "dimensionless"),
                    },
                )
            ])

        return tcsims


    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for intervention in self.interventions:
            fed = "fed" in intervention
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
                        yid_sd= None,
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=f"task_{intervention}",
                        xid="time",
                        yid=sid,
                    ),
                    metadata=EdoxabanMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FED if fed else Fasting.FASTED,
                    ),
                )

        return mappings

    def figures(self) -> Dict[str, Figure]:
        return {
            **self.figure_pk()
        }

    def figure_pk(self) -> Dict[str, Figure]:

        fig = Figure(
            experiment=self,
            sid="Fig_PK",
            name=f"{self.__class__.__name__} (healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_edo_plasma, unit=self.unit_edo)

        for intervention in self.interventions:
            for k, sid in enumerate(self.infos_pk.keys()):
                name = self.infos_pk[sid]
                # simulation
                plots[k].add_data(
                    task=f"task_{intervention}",
                    xid="time",
                    yid=sid,
                    label=intervention,
                    color=self.colors[intervention],
                )

                # data
                plots[k].add_data(
                    dataset=f"{name}_{intervention}",
                    xid="time",
                    yid="mean",
                    yid_sd="mean_sd",
                    count="count",
                    label=intervention,
                    color=self.colors[intervention]
                )

        return {
            fig.sid: fig,
        }


if __name__ == "__main__":
    run_experiments(Liu2022, output_dir=Liu2022.__name__)