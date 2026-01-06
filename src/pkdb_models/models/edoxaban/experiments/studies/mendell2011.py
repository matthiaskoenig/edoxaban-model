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

class Mendell2011(EdoxabanSimulationExperiment):
    """Simulation experiment of Mendell2011."""

    fasting_states = ["fasted", "fed"]
    fraction_absorbed = {
        "Fasted": EdoxabanSimulationExperiment.fasting_map["fasted"],
        "Fed": EdoxabanSimulationExperiment.fasting_map["fed"]
    }

    #fasting_states = ["Fasting", "Fed"]
    colors = {
        "japanFed": "tab:blue",
        "caucasFed": "tab:blue",
        "japanFast": "black",
        "caucasFast": "black",
    }
    groups = list(colors.keys())

    bodyweights = {
        "japanFed": 67.2,
        "caucasFed": 75.6,
        "japanFast": 67.2,
        "caucasFast": 75.6,
    }
    aptts = {
        "japanFed": 34.8,
        "caucasFed": 37.4,
        "japanFast": 35.5,
        "caucasFast": 33.7,
    }
    pts = {
        "japanFed": 12.4,
        "caucasFed": 12.5,
        "japanFast": 12.8,
        "caucasFast": 12.5,
    }

    infos_pk = {
        "[Cve_edo]": "edoxaban",
        "Aurine_edo": "edoxaban_urine"
    }

    infos_pd = {
        "aPTT": "aPTT",
        "PT": "PT"
    }
    ethnicities = {
        "japan": "Japanese",
        "caucas": "Caucasian",
    }
    fasting_states = ["Fast", "Fed"]

    aptt_baseline = 33


    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        for fig_id in ["Fig1", "Fig2", "Fig3", "Tab4A"]:
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

        for group in self.groups:
            bw = self.bodyweights[group]
            fed = "Fed" in group
            fasting_key = "fed" if fed else "fasted"
            fraction_absorbed = self.fasting_map[fasting_key]

            tcsims[group] = TimecourseSim([
                Timecourse(
                    start=0,
                    end=25 * 60,  # 30 hours in minutes
                    steps=500,
                    changes={
                        **self.default_changes(),
                        "BW": Q_(bw, "kg"),
                        "PODOSE_edo": Q_(60, "mg"),
                        "GU__F_edo_abs": Q_(fraction_absorbed, "dimensionless"),
                        "PT_ref": Q_(self.pts[group], "s"),
                        "aPTT_ref": Q_(self.aptts[group], "s"),
                    },
                )
            ])

        return tcsims


    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for group in self.groups:
            # PK
            for k, sid in enumerate(self.infos_pk.keys()):
                name = self.infos_pk[sid]
                mappings[f"fm_{name}_{group}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{name}_{group}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd" if "urine" in name else None,
                        count="count",
                    ),
                    observable=FitData(
                        self,
                        task=f"task_{group}",
                        xid="time",
                        yid=sid,
                    ),
                    metadata=EdoxabanMappingMetaData(
                        tissue=Tissue.URINE if "urine" in name else Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FED if "Fed" in group else Fasting.FASTED,
                    ),
                )
            # PD
            for ks, sid in enumerate(self.infos_pd):
                name = self.infos_pd[sid]
                mappings[f"fm_{group}_{name}"] = FitMapping(
                    self,
                    reference=FitData(
                        self,
                        dataset=f"{name}_{group}",
                        xid="time",
                        yid="mean",
                        yid_sd=None,
                        count="count",
                    ),
                    observable=FitData(
                        self, task=f"task_{group}", xid="time", yid=sid,
                    ),
                    metadata=EdoxabanMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FED if "Fed" in group else Fasting.FASTED,
                        coadministration=Coadministration.NONE
                    ),
                )

        return mappings

    def figures(self) -> Dict[str, Figure]:
        return {
            **self.figure_pk(),
            **self.figure_pd(),
        }

    def figure_pk(self) -> Dict[str, Figure]:

        figures = {}
        for ethnicity, ethnicity_label in self.ethnicities.items():
            fig = Figure(
                experiment=self,
                sid=f"Fig_PK_{ethnicity}",
                num_cols=2,
                name=f"{self.__class__.__name__} (healthy, {ethnicity_label})",
            )
            plots = fig.create_plots(
                xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
            )
            plots[0].set_yaxis(self.label_edo_plasma, unit=self.unit_edo)
            plots[1].set_yaxis(self.label_edo_urine, unit=self.unit_edo_urine)

            for fasting_state in self.fasting_states:
                group = f"{ethnicity}{fasting_state}"
                for k, sid in enumerate(self.infos_pk.keys()):
                    name = self.infos_pk[sid]

                    # simulation
                    plots[k].add_data(
                        task=f"task_{group}",
                        xid="time",
                        yid=sid,
                        label=fasting_state,
                        color=self.colors[f"{group}"],
                    )

                    # data
                    plots[k].add_data(
                        dataset=f"{name}_{group}",
                        xid="time",
                        yid="mean",
                        yid_sd="mean_sd" if "urine" in name else None,
                        count="count",
                        label=fasting_state,
                        color=self.colors[f"{group}"],
                        linestyle="" if "urine" in name else "--",
                    )

            figures[fig.sid] = fig

        return figures

    def figure_pd(self) -> Dict[str, Figure]:
        figures = {}
        for ethnicity, ethnicity_label in self.ethnicities.items():

            fig = Figure(
                experiment=self,
                sid=f"Fig_PD_{ethnicity}",
                num_cols=2,
                name=f"{self.__class__.__name__} (healthy, {ethnicity_label})",
            )
            plots = fig.create_plots(
                xaxis=Axis(self.label_time, unit=self.unit_time),
                legend=True,
            )
            plots[0].set_yaxis(self.labels["aPTT"], unit=self.units["aPTT"])
            plots[1].set_yaxis(self.labels["PT"], unit=self.units["PT"])

            for fasting_state in self.fasting_states:
                group = f"{ethnicity}{fasting_state}"
                for k, sid in enumerate(self.infos_pd):
                    name = self.infos_pd[sid]

                    plots[k].add_data(
                        task=f"task_{group}",
                        xid="time",
                        yid=sid,
                        label=fasting_state,
                        color=self.colors[group],
                    )
                    plots[k].add_data(
                        dataset=f"{name}_{group}",
                        xid="time",
                        yid="mean",
                        yid_sd=None,
                        count="count",
                        label=fasting_state,
                        color=self.colors[group],
                    )
            figures[fig.sid] = fig

        return figures

if __name__ == "__main__":
    run_experiments(Mendell2011, output_dir=Mendell2011.__name__)
