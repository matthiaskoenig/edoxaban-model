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


class Chen2017b(EdoxabanSimulationExperiment):
    """Simulation experiment of Chen2017b."""

    bodyweight =  61.2  # kg

    #fraction_absorbed = {"Fasting": EdoxabanSimulationExperiment.fasting_map["fasted"]}
    colors = {
        "ED30": EdoxabanSimulationExperiment.dose_colors[30],
        "ED60": EdoxabanSimulationExperiment.dose_colors[60],
        "ED90": EdoxabanSimulationExperiment.dose_colors[90],
    }

    interventions = list(colors.keys())
    doses = {
        "ED30": 30,
        "ED60": 60,
        "ED90": 90
    }
    infos_pk = {
        "[Cve_edo]": "edoxaban",
        "[Cve_m4]": "M4",
        "Aurine_edo": "edoxaban_urine",
    }
    infos_pd = {
        "aPTT_change": "aPTT_change",
        "PT_change": "PT_change",
    }
    infos_pk_pd = {
        "aPTT_change": "aPTT_change",
        "PT_change": "PT_change",
    }

    def datasets(self) -> Dict[str, DataSet]:
        dsets = {}
        # timecourses
        for fig_id in ["Fig1", "Fig2", "Tab2A"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            label: str
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                if label.startswith("edoxaban"):
                    dset.unit_conversion("mean", 1 / self.Mr.edo)
                if label.startswith("M4_"):
                    dset.unit_conversion("mean", 1 / self.Mr.m4)
                dsets[label] = dset

        # scatter individual
        for fig_id in ["Fig3"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            y_label: str
            for y_label, df_label in df.groupby("y_label"):
                dset = DataSet.from_df(df_label, self.ureg)
                dset.unit_conversion("x", 1 / self.Mr.edo)
                dsets[y_label] = dset

        # scatter mean
        for fig_id in ["Fig3A"]:
            df = load_pkdb_dataframe(f"{self.sid}_{fig_id}", data_path=self.data_path)
            label: str
            for label, df_label in df.groupby("label"):
                dset = DataSet.from_df(df_label, self.ureg)
                x_label = df_label["x_label"].unique()[0]

                if x_label.startswith("edoxaban"):
                    dset.unit_conversion("x", 1 / self.Mr.edo)
                if x_label.startswith("M4_"):
                    dset.unit_conversion("x", 1 / self.Mr.m4)
                dsets[label] = dset

        return dsets

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}
        for intervention in self.interventions:
            dose = self.doses[intervention]
            tcsims[f"po_{intervention}"] = TimecourseSim(
                [Timecourse(
                    start=0,
                    end=50 * 60,  # [min]
                    steps=500,
                    changes={
                        **self.default_changes(),
                        "BW": Q_(self.bodyweight, "kg"),
                        "PODOSE_edo": Q_(dose, "mg"),
                    },
                )]
            )
        return tcsims

    def fit_mappings(self) -> Dict[str, FitMapping]:
        mappings = {}
        for intervention in self.interventions:
            # PK
            for ks, sid in enumerate(self.infos_pk):
                name = self.infos_pk[sid]
                mappings[f"fm_po_{intervention}_{name}"] = FitMapping(
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
                        self, task=f"task_po_{intervention}", xid="time", yid=sid,
                    ),
                    metadata=EdoxabanMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                        coadministration=Coadministration.NONE
                    ),
                )
            # PD
            for ks, sid in enumerate(self.infos_pd):
                name = self.infos_pd[sid]
                mappings[f"fm_po_{intervention}_{name}"] = FitMapping(
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
                        self, task=f"task_po_{intervention}", xid="time", yid=sid,
                    ),
                    metadata=EdoxabanMappingMetaData(
                        tissue=Tissue.PLASMA,
                        route=Route.PO,
                        application_form=ApplicationForm.TABLET,
                        dosing=Dosing.SINGLE,
                        health=Health.HEALTHY,
                        fasting=Fasting.FASTED,
                        coadministration=Coadministration.NONE
                    ),
                )

        return mappings

    def figures(self) -> Dict[str, Figure]:
        return {
            **self.figure_pk(),
            **self.figure_pd(),
            **self.figure_pk_pd(),
        }

    def figure_pk(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig1",
            num_cols=3,
            name=f"{self.__class__.__name__} (healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time), legend=True
        )
        plots[0].set_yaxis(self.label_edo_plasma, unit=self.unit_edo)
        plots[1].set_yaxis(self.label_m4_plasma, unit=self.unit_m4)
        plots[2].set_yaxis(self.label_edo_urine, unit=self.unit_edo_urine)

        for intervention in self.interventions:
            for ks, sid in enumerate(self.infos_pk):
                name = self.infos_pk[sid]
                # simulation
                plots[ks].add_data(
                    task=f"task_po_{intervention}",
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
                    linestyle="" if "urine" in name else "--",
                )

        return {
            fig.sid: fig,
        }

    def figure_pd(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig2",
            num_cols=2,
            name=f"{self.__class__.__name__} (healthy)",
        )
        plots = fig.create_plots(
            xaxis=Axis(self.label_time, unit=self.unit_time),
            legend=True,
        )
        plots[0].set_yaxis(self.labels["aPTT_change"], unit=self.units["aPTT_change"])
        plots[1].set_yaxis(self.labels["PT_change"], unit=self.units["PT_change"])

        for intervention in self.interventions:
            for ks, sid in enumerate(self.infos_pd):
                name = self.infos_pd[sid]
                # simulation
                plots[ks].add_data(
                    task=f"task_po_{intervention}",
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


    def figure_pk_pd(self) -> Dict[str, Figure]:
        fig = Figure(
            experiment=self,
            sid="Fig3",
            num_cols=2,
            name=f"{self.__class__.__name__} (healthy)",
        )
        plots = fig.create_plots()
        for k in range(2):
            plots[k].set_xaxis(self.label_edo_plasma, unit=self.unit_edo)

        plots[0].set_yaxis(self.labels["aPTT_change"], unit=self.units["aPTT_change"])
        plots[1].set_yaxis(self.labels["PT_change"], unit=self.units["PT_change"])

        for ks, sid in enumerate(self.infos_pk_pd):
            name = self.infos_pk_pd[sid]

            # individual data
            plots[ks].add_data(
                dataset=f"{name}",
                xid="x",
                yid="y",
                xid_sd=None,
                yid_sd=None,
                # count="count",
                label="ED",
                color="white",
                markeredgecolor="black",
                marker="o",
                linestyle="",
            )

            # simulation
            for ki, intervention in enumerate(self.interventions):
                plots[ks].add_data(
                    task=f"task_po_{intervention}",
                    xid="[Cve_edo]",
                    yid=sid,
                    label=intervention,
                    color=self.colors[intervention],
                )
                # mean data
                plots[ks].add_data(
                    dataset=f"edoxaban_{name}_{intervention}",
                    xid="x",
                    yid="y",
                    xid_sd="x_sd",
                    yid_sd="y_sd",
                    count="count",
                    label=intervention,
                    color=self.colors[intervention],
                    markeredgecolor="black",
                    marker="o",
                    linestyle="",
                )

        return {
            fig.sid: fig,
        }

if __name__ == "__main__":
    run_experiments(Chen2017b, output_dir=Chen2017b.__name__)
