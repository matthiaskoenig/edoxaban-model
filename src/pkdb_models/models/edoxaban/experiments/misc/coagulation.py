from typing import Dict

import numpy as np
from sbmlsim.plot import Axis, Figure, Plot
from sbmlsim.simulation import Timecourse, TimecourseSim

from pkdb_models.models.edoxaban.experiments.coagulation_experiment import CoagulationSimulationExperiment
from pkdb_models.models.edoxaban.helpers import run_experiments


class CoagulationExperiment(CoagulationSimulationExperiment):
    """Effect of edoxaban."""

    edo_values = np.linspace(0.0, 1, num=6)  # [µmole/l]  200-400 ng/ml = 0.36-0.72 µmole/l
    colors = ['black', '#eff3ff','#bdd7e7','#6baed6','#3182bd','#08519c']

    def simulations(self) -> Dict[str, TimecourseSim]:
        Q_ = self.Q_
        tcsims = {}

        baseline_changes = {}

        for k, edo in enumerate(self.edo_values):
            tcsims[f"coagulation_{k}"] = TimecourseSim([
                Timecourse(
                    start=0,
                    end=5 * 60,  # [min] # simulate 1 day
                    steps=500,
                    changes={
                        **self.default_changes(),
                        "[Cve_edo]": Q_(0, "µM"),
                        **baseline_changes,
                    },
                ),
                Timecourse(
                    start=0,
                    end=24 * 60,  # [min] # simulate 1 day
                    steps=2000,
                    changes={
                        "[Cve_edo]": Q_(edo, "µM")
                    },
                ),
                Timecourse(
                    start=0,
                    end=24 * 60,  # [min] # simulate 1 day
                    steps=2000,
                    changes={
                        "[Cve_edo]": Q_(0, "µM")
                    },
                ),
                ],
                time_offset=-5*60
            )

        return tcsims

    def figures(self) -> Dict[str, Figure]:

        info = [
            ("PT", 0),
            ("PT_change", 1),
            ("PT_ratio", 2),

            ("aPTT", 3),
            ("aPTT_change", 4),
            ("aPTT_ratio", 5),

            ("Xa_inhibition", 6),
        ]

        fig = Figure(
            experiment=self,
            sid=f"Fig_coagulation",
            num_rows=3,
            num_cols=3,
            name=f"Coagulation",
        )
        plots = fig.create_plots(xaxis=Axis("time", unit="hour"), legend=True)

        for sid, ksid in info:
            plots[ksid].set_yaxis(label=self.labels[sid], unit=self.units[sid])

            for k, edo in enumerate(self.edo_values):
                plots[ksid].add_data(
                    task=f"task_coagulation_{k}",
                    xid="time",
                    yid=sid,
                    label=f"{edo:.1f} nM",
                    color=self.colors[k],
                )

        return {fig.sid: fig}


if __name__ == "__main__":
    run_experiments(
       CoagulationExperiment, output_dir=CoagulationExperiment.__name__
    )
