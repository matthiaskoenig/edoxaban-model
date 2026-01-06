"""Run all simulation experiments."""
import shutil
from pathlib import Path

from pkdb_models.models.edoxaban.helpers import run_experiments
from pkdb_models.models.edoxaban.experiments.studies import *
from pkdb_models.models.edoxaban.experiments.misc import *
from pkdb_models.models.edoxaban.experiments.scans.scan_parameters import EdoxabanParameterScan
from pkdb_models.models import edoxaban
from pymetadata.console import console
from pymetadata import log
from sbmlsim.plot import Figure

Figure.legend_fontsize = 8
Figure.fig_dpi = 300


logger = log.get_logger(__name__)

EXPERIMENTS = {
        "studies": [
            Bathala2012,
            Brown2015,
            Chen2017a,
            Chen2017b,
            Lenard2024,
            Lenard2025,
            Liu2022,
            Matsushima2013,
            Mendell2011,
            Mendell2015,
            Mendell2015a,
            Ogata2010,
            Parasrampuria2016,
            Parasrampuria2016b,
            Rohr2024,
        ],
        "bodyweight": [
        ],
        "dose_dependency": [
            Brown2015,
            Chen2017b,
            Ogata2010,
        ],
        "iv": [
            Matsushima2013,
        ],
        "ddi": [
            Matsushima2013,
            Lenard2024,
            Lenard2025,
            Mendell2015a,
        ],
        "multi": [
            Chen2017a,
            Ogata2010,
            Parasrampuria2016b,
        ],
        "food": [
            Mendell2011,
            Liu2022,
        ],
        "hepatic_impairment": [
            Mendell2015,
        ],
        "renal_impairment": [
        ],
        "misc": [
            DoseDependencyExperiment,
        ],
        "scan": [
            EdoxabanParameterScan,
        ]

}
EXPERIMENTS["all"] = EXPERIMENTS["studies"] + EXPERIMENTS["misc"] + EXPERIMENTS["scan"]


def run_simulation_experiments(
    selected: str = None,
    experiment_classes: list = None,
    output_dir: Path = None
) -> None:
    """Run edoxaban simulation experiments."""

    Figure.fig_dpi = 600
    Figure.legend_fontsize = 10

    # Determine which experiments to run
    if experiment_classes is not None:
        experiments_to_run = experiment_classes
        if output_dir is None:
            output_dir = edoxaban.RESULTS_PATH_SIMULATION / "custom_selection"
    elif selected:
        # Using the 'selected' parameter
        if selected not in EXPERIMENTS:
            console.rule(style="red bold")
            console.print(
                f"[red]Error: Unknown group '{selected}'. Valid groups: {', '.join(EXPERIMENTS.keys())}[/red]"
            )
            console.rule(style="red bold")
            return
        experiments_to_run = EXPERIMENTS[selected]
        if output_dir is None:
            output_dir = edoxaban.RESULTS_PATH_SIMULATION / selected
    else:
        console.print("\n[red bold]Error: No experiments specified![/red bold]")
        console.print("[yellow]Use selected='all' or selected='studies' or provide experiment_classes=[...][/yellow]\n")
        return

    # Run the experiments
    run_experiments(experiment_classes=experiments_to_run, output_dir=output_dir)

    # Collect figures into one folder
    figures_dir = output_dir / "_figures"
    figures_dir.mkdir(parents=True, exist_ok=True)
    for f in output_dir.glob("**/*.png"):
        if f.parent == figures_dir:
            continue
        try:
            shutil.copy2(f, figures_dir / f.name)
        except Exception as err:
            print(f"file {f.name} in {f.parent} fails, skipping. Error: {err}")
    console.print(f"Figures copied to: file://{figures_dir}", style="info")


if __name__ == "__main__":
    """Run experiments."""

    run_simulation_experiments(selected="all")
