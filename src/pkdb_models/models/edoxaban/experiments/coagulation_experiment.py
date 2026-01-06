"""Reusable functionality for multiple simulation experiments."""

from typing import Dict

from sbmlsim.experiment import SimulationExperiment
from sbmlsim.model import AbstractModel
from sbmlsim.task import Task

from pkdb_models.models.edoxaban import MODEL_BASE_PATH


class CoagulationSimulationExperiment(SimulationExperiment):
    """Base class for all SimulationExperiments."""

    font = {"weight": "bold", "size": 22}
    scan_font = {"weight": "bold", "size": 15}
    tick_font_size = 15
    legend_font_size = 13
    suptitle_font_size = 40

    units: Dict[str, str] = {
        "time": "hr",
        "PT": "s",
        "PT_change": "s",
        "PT_ratio": "dimensionless",
        "aPTT": "s",
        "aPTT_change": "s",
        "aPTT_ratio": "dimensionless",
        "Xa_inhibition": "dimensionless",
    }
    labels: Dict[str, str] = {
        "time": "time",
        "PT": "prothrombin time",
        "PT_change": "prothrombin time change",
        "PT_ratio": "prothrombin time ratio",
        "aPTT": "aPTT",
        "aPTT_change": "aPTT",
        "aPTT_ratio": "aPTT ratio",
        "Xa_inhibition": "Xa inhibition",
    }

    def models(self) -> Dict[str, AbstractModel]:
        Q_ = self.Q_
        return {
            "model": AbstractModel(
                source=MODEL_BASE_PATH / "edoxaban_coagulation.xml",
                language_type=AbstractModel.LanguageType.SBML,
                changes={},
            )
        }

    @staticmethod
    def _default_changes(Q_):
        """Default changes to simulations."""

        changes = {}

        return changes

    def default_changes(self: SimulationExperiment) -> Dict:
        """Default changes to simulations."""
        return CoagulationSimulationExperiment._default_changes(Q_=self.Q_)

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
                "PT",
                "PT_change",
                "PT_ratio",
                "aPTT",
                "aPTT_change",
                "aPTT_ratio",
                "Xa_inhibition",
            ]
        )
        return {}
