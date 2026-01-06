"""Parameter fit problems"""
from typing import Dict, List
from sbmlsim.fit.helpers import f_fitexp, filter_empty
from sbmlutils.console import console
from sbmlutils.log import get_logger

from sbmlsim.fit import FitExperiment, FitMapping

from pkdb_models.models.edoxaban import EDOXABAN_PATH, DATA_PATHS
from pkdb_models.models.edoxaban.experiments.metadata import (
    Tissue, Route, Dosing, ApplicationForm, Health,
    Fasting, EdoxabanMappingMetaData, Coadministration
)
from pkdb_models.models.edoxaban.experiments.studies import *


logger = get_logger(__name__)


# --- Filters ---
def filter_baseline(fit_mapping_key: str, fit_mapping: FitMapping) -> bool:
    """Return baseline experiments/mappings for reference model."""

    metadata: EdoxabanMappingMetaData = fit_mapping.metadata

    # only PO and IV (no SL, MU, RE)
    if metadata.route not in {Route.PO, Route.IV}:
        return False

    # filter coadminstration
    if metadata.coadministration != Coadministration.NONE:
        return False

    # filter health (no renal, cardiac impairment, ...)
    if metadata.health not in {Health.HEALTHY, Health.T2DM, Health.HYPERTENSION}:
        return False

    # filter multiple dosing (only single dosing)
    if metadata.dosing == Dosing.MULTIPLE:
        return False

    # # only fasted subjects
    # if metadata.fasting not in {Fasting.FASTED, Fasting.NR}:
    #     return False

    # remove outliers
    if metadata.outlier is True:
        return False

    return True


# --- Fit experiments ---
f_fitexp_kwargs = dict(
    experiment_classes  = [
        # all
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
    base_path=EDOXABAN_PATH,
    data_path=DATA_PATHS,
)
# --- Experiment classes ---

def filter_pharmacokinetics(fit_mapping_key: str, fit_mapping: FitMapping) -> bool:
    """Only pk data."""
    yid = "__".join(fit_mapping.observable.y.sid.split("__")[1:])
    if yid not in {
        "Cve_edo",
        "Cve_m4",
        "Cve_m6",
        "Cve_mx",
        "Cve_edo_tot",
        "Aurine_edo",
        "Aurine_m4",
        "Aurine_m6",
        "Aurine_mx",
        "Aurine_edo_tot",
        "Afeces_edo",
        "Afeces_m4",
        "Afeces_m6",
        "Afeces_mx",
        "Afeces_edo_tot",
    }:
        return False
    return True

def filter_pharmacodynamics(fit_mapping_key: str, fit_mapping: FitMapping) -> bool:
    """Only pharmacodynamics data."""
    yid = "__".join(fit_mapping.observable.y.sid.split("__")[1:])
    if yid not in {
        "PT",
        "PT_change",
        "PT_ratio",
        "aPTT",
        "aPTT_change",
        "aPTT_ratio",
        "Xa_inhibition",
    }:
        return False

    return True

def f_fitexp_all():
    """All data."""
    return f_fitexp(metadata_filters=filter_empty, **f_fitexp_kwargs)

def f_fitexp_control() -> Dict[str, List[FitExperiment]]:
    """Control data."""
    return f_fitexp(metadata_filters=[filter_baseline], **f_fitexp_kwargs)

def f_fitexp_pharmacokinetics() -> Dict[str, List[FitExperiment]]:
    """Pharmacodynamics data."""
    return f_fitexp(metadata_filters=[filter_baseline, filter_pharmacokinetics], **f_fitexp_kwargs)

def f_fitexp_pharmacodynamics() -> Dict[str, List[FitExperiment]]:
    """Pharmacodynamics data."""
    return f_fitexp(metadata_filters=[filter_baseline, filter_pharmacodynamics], **f_fitexp_kwargs)


if __name__ == "__main__":
    """Test construction of FitExperiments."""

    for f in [
        f_fitexp_all,
        f_fitexp_control,
        f_fitexp_pharmacokinetics,
        f_fitexp_pharmacodynamics,
    ]:
        console.rule(style="white")
        console.print(f"{f.__name__}")
        fitexp = f()
