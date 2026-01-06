"""Model for coagulation factors and readouts."""
from dataclasses import dataclass

import numpy as np

from sbmlutils import cytoscape as cyviz
from sbmlutils.converters import odefac
from sbmlutils.factory import *
from sbmlutils.metadata import *

from pkdb_models.models.edoxaban.models import annotations
from pkdb_models.models.edoxaban.models import templates


class U(templates.U):
    """UnitDefinitions"""
    pass


mid = "edoxaban_coagulation"
version = 1

_m = Model(
    sid=mid,
    name="Model for coagulation factors and readouts.",
    notes=f"""
    Model for coagulation factors and readouts.        
    """ + templates.terms_of_use,
    creators=templates.creators,
    units=U,
    model_units=templates.model_units,
    annotations=[
        (BQB.HAS_PROPERTY, "go/GO:0007596"),  # blood coagulation
    ] + templates.model_annotations
)

_m.compartments = [
    Compartment(
        "Vplasma",
        value=5.0,
        unit=U.liter,
        name="plasma",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["plasma"],
        port=True
    ),
]

_m.species = [
    Species(
        "Cve_edo",
        name="edoxaban",
        initialConcentration=0,
        compartment="Vplasma",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["edo"],
    )
]

names = {
    "PT": "prothrombin time",
    "aPTT": "activated partial thromboplastin",
    "Xa_inhibition": "inhibition of factor Xa",
}
reference_values = {
    "PT": 12.5,  # [s]
    "aPTT": 28.4,  # [s]
    "Xa_inhibition": 0, # [-] not inhibited
}
unit_values = {
    "PT": U.second,
    "aPTT": U.second,
    "Xa_inhibition": U.dimensionless,
}

# Xa, PT aPTT
for sid in ["PT", "aPTT"]:
    _m.parameters.extend([
        Parameter(
            sid=f"{sid}_ref",
            name=f"{names[sid]} reference",
            value=reference_values[sid],
            unit=unit_values[sid],
            constant=True,
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        ),
    ])

for sid in ["PT", "aPTT", "Xa_inhibition"]:
    _m.parameters.extend([
        Parameter(
            sid=sid,
            name=names[sid],
            value=reference_values[sid],
            unit=unit_values[sid],
            constant=False,
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        )
    ])

for sid in ["PT", "aPTT"]:
    # absolute change
    _m.parameters.append(
        Parameter(
            f"{sid}_change",
            name=f"{names[sid]} change",
            value=np.nan,
            unit=U.second,
            constant=False,
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes=f"Absolute change to baseline {names[sid]}",
        )
    )
    _m.rules.append(
        AssignmentRule(
            f"{sid}_change", f"{sid}-{sid}_ref", unit=U.second
        )
    )
    # ratio to baseline
    _m.parameters.append(
        Parameter(
            f"{sid}_ratio",
            name=f"{names[sid]} ratio",
            value=np.nan,
            unit=U.dimensionless,
            constant=False,
            sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
            notes=f"Ratio relative to baseline {names[sid]}",
        )
    )
    _m.rules.append(
        AssignmentRule(
            f"{sid}_ratio", f"{sid}/{sid}_ref", unit=U.dimensionless
        )
    )


# Effect of edoxaban on PT, aPTT and Xa_inhibition
_m.parameters.extend([
    Parameter(
        sid="Emax_PT",
        name="Effect constant for PT inhibition",
        value=3.55309,
        unit=U.dimensionless,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
    ),
    Parameter(
        sid="EC50_edo_PT",
        name="Effect constant for PT",
        value=0.00358,
        unit=U.mM,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
    ),
    Parameter(
        sid="Emax_aPTT",
        name="Effect constant for aPTT inhibition",
        value=0.950019,
        unit=U.dimensionless,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
    ),
    Parameter(
        sid="EC50_edo_aPTT",
        name="Effect constant for aPTT",
        value=0.0004089,
        unit=U.mM,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
    ),
    Parameter(
        sid="Emax_Xa",
        name="Effect constant for Xa inhibition",
        value=0.686093,
        unit=U.dimensionless,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
    ),
    Parameter(
        sid="EC50_edo_Xa",
        name="Effect constant for Xa inhibition",
        value=0.0002959,
        unit=U.mM,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
    ),
])

_m.rules.extend([
    AssignmentRule(
        "PT", "PT_ref * (1 dimensionless + Emax_PT * Cve_edo / (Cve_edo + EC50_edo_PT))", unit=U.second
    ),
    AssignmentRule(
        "aPTT", "aPTT_ref * (1 dimensionless + Emax_aPTT * Cve_edo / (Cve_edo + EC50_edo_aPTT))", unit=U.second
    ),
    AssignmentRule(
        "Xa_inhibition", "Emax_Xa * Cve_edo / (Cve_edo + EC50_edo_Xa)", unit=U.dimensionless
    )
])

model_coagulation = _m


if __name__ == "__main__":
    from pkdb_models.models.edoxaban import MODEL_BASE_PATH

    results: FactoryResult = create_model(
        model=model_coagulation,
        filepath=MODEL_BASE_PATH / f"{model_coagulation.sid}.xml",
        sbml_level=3, sbml_version=2,
        validation_options=ValidationOptions(units_consistency=True)
    )
    # create differential equations
    md_path = MODEL_BASE_PATH / f"{model_coagulation.sid}.md"
    ode_factory = odefac.SBML2ODE.from_file(sbml_file=results.sbml_path)
    ode_factory.to_markdown(md_file=md_path)

    cyviz.visualize_sbml(results.sbml_path, delete_session=True)
