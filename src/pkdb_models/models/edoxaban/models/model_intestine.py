"""Edoxaban intestine model."""

from pathlib import Path
import pandas as pd
from sbmlutils.factory import create_model, FactoryResult
from sbmlutils.converters import odefac
from sbmlutils import cytoscape as cyviz
import numpy as np
from sbmlutils.converters import odefac
from sbmlutils.cytoscape import visualize_sbml
from sbmlutils.factory import *
from sbmlutils.metadata import *

from pkdb_models.models.edoxaban.models import annotations
from pkdb_models.models.edoxaban.models import templates
from pkdb_models.models.edoxaban import MODEL_BASE_PATH

class U(templates.U):
    """UnitDefinitions"""

    per_hr = UnitDefinition("per_hr", "1/hr")
    mg_per_min = UnitDefinition("mg_per_min", "mg/min")


_m = Model(
    "edoxaban_intestine",
    name="Model for edoxaban absorption",
    notes="""
    # Model for edoxaban absorption
    """
    + templates.terms_of_use,
    creators=templates.creators,
    units=U,
    model_units=templates.model_units,
    annotations=[
        # tissue
        (BQB.OCCURS_IN, "fma/FMA:45615"),  # gut
        (BQB.OCCURS_IN, "bto/BTO:0000545"),  # gut
        (BQB.OCCURS_IN, "NCIT:C12736"),  # intestine
        (BQB.OCCURS_IN, "fma/FMA:7199"),  # intestine
        (BQB.OCCURS_IN, "bto/BTO:0000648"),  # intestine

        (BQB.HAS_PROPERTY, "NCIT:C79369"),  # Pharmacokinetics: Absorption
        (BQB.HAS_PROPERTY, "NCIT:C79372"),  # Pharmacokinetics: Excretion
    ] + templates.model_annotations
)

_m.compartments = [
    Compartment(
        "Vext",
        1.0,
        name="plasma",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        port=True,
        annotations=annotations.compartments["plasma"],
    ),
    Compartment(
        "Vlumen",
        value=1.0,
        name="lumen",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        constant=False,
        port=True,
        annotations=annotations.compartments["gu_lumen"],
    ),
    Compartment(
        "Vfeces",
        metaId="meta_Vfeces",
        value=1.0,
        unit=U.liter,
        constant=True,
        name="feces",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        port=True,
        annotations=annotations.compartments["feces"],
    ),
    Compartment(
        "Vstomach",
        metaId="meta_Vstomach",
        value=1,
        unit=U.liter,
        constant=True,
        name="stomach",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        port=True,
        annotations=annotations.compartments["stomach"],
    ),
    Compartment(
        "Vapical",
        np.nan,
        name="apical membrane (intestinal membrane enterocytes)",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.m2,
        spatialDimensions=2,
        annotations=annotations.compartments["apical"],
    ),
]

_m.species = [
    Species(
        f"edo_stomach",
        metaId=f"meta_edo_stomach",
        initialConcentration=0.0,
        compartment="Vstomach",
        substanceUnit=U.mmole,
        name=f"edoxaban (stomach)",
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["edo"],
        boundaryCondition=True,
    ),
    Species(
        "edo_ext",
        initialConcentration=0.0,
        name="edoxaban (plasma)",
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["edo"],
        port=True,
    ),
    Species(
        "edo_lumen",
        initialConcentration=0.0,
        name="edoxaban (intestinal volume)",
        compartment="Vlumen",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["edo"],
        port=True,
    ),
    Species(
        "edo_feces",
        initialConcentration=0.0,
        name="edoxaban (feces)",
        compartment="Vfeces",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["edo"],
        port=True,
    ),
    Species(
        "m4_lumen",
        initialConcentration=0.0,
        name="M4 (intestinal volume)",
        compartment="Vlumen",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["m4"],
        port=True,
    ),
    Species(
        "m6_lumen",
        initialConcentration=0.0,
        name="M6 (intestinal volume)",
        compartment="Vlumen",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["m6"],
        port=True,
    ),
    Species(
        "mx_lumen",
        initialConcentration=0.0,
        name="Mx (intestinal volume)",
        compartment="Vlumen",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["mx"],
        port=True,
    ),
    Species(
        "m4_feces",
        initialConcentration=0.0,
        name="M4 (feces)",
        compartment="Vfeces",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["m4"],
        port=True,
    ),
    Species(
        "m6_feces",
        initialConcentration=0.0,
        name="M6 (feces)",
        compartment="Vfeces",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["m6"],
        port=True,
    ),
    Species(
        "mx_feces",
        initialConcentration=0.0,
        name="Mx (feces)",
        compartment="Vfeces",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["mx"],
        port=True,
    )
]

_m.parameters = [
    Parameter(
        f"F_edo_abs",
        0.82,  # [0.82 - 0.88]
        U.dimensionless,
        constant=True,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        name=f"fraction absorbed edoxaban",
    ),
    Parameter(
        "EDOABS_k",
        0.05,
        unit=U.per_min,
        name="rate of edoxaban absorption",
        sboTerm=SBO.KINETIC_CONSTANT,
    ),
    Parameter(
        "f_absorption",
        1,
        unit=U.dimensionless,
        name="scaling factor for absorption rate",
        sboTerm=SBO.KINETIC_CONSTANT,
        notes="""1.0: normal absorption corresponding to tablet under fasting conditions.

        allows to change the velocity of absorption
        """
    ),
]

_m.rules.append(
    AssignmentRule(
        "absorption_edo",
        value="f_absorption * EDOABS_k * Vlumen * edo_lumen",
        unit=U.mmole_per_min,
        name="absorption edoxaban",
    ),
)

_m.reactions = [
    Reaction(
        "EDOABS",
        name="absorption edoxaban",
        equation="edo_lumen -> edo_ext",
        sboTerm=SBO.TRANSPORT_REACTION,
        compartment="Vapical",
        formula=("F_edo_abs * absorption_edo", U.mmole_per_min),
    ),

    Reaction(
        sid="EDOEXC",
        name=f"excretion edoxaban (feces)",
        compartment="Vlumen",
        equation=f"edo_lumen -> edo_feces",
        sboTerm=SBO.TRANSPORT_REACTION,
        formula=(
            f"(1 dimensionless - F_edo_abs) * absorption_edo",
            U.mmole_per_min,
        )
    ),
    Reaction(
        sid="M4EXC",
        name=f"excretion M4 (feces)",
        compartment="Vlumen",
        equation=f"m4_lumen -> m4_feces",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "M4EXC_k",
                0.10,
                unit=U.per_min,
                name="rate of M4 excretion",
                sboTerm=SBO.KINETIC_CONSTANT,
            ),
        ],
        formula=(
            f"f_absorption * M4EXC_k * Vlumen * m4_lumen",
            U.mmole_per_min,
        )
    ),
    Reaction(
            sid="M6EXC",
            name=f"excretion M6 (feces)",
            compartment="Vlumen",
            equation=f"m6_lumen -> m6_feces",
            sboTerm=SBO.TRANSPORT_REACTION,
            pars=[
                Parameter(
                    "M6EXC_k",
                    0.10,
                    unit=U.per_min,
                    name="rate of M6 excretion",
                    sboTerm=SBO.KINETIC_CONSTANT,
                ),
            ],
            formula=(
                f"f_absorption * M6EXC_k * Vlumen * m6_lumen",
                U.mmole_per_min,
            )
        ),
    Reaction(
            sid="MXEXC",
            name=f"excretion Mx (feces)",
            compartment="Vlumen",
            equation=f"mx_lumen -> mx_feces",
            sboTerm=SBO.TRANSPORT_REACTION,
            pars=[
                Parameter(
                    "MXEXC_k",
                    0.10,
                    unit=U.per_min,
                    name="rate of Mx excretion",
                    sboTerm=SBO.KINETIC_CONSTANT,
                ),
            ],
            formula=(
                f"f_absorption * MXEXC_k * Vlumen * mx_lumen",
                U.mmole_per_min,
            )
        )
]

_m.parameters.extend([
    Parameter(
        f"PODOSE_edo",
        0,
        U.mg,
        constant=False,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        name=f"oral dose edoxaban [mg]",
        port=True,
    ),
    Parameter(
        f"Ka_dis_edo",
        0.15,
        U.per_hr,
        constant=True,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        name=f"Ka_dis [1/hr] dissolution edoxaban",
        port=True
    ),
    Parameter(
        f"Mr_edo",
        548.058,
        U.g_per_mole,
        constant=True,
        name=f"Molecular weight edoxaban [g/mole]",
        sboTerm=SBO.MOLECULAR_MASS,
        port=True,
    ),
])

# -------------------------------------
# Dissolution of tablet/dose in stomach
# -------------------------------------
_m.reactions.extend(
    [
        # fraction dose available for absorption from stomach
        Reaction(
            sid=f"dissolution_edo",
            name=f"dissolution edoxaban",
            formula=(
                f"Ka_dis_edo/60 min_per_hr * PODOSE_edo/Mr_edo",
                U.mmole_per_min,
            ),
            equation=f"edo_stomach -> edo_lumen",
            compartment="Vlumen",
            notes="""Swallowing, dissolution of tablet, and transport into intestine.
            Overall process describing the rates of this processes.
            """
        ),
    ]
)
_m.rate_rules.append(
    RateRule(f"PODOSE_edo", f"-dissolution_edo * Mr_edo", U.mg_per_min),
)

model_intestine = _m


if __name__ == "__main__":
    results: FactoryResult = create_model(
        model=model_intestine,
        filepath=MODEL_BASE_PATH / f"{model_intestine.sid}.xml",
        sbml_level=3,
        sbml_version=2,
    )

    # ODE equations
    ode_factory = odefac.SBML2ODE.from_file(sbml_file=results.sbml_path)
    ode_factory.to_markdown(md_file=results.sbml_path.parent / f"{results.sbml_path.stem}.md")

    # visualization in Cytoscape
    cyviz.visualize_sbml(results.sbml_path, delete_session=False)