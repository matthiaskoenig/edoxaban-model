"""Liver model for edoxaban."""

from sbmlutils import cytoscape as cyviz
from pathlib import Path
import pandas as pd
import numpy as np
from sbmlutils.converters import odefac
from sbmlutils.factory import *
from sbmlutils.metadata import *
from pkdb_models.models.edoxaban.models import annotations
from pkdb_models.models.edoxaban.models import templates
from pkdb_models.models.edoxaban import MODEL_BASE_PATH


class U(templates.U):
    """UnitDefinitions"""

    pass


mid = "edoxaban_liver"
version = 1

_m = Model(
    sid=mid,
    name="Model for hepatic edoxaban metabolism.",
    notes=f"""
    Model for hepatic edoxaban metabolism.
    
    """ + templates.terms_of_use,
    creators=templates.creators,
    units=U,
    model_units=templates.model_units,
    annotations=[
        # tissue
        (BQB.OCCURS_IN, "fma/FMA:7197"),
        (BQB.OCCURS_IN, "bto/BTO:0000759"),
        (BQB.OCCURS_IN, "NCIT:C12392"),

        (BQB.HAS_PROPERTY, "NCIT:C79371"),  # Pharmacokinetics: Metabolism
        (BQB.HAS_PROPERTY, "NCIT:C79372"),  # Pharmacokinetics: Excretion
    ] + templates.model_annotations
)

_m.compartments = [
    Compartment(
        "Vext",
        value=1.5,
        unit=U.liter,
        name="plasma",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["plasma"],
        port=True
    ),
    Compartment(
        "Vli",
        value=1.5,
        unit=U.liter,
        name="liver",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["li"],
        port=True
    ),
    Compartment(
        "Vmem",
        value=np.nan,
        unit=U.m2,
        name="membrane",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["basolateral"],
        spatialDimensions=2,
        port=True
    ),
    Compartment(
        "Vbi",
        1.0,
        name="bile",
        unit=U.liter,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["bi"],
        port=True,
    ),
    Compartment(
        "Vlumen",
        1.0,
        name="gut lumen",
        unit=U.liter,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["gu"],
        port=True,
    ),
]

_m.species = [
    Species(
        "edo_ext",
        name="edoxaban (plasma)",
        initialConcentration=0.0,
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["edo"],
        port=True
    ),
    Species(
        "edo",
        name="edoxaban (liver)",
        initialConcentration=0.0,
        compartment="Vli",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["edo"]
    ),
    Species(
        "m4_ext",
        name="M4 (plasma)",
        initialConcentration=0.0,
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["m4"],
        port=True
    ),
    Species(
        "m4",
        name="M4 (liver)",
        initialConcentration=0.0,
        compartment="Vli",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["m4"]
    ),
    Species(
        "m4_bi",
        initialConcentration=0.0,
        name="M4 (bile)",
        compartment="Vbi",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["m4"]
    ),
    Species(
        "m4_lumen",
        initialConcentration=0.0,
        name="M4 (gut)",
        compartment="Vlumen",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["m4"],
        port=True
    ),
    Species(
        "m6_ext",
        name="M6 (plasma)",
        initialConcentration=0.0,
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["m6"],
        port=True
    ),
    Species(
        "m6",
        name="M6 (liver)",
        initialConcentration=0.0,
        compartment="Vli",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["m6"],
    ),
    Species(
        "m6_bi",
        initialConcentration=0.0,
        name="M6 (bile)",
        compartment="Vbi",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["m6"],
    ),
    Species(
        "m6_lumen",
        initialConcentration=0.0,
        name="M6 (gut)",
        compartment="Vlumen",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["m6"],
        port=True
    ),
    Species(
        "mx_ext",
        name="Mx (plasma)",
        initialConcentration=0.0,
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["mx"],
        port=True
    ),
    Species(
        "mx",
        name="Mx (liver)",
        initialConcentration=0.0,
        compartment="Vli",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["mx"],
    ),
    Species(
        "mx_bi",
        initialConcentration=0.0,
        name="Mx (bile)",
        compartment="Vbi",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["mx"],
    ),
    Species(
        "mx_lumen",
        initialConcentration=0.0,
        name="Mx (gut)",
        compartment="Vlumen",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["mx"],
        port=True
    )
]

_m.reactions = [
    Reaction(
        sid="EDOIM",
        name="edoxaban import (EDOIM)",
        equation="edo_ext <-> edo",

        compartment="Vmem",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "EDOIM_Vmax",
                1000.0,
                U.per_min,
                name="Vmax edoxaban import",
                sboTerm=SBO.MAXIMAL_VELOCITY,
            ),
        ],
        formula=(
            "EDOIM_Vmax * Vli * (edo_ext - edo)"
        ),
        notes="""Assuming fast equilibration."""
    ),
    Reaction(
        sid="EDO2M4",
        name="edoxaban -> M4 (CES1)",
        equation="edo-> m4",
        compartment="Vli",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        pars=[
            Parameter(
                "EDO2M4_Vmax",
                0.1,
                U.per_min,
                name="Vmax edo-> m4 conversion",
                sboTerm=SBO.MAXIMAL_VELOCITY,
            ),
        ],
        formula=(
            "EDO2M4_Vmax * Vli * edo"
        ),
        notes="""Main metabolite M4."""
    ),
    Reaction(
        sid="EDO2M6",
        name="edoxaban -> M6 (CYP3A4/5)",
        equation="edo -> m6",
        compartment="Vli",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        pars=[
            Parameter(
                "EDO2M6_f",
                55.1/145.5,
                U.dimensionless,
                name="Vmax scaling edo -> m6 conversion",
                sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                notes="""
                Setting relative rates based on abundance of the metabolites:
                
                f = AUC_M6/AUC_M4 = 55.1/145.5
                """
            ),
        ],
        formula=(
            "EDO2M6_f * EDO2M4_Vmax * Vli * edo"
        ),
    ),
    Reaction(
        sid="EDO2MX",
        name="edoxaban -> Mx",
        equation="edo -> mx",
        compartment="Vli",
        sboTerm=SBO.BIOCHEMICAL_REACTION,
        pars=[
            Parameter(
                "EDO2MX_f",
                (97.6 - 72.8 - 3.49 - 2.11)/3.49,
                U.dimensionless,
                name="Vmax scaling edo -> mx conversion",
                sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
                notes="""
                Setting relative rates based on recovery of metabolites
                
                f_mx = (97.6 - 72.8 - 3.49 - 2.11)/3.49 = 5.50                
                """
            ),
        ],
        formula=(
            "EDO2MX_f * EDO2M4_Vmax * Vli * edo"
        ),
    ),
    Reaction(
        sid="M4EX",
        name="M4 export (M4EX)",
        equation="m4 <-> m4_ext",
        compartment="Vmem",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "M4EX_Vmax",
                1000.0,
                U.per_min,
                name="Vmax M4 export",
                sboTerm=SBO.MAXIMAL_VELOCITY,
            ),
        ],
        formula=(
            "M4EX_Vmax * Vli * (m4 - m4_ext)"
        ),
        notes="""Assuming fast equilibration."""
    ),
    Reaction(
        sid="M6EX",
        name="M6 export (M6EX)",
        equation="m6 <-> m6_ext",
        compartment="Vmem",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "M6EX_Vmax",
                1000.0,
                U.per_min,
                name="Vmax M6 export",
                sboTerm=SBO.MAXIMAL_VELOCITY,
            ),
        ],
        formula=(
            "M6EX_Vmax * Vli * (m6 - m6_ext)"
        ),
        notes="""Assuming fast equilibration."""
    ),
    Reaction(
            sid="MXEX",
            name="Mx export (MXEX)",
            equation="mx <-> mx_ext",
            compartment="Vmem",
            sboTerm=SBO.TRANSPORT_REACTION,
            pars=[
                Parameter(
                    "MXEX_Vmax",
                    1000.0,
                    U.per_min,
                    name="Vmax Mx export",
                    sboTerm=SBO.MAXIMAL_VELOCITY,
                ),
            ],
            formula=(
                "MXEX_Vmax * Vli * (mx - mx_ext)"
            ),
            notes="""Assuming fast equilibration."""
    ),
]

# biliary excretion
_m.parameters.append(
    Parameter(
       "MXEXBI_k",
       9.99999999999998e-05,
        unit=U.per_min,
        name="rate for edoxaban metabolites export in bile",
        sboTerm=SBO.KINETIC_CONSTANT,
        ),
)

for sid in ["m4", "m6", "mx"]:
    _m.reactions.extend([
        Reaction(
            f"{sid.upper()}EXBI",
            name=f"{sid.upper()} bile export",
            equation=f"{sid} -> {sid}_bi",
            sboTerm=SBO.TRANSPORT_REACTION,
            compartment="Vmem",
            formula=(
                f"MXEXBI_k * Vli * {sid}",
                U.mmole_per_min,
            ),
        ),
        Reaction(
            f"{sid.upper()}EXEHC",
            name=f"{sid.upper()} enterohepatic circulation",
            equation=f"{sid}_bi -> {sid}_lumen",
            sboTerm=SBO.TRANSPORT_REACTION,
            compartment="Vlumen",
            formula=(
                f"{sid.upper()}EXBI",
                U.mmole_per_min,
            ),
        ),
    ])

model_liver = _m


if __name__ == "__main__":
    results: FactoryResult = create_model(
        model=model_liver,
        filepath=MODEL_BASE_PATH / f"{model_liver.sid}.xml",
        sbml_level=3,
        sbml_version=2,
    )

    # ODE equations
    ode_factory = odefac.SBML2ODE.from_file(sbml_file=results.sbml_path)
    ode_factory.to_markdown(md_file=results.sbml_path.parent / f"{results.sbml_path.stem}.md")

    # visualization in Cytoscape
    cyviz.visualize_sbml(results.sbml_path, delete_session=True)
