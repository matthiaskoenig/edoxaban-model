"""Kidney model for edoxaban."""
from pathlib import Path

from sbmlutils import cytoscape as cyviz
import numpy as np
from sbmlutils.console import console
from sbmlutils.converters import odefac
from sbmlutils.factory import *
from sbmlutils.metadata import *
from pkdb_models.models.edoxaban.models import annotations
from pkdb_models.models.edoxaban.models import templates
from pkdb_models.models.edoxaban import MODEL_BASE_PATH


class U(templates.U):
    """UnitDefinitions"""
    mg_per_g = UnitDefinition("mg_per_g", "mg/g")
    ml_per_l = UnitDefinition("ml_per_l", "ml/l")
    ml_per_min = UnitDefinition("ml_per_min", "ml/min")
    ml_per_min_bsa = UnitDefinition("ml_per_min_bsa", "ml/min/(1.73 * m^2)")


mid = "edoxaban_kidney"
version = 1

_m = Model(
    sid=mid,
    name="Model for renal edoxaban excretion.",
    notes=f"""
    Model for renal edoxaban excretion:
     - ~35 % in urine [Bathala2012]
     - ~62 % in feces [predominant- edoxaban (49.1%), other metabolites (M1, M4, M6) < 2%] [Bathala2012]
     
     The main elimination routes for edoxaban are through renal excretion (about 50%) 
     and hepatobiliary pathways (about 50%), with minimal overall metabolism
    """ + templates.terms_of_use,
    creators=templates.creators,
    units=U,
    model_units=templates.model_units,
    annotations=[
        # tissue
        (BQB.OCCURS_IN, "fma/FMA:7203"),  # kidney
        (BQB.OCCURS_IN, "bto/BTO:0000671"),
        (BQB.OCCURS_IN, "NCIT:C12415"),

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
        "Vki",
        value=0.3,
        unit=U.liter,
        name="kidney",
        constant=True,
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        annotations=annotations.compartments["ki"],
        port=True
    ),
    Compartment(
        "Vurine",
        1.0,
        name="urine",
        sboTerm=SBO.PHYSICAL_COMPARTMENT,
        unit=U.liter,
        port=True,
        annotations=annotations.compartments["urine"],
    ),
]

_m.species = [
    Species(
        "edo_ext",
        name="edoxaban (plasma)",
        initialConcentration=0.0,
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["edo"],
        port=True
    ),
    Species(
        "edo_urine",
        name="edoxaban (urine)",
        initialConcentration=0.0,
        compartment="Vurine",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["edo"],
        port=True
    ),
    Species(
        "m4_ext",
        name="M4 (plasma)",
        initialConcentration=0.0,
        compartment="Vext",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["m4"],
        port=True
    ),
    Species(
        "m4_urine",
        name="M4 (urine)",
        initialConcentration=0.0,
        compartment="Vurine",
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
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["m6"],
        port=True
    ),
    Species(
        "m6_urine",
        name="M6 (urine)",
        initialConcentration=0.0,
        compartment="Vurine",
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
        hasOnlySubstanceUnits=False,  # this is a concentration
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["mx"],
        port=True
    ),
    Species(
        "mx_urine",
        name="Mx (urine)",
        initialConcentration=0.0,
        compartment="Vurine",
        substanceUnit=U.mmole,
        hasOnlySubstanceUnits=True,
        sboTerm=SBO.SIMPLE_CHEMICAL,
        annotations=annotations.species["mx"],
        port=True
    )
]

# glomerular filtration rate, creatinine clearance & kidney function
_m.parameters.extend([
    Parameter(
        "BSA",
        1.73,
        U.m2,
        constant=False,
        name="body surface area [m^2]",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        port=True
    ),
    Parameter(
        "f_renal_function",
        name="parameter for renal function",
        value=1.0,
        unit=U.dimensionless,
        sboTerm=SBO.KINETIC_CONSTANT,
        notes="""scaling factor for renal function. 1.0: normal renal function; 
        <1.0: reduced renal function
        """
    ),
    Parameter(
        "egfr",
        np.nan,
        unit=U.ml_per_min_bsa,
        constant=False,
        name="estimated GFR [ml/min/(1.73 m^2)]",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        annotations=[
            (BQB.IS, "NCIT:C110935"),
        ],
        port=True,
        notes="""A laboratory test that estimates kidney function. It is calculated using an individual's 
        serum creatinine measurement, age, gender, and race."""
    ),
    Parameter(
        "crcl",
        np.nan,
        U.ml_per_min,
        constant=False,
        name="creatinine clearance [ml/min]",
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        annotations=[
            (BQB.IS, "CMO:0000765"),
        ],
        port=True,
        notes="""The clearance rate of creatinine, that is, the volume of plasma that is cleared of creatinine 
        by the kidneys per unit time. Creatinine clearance is calculated using the level of creatinine in a 
        sample of urine, usually one collected over a period of 24 hours, the corresponding plasma creatinine 
        level, and the volume of urine excreted. It is used as an approximation of the glomerular 
        filtration rate (GFR)."""
    ),
    Parameter(
        sid="egfr_healthy",
        name="estimated glomerular filtration (eGFR) rate healthy",
        value=100.0,
        unit=U.ml_per_min_bsa,
        sboTerm=SBO.QUANTITATIVE_SYSTEMS_DESCRIPTION_PARAMETER,
        annotations=[
        ],
        notes="""eGFR is often estimated via creatinine clearance, creatinine urinary amount,
        or creatinine plasma amounts.

        CKD-EPI Creatinine Equation (2021): This is the most recent and recommended equation for 
        estimating GFR35. It is more accurate than previous formulas and does not include 
        race as a variable.

        MDRD (Modification of Diet in Renal Disease) Study Equation
    """,
    )
])
_m.rules.extend([
    AssignmentRule(
        "egfr", "f_renal_function * egfr_healthy", unit=U.ml_per_min_bsa,
        name="estimated eGFR"
    ),
    AssignmentRule(
        "crcl", "egfr * BSA/(1.73 dimensionless) * 1.1 dimensionless", unit=U.ml_per_min,
        name="creatinine clearance",
        notes="CrCl typically overestimates GFR by 10-20% due to the active secretion"
              "of creatinine in the proximal tubules."
              "Using a 10 % overestimation in the model."
    ),
])

_m.reactions = [
    Reaction(
        sid="EDOEX",
        name="edoxaban excretion (EDOEX)",
        equation="edo_ext -> edo_urine",
        compartment="Vki",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "EDOEX_k",
                0.16,
                U.per_min,
                name="rate of edoxaban urinary excretion",
                sboTerm=SBO.KINETIC_CONSTANT,
            ),
        ],
        formula=(
            "f_renal_function * Vki * EDOEX_k * edo_ext"
        )
    ),
    Reaction(
        sid="M4EX",
        name="M4 excretion (M4EX)",
        equation="m4_ext -> m4_urine",
        compartment="Vki",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "M4EX_k",
                0.15,
                U.per_min,
                name="rate of M4 urinary excretion",
                sboTerm=SBO.KINETIC_CONSTANT,
            ),
        ],
        formula=(
            "f_renal_function * Vki * M4EX_k * m4_ext"
        )
    ),
    Reaction(
        sid="M6EX",
        name="M6 excretion (M6EX)",
        equation="m6_ext -> m6_urine",
        compartment="Vki",
        sboTerm=SBO.TRANSPORT_REACTION,
        pars=[
            Parameter(
                "M6EX_k",
                0.15,
                U.per_min,
                name="rate of M6 urinary excretion",
                sboTerm=SBO.KINETIC_CONSTANT,
            ),
        ],
        formula=(
            "f_renal_function * Vki * M6EX_k * m6_ext"
        )
    ),
    Reaction(
            sid="MXEX",
            name="Mx excretion (MXEX)",
            equation="mx_ext -> mx_urine",
            compartment="Vki",
            sboTerm=SBO.TRANSPORT_REACTION,
            pars=[
                Parameter(
                    "MXEX_k",
                    0.15,
                    U.per_min,
                    name="rate of Mx urinary excretion",
                    sboTerm=SBO.KINETIC_CONSTANT,
                ),
            ],
            formula=(
                "f_renal_function * Vki * MXEX_k * mx_ext"
            )
        )
]

model_kidney = _m


if __name__ == "__main__":
    results: FactoryResult = create_model(
        model=model_kidney,
        filepath=MODEL_BASE_PATH / f"{model_kidney.sid}.xml",
        sbml_level=3, sbml_version=2,
    )

    # ODE equations
    ode_factory = odefac.SBML2ODE.from_file(sbml_file=results.sbml_path)
    ode_factory.to_markdown(md_file=results.sbml_path.parent / f"{results.sbml_path.stem}.md")

    # visualization in Cytoscape
    cyviz.visualize_sbml(sbml_path=results.sbml_path, delete_session=False)


