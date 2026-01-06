"""Template definition."""
from datetime import datetime
from pymetadata.identifiers.miriam import BQB
from sbmlutils.factory import *


class U(Units):
    """UnitDefinitions."""

    mmole = UnitDefinition("mmole")
    min = UnitDefinition("min")
    mg = UnitDefinition("mg")
    m2 = UnitDefinition("m2", "meter^2")
    mM = UnitDefinition("mM", "mmole/liter")
    mmole_per_min = UnitDefinition("mmole_per_min", "mmole/min")
    mmole_per_min_l = UnitDefinition("mmole_per_min_l", "mmole/min/l")
    l_per_min = UnitDefinition("l_per_min", "l/min")
    per_min = UnitDefinition("per_min", "1/min")
    per_mmole = UnitDefinition("per_mmole", "1/mmole")
    mg_per_s_l = UnitDefinition("mg_per_s_l", "mg/s/l")
    g_per_mole = UnitDefinition("g_per_mole", "g/mole")
    min_per_hr = UnitDefinition("min_per_hr", "min/hour")


model_units = ModelUnits(
    time=U.min,
    extent=U.mmole,
    substance=U.mmole,
    length=U.meter,
    area=U.m2,
    volume=U.liter,
)

model_annotations = [
    # taxonomy
    (BQB.HAS_TAXON, "taxonomy/9606"),  # human
    (BQB.HAS_TAXON, "snomedct/337915000"),  # human

    # modelling approach
    (BQB.HAS_PROPERTY, "mamo/MAMO_0000046"),  # ordinary differential equation model

    # biological process explained by model (GO/NCIT)
    (BQB.HAS_PROPERTY, "NCIT:C96539"),  # edoxaban
    (BQB.HAS_PROPERTY, "CHEBI:85973"),  # edoxaban
    (BQB.HAS_PROPERTY, "pubchem.compound/10280735"),  # edoxaban
    (BQB.HAS_PROPERTY, "inchikey/HGVDHZBSSITLCT-JLJPHGGASA-N"),  # edoxaban
    (BQB.HAS_PROPERTY, "snomedct/712778008"),  # edoxaban
    (BQB.HAS_PROPERTY, "NCIT:C180619"),  # Direct Oral Anticoagulant

    (BQB.HAS_PROPERTY, "NCIT:C15299"),  # pharmacokinetics
    (BQB.HAS_PROPERTY, "NCIT:C15720"),  # pharmacodynamics

    # model relevance to a particular area
    (BQB.HAS_PROPERTY, "go/GO:0007596"),  # blood coagulation

    # model reference
    # (BQB.IS_DESCRIBED_BY, "doi/10.20944/preprints202507.1945.v1"),  # preprint
    (BQB.IS_DESCRIBED_BY, "doi/10.5281/zenodo.18165479"),
]

creators = [
    Creator(
        familyName="Babaeva",
        givenName="Mariia",
        email="mariia.babaeva@charite.de",
        organization="Charité University Medicine Berlin, Master Molecular Medicine",
    ),
    Creator(
       familyName="König",
       givenName="Matthias",
       email="koenigmx@hu-berlin.de",
       organization="Humboldt-University Berlin, Institute for Theoretical Biology",
       site="https://livermetabolism.com",
    ),
]

terms_of_use = """
    The content of this model has been carefully created in a manual research effort.

    ## Terms of use
    Copyright © {year} <a href="{site}" title="{given_name} {family_name}" target="_blank">{given_name} {family_name}</a>.

    <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">
    <img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

    Redistribution and use of any part of this model, with or without modification, are permitted provided
    that the following conditions are met:

    - Redistributions of this SBML file must retain the above copyright notice, this list of conditions and the following disclaimer.
    - Redistributions in a different form must reproduce the above copyright notice, this list of conditions
      and the following disclaimer in the documentation and/or other materials provided
      with the distribution.

    This model is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
    implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
""".format(
    year=datetime.now().year,
    given_name=creators[-1].givenName,
    family_name=creators[-1].familyName,
    site=creators[-1].site,
)