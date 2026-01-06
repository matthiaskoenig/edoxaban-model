[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13987865.svg)](https://doi.org/10.5281/zenodo.13987865)
[![GitHub Action](https://github.com/matthiaskoenig/edoxaban-model/actions/workflows/python.yml/badge.svg)](https://github.com/matthiaskoenig/edoxaban-model/actions/workflows/python.yml)
[![GitHub Action](https://github.com/matthiaskoenig/edoxaban-model/actions/workflows/docker.yml/badge.svg)](https://github.com/matthiaskoenig/edoxaban-model/actions/workflows/docker.yml)

# edoxaban model
This repository provides the edoxaban physiologically based pharmacokinetics/ pharmacodynamics (PBPK/PD) model.

The model is distributed as [SBML](http://sbml.org) available from [`edoxaban_body_flat.xml`](./models/edoxaban_body_flat.xml) with 
corresponding SBML4humans model report at [https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/edoxaban-model/main/models/edoxaban_body_flat.xml](https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/edoxaban-model/main/models/edoxaban_body_flat.xml) and equations from [`edoxaban_body_flat.md`](./models/edoxaban_body_flat.md).

The COMBINE archive is available from [`edoxaban_model.omex`](./edoxaban_model.omex).

![model overview](./figures/edoxaban_model.png)

### Comp submodels
The liver submodel is available from [`edoxaban_liver.xml`](./models/edoxaban_liver.xml) with corresponding SBML4humans report at
[https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/edoxaban-model/main/models/edoxaban_liver.xml](https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/edoxaban-model/main/models/edoxaban_liver.xml) and equations from [`edoxaban_liver.md`](./models/edoxaban_liver.md).

The kidney submodel is available from [`edoxaban_kidney.xml`](./models/edoxaban_kidney.xml) with corresponding SBML4humans report at
[https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/edoxaban-model/main/models/edoxaban_kidney.xml](https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/edoxaban-model/main/models/edoxaban_kidney.xml) and equations from [`edoxaban_kidney.md`](./models/edoxaban_kidney.md).

The intestine submodel is available from [`edoxaban_intestine.xml`](./models/edoxaban_intestine.xml) with corresponding SBML4humans report at
[https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/edoxaban-model/main/models/edoxaban_intestine.xml](https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/edoxaban-model/main/models/edoxaban_intestine.xml) and equations from [`edoxaban_intestine.md`](./models/edoxaban_intestine.md).

The whole-body submodel is available from [`edoxaban_body.xml`](./models/edoxaban_body.xml) with corresponding SBML4humans report at
[https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/edoxaban-model/main/models/edoxaban_body.xml](https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/edoxaban-model/main/models/edoxaban_body.xml) and equations from [`edoxaban_body.md`](./models/edoxaban_body.md).

The coagulation submodel is available from [`edoxaban_coagulation.xml`](./models/edoxaban_body.xml) with corresponding SBML4humans report at
[https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/edoxaban-model/main/models/edoxaban_coagulation.xml](https://sbml4humans.de/model_url?url=https://raw.githubusercontent.com/matthiaskoenig/edoxaban-model/main/models/edoxaban_body.xml) and equations from [`edoxaban_coagulation.md`](./models/edoxaban_body.md).


## How to cite
To cite the model repository

> Babaeva, M., Myshkina, M. & König, M. (2025).
> *Physiologically based pharmacokinetic/ pharmacodynamic (PBPK/PD) model of edoxaban.*   
> Zenodo. [https://doi.org/10.5281/zenodo.13987865](https://doi.org/10.5281/zenodo.13987865)

## License

* Source Code: [MIT](https://opensource.org/license/MIT)
* Documentation: [CC BY-SA 4.0](http://creativecommons.org/licenses/by-sa/4.0/)
* Models: [CC BY-SA 4.0](http://creativecommons.org/licenses/by-sa/4.0/)

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.

## Run simulations
### python
Clone the repository 
```bash
git clone https://github.com/matthiaskoenig/edoxaban-model.git
cd edoxaban-model
```

#### uv
Setup environment with uv (https://docs.astral.sh/uv/getting-started/installation/)
```bash
uv sync
```
Run the complete analysis:
```bash
uv run run_edoxaban -a all -r results
```

#### pip
If you use pip install the package via
```bash
pip install -e .
```
Run the complete analysis in the environment via:
```bash
run run_edoxaban -a all -r results
```

### docker
Simulations can also be run within a docker container:

```bash
docker run -v "${PWD}/results:/results" -it matthiaskoenig/edoxaban:latest /bin/bash
```

Run the complete analysis:
```bash
uv run run_edoxaban -a all -r /results
```
The results are written into the mounted `/results` folder on the host.

In case of permission issues with the mounted folder, adjust ownership and access rights with:
```bash
sudo chown $(id -u):$(id -g) -R "${PWD}/results"
sudo chmod 775 "${PWD}/results"
```

## Funding

Matthias König was supported by the Federal Ministry of Education and Research (BMBF, Germany) within LiSyM by grant number 031L0054 and ATLAS by grant number 031L0304B and by the German Research Foundation (DFG) within the Research Unit Program FOR 5151 QuaLiPerF (Quantifying Liver Perfusion-Function Relationship in Complex Resection - A Systems Medicine Approach) by grant number 436883643 and by grant number 465194077 (Priority Programme SPP 2311, Subproject SimLivA). This work was supported by the BMBF-funded de.NBI Cloud within the German Network for Bioinformatics Infrastructure (de.NBI) (031A537B, 031A533A, 031A538A, 031A533B, 031A535A, 031A537C, 031A534A, 031A532B). 

© 2025-2026 Mariia Babaeva, Mariia Myshkina & Matthias König, [Systems Medicine of the Liver](https://livermetabolism.com)
