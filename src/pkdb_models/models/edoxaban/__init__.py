from pathlib import Path

EDOXABAN_PATH = Path(__file__).parent

MODEL_BASE_PATH = EDOXABAN_PATH / "models" / "results" / "models"
MODEL_PATH = MODEL_BASE_PATH / "edoxaban_body_flat.xml"

RESULTS_PATH = EDOXABAN_PATH / "results"
RESULTS_PATH_SIMULATION = RESULTS_PATH / "simulation"
RESULTS_PATH_FIT = RESULTS_PATH / "fit"

# DATA_PATH_BASE = EDOXABAN_PATH.parents[3] / "pkdb_data" / "studies"
DATA_PATH_BASE = EDOXABAN_PATH / "data"

DATA_PATH_EDOXABAN = DATA_PATH_BASE / "edoxaban"
DATA_PATHS = [
     DATA_PATH_EDOXABAN,
]
