import hashlib
import json

import numpy as np
import pytest

from grendel.models.alp import decay_matrix_elements as matrix_elements
from grendel.models.alp import exclusive_decays
from grendel.models.alp import model
from grendel.models.alp.decay_models import decay_data_dir, validate_decay_model
from grendel.models.alp.templates_pythia import (
    _masses_from_csv,
    _validate_resumed_template,
)


def _sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_unknown_decay_model_is_rejected():
    with pytest.raises(ValueError, match="unknown decay model"):
        validate_decay_model("approximate")


def test_template_resume_is_model_count_and_surrogate_specific(tmp_path):
    path = tmp_path / "template.npz"
    np.savez_compressed(
        path,
        decay_model=np.array("2501"),
        n_templates=np.int32(20_000),
        gluon_surrogate=np.array("uds"),
    )
    _validate_resumed_template(path, 20_000, "uds", "2501")
    with pytest.raises(RuntimeError, match="model/count/surrogate"):
        _validate_resumed_template(path, 10_000, "uds", "2501")
    with pytest.raises(RuntimeError, match="model/count/surrogate"):
        _validate_resumed_template(path, 20_000, "gg", "2501")


def test_template_generator_accepts_exact_mass_grid_csv(tmp_path):
    path = tmp_path / "mass_grid.csv"
    path.write_text("mass_GeV\n1.18\n1.19\n1.20\n")
    assert _masses_from_csv(path) == [1.18, 1.19, 1.20]

    path.write_text("mass_GeV\n1.18\n1.20\n1.19\n")
    with pytest.raises(ValueError, match="strictly increasing"):
        _masses_from_csv(path)
