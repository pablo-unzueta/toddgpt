import numpy as np
import pytest

from toddgpt.calculators.mace_calc import MaceCalculator
from toddgpt.tools.datatypes import AtomsDict


def mace_calculator():
    return MaceCalculator()


def test_mace_calculator_initialization():
    calc = MaceCalculator()
    assert calc.name == "mace_calculator"


@pytest.mark.parametrize(
    "run_type,expected_result",
    [
        ("sp_energy", {"sp_energy": -31.847604599}),
        ("forces", {"forces": [[0.0, 0.0, 0.36013065], [0.0, 0.0, -0.36013065]]}),
        (
            "minimize_positions",
            {
                "sp_energy": -31.851215287218437,
                "minimize_positions": np.array(
                    [[0.0, 0.0, 0.00995981], [0.0, 0.0, 0.73004019]]
                ),
            },
        ),
    ],
)
def test_mace_calculator_run(run_type, expected_result):
    mace_calculator = MaceCalculator()
    atoms_dict = AtomsDict(numbers=[1, 1], positions=[[0, 0, 0], [0, 0, 0.74]])

    result = mace_calculator._run(atoms_dict, run_type)
    print(result)

    assert isinstance(result, dict)

    if run_type == "sp_energy":
        assert "sp_energy" in result
        assert np.allclose(result["sp_energy"], expected_result["sp_energy"], atol=1e-7)
    elif run_type == "forces":
        assert "forces" in result
        assert np.allclose(result["forces"], expected_result["forces"], atol=1e-7)
    elif run_type == "minimize_positions":
        assert "sp_energy" in result and "minimize_positions" in result
        assert np.allclose(result["sp_energy"], expected_result["sp_energy"], atol=1e-7)
        assert np.allclose(result["minimize_positions"], expected_result["minimize_positions"], atol=1e-7)
    else:
        assert result == {}


def test_mace_calculator_with_invalid_atoms_dict():
    mace_calculator = MaceCalculator()
    invalid_atoms_dict = AtomsDict(numbers=[1], positions=[[0, 0, 0], [0, 0, 1]])

    with pytest.raises(ValueError):
        mace_calculator._run(invalid_atoms_dict, "sp_energy")
