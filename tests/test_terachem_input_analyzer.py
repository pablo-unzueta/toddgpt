import pytest

from toddgpt.terachem_input_analyzer import SuggestTerachemRunType


@pytest.fixture
def suggest_terachem_run_type():
    return SuggestTerachemRunType()

@pytest.mark.parametrize("description, expected_output", [
    ("Perform a single point energy calculation", "energy"),
    ("Perform a singel-point energy job", "energy"),  # noqa: E501
    ("Perform a single point eneryg calculation", "energy"),  # noqa: E501
    ("Perform a eneryg calculation", "energy"),  # noqa: E501
    ("Run a geometry optimization", "minimize"),
    ("Run a geometry optimisation", "minimize"),  # noqa: E501
    ("Rnu a geometry optimization", "minimize"),  # noqa: E501
    ("Run a geometry optimization", "minimize"),  # noqa: E501
    ("Conduct a molecular dynamics simulation", "md"),
    ("Condcut a molecualr dynmaics simulaton", "md"),  # noqa: E501
    ("Run a MD simulation", "md"),
    ("Find the minimal energy conical intersection", "conical"),
    ("Find the minimal energy conical intersection", "conical"),  # noqa: E501
    ("Locate MECI", "conical"),
    ("Locate a transition state", "ts"),
    ("Locat a transition state", "ts"),  # noqa: E501
    ("Find TS", "ts"),
    ("Perform a nudged elastic band calculation", "neb"),
    ("Perform a nudgd elastc band calculaton", "neb"),  # noqa: E501
    ("Run NEB", "neb"),
    ("Calculate vibrational frequencies", "frequencies"),
    ("Calculate vibrational frequencies", "frequencies"),  # noqa: E501
    ("Compute Hessian", "frequencies"),
    ("Generate initial conditions", "initcond"),
    ("Generate initial conditions", "initcond"),  # noqa: E501
    ("Create starting conditions", "initcond"),
    ("Compute coupling elements", "coupling"),
    ("Compue couplig elements", "coupling"),  # noqa: E501
    ("Calculate nonadiabatic coupling", "coupling"),
    ("Calculate forces and gradients", "gradient"),
    ("Calculate forcs and gradents", "gradient"),  # noqa: E501
])
def test_suggest_terachem_run_type(suggest_terachem_run_type, description, expected_output):
    assert suggest_terachem_run_type.suggest_input(description) == expected_output

def test_multiple_matches(suggest_terachem_run_type):
    description = "Perform an energy calculation and optimize the geometry"
    result = suggest_terachem_run_type.suggest_input(description)
    assert result in ["energy", "minimize"]
    
    # Test that the result is the one with the highest match probability
    energy_result = suggest_terachem_run_type.suggest_input("Perform an energy calculation")
    minimize_result = suggest_terachem_run_type.suggest_input("Optimize the geometry")
    assert result == (energy_result if energy_result == result else minimize_result)

def test_unknown_calculation_type(suggest_terachem_run_type):
    description = "Perform a quantum teleportation experiment"
    result = suggest_terachem_run_type.suggest_input(description)
    assert "Unknown calculation type" in result
    assert "Please provide more specific information" in result

def test_case_insensitivity(suggest_terachem_run_type):
    assert suggest_terachem_run_type.suggest_input("SINGLE POINT ENERGY") == "energy"
    assert suggest_terachem_run_type.suggest_input("MOLECULAR DYNAMICS") == "md"
    assert suggest_terachem_run_type.suggest_input("Gradient Calculation") == "gradient"

def test_hyphen_underscore_handling(suggest_terachem_run_type):
    assert suggest_terachem_run_type.suggest_input("single-point energy") == "energy"
    assert suggest_terachem_run_type.suggest_input("single_point_energy") == "energy"
    assert suggest_terachem_run_type.suggest_input("nudged-elastic-band") == "neb"
    assert suggest_terachem_run_type.suggest_input("transition_state_search") == "ts"

def test_spell_check(suggest_terachem_run_type):
    assert suggest_terachem_run_type.suggest_input("singel point energy") == "energy"
    assert suggest_terachem_run_type.suggest_input("moleclar dynamics") == "md"
    assert suggest_terachem_run_type.suggest_input("geometry optimization") == "minimize"

def test_multiple_keywords(suggest_terachem_run_type):
    result = suggest_terachem_run_type.suggest_input("Calculate energy and optimize geometry")
    assert result in ["energy", "minimize"]
