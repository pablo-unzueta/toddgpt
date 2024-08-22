import pytest
from ase import Atoms

from toddgpt.tools.datatypes import AtomsDict
from toddgpt.tools.geom_reporter import GeomReporter


@pytest.fixture
def sample_structures():
    # Create two sample structures
    atoms1 = Atoms("H2O", positions=[(0, 0, 0), (0, 0, 1), (1, 0, 0)])
    atoms2 = Atoms("H2O", positions=[(0, 0, 0), (0, 0, 1.1), (1.1, 0, 0)])

    return [
        AtomsDict(
            numbers=atoms1.get_atomic_numbers().tolist(), positions=atoms1.positions.tolist()
        ),
        AtomsDict(
            numbers=atoms2.get_atomic_numbers().tolist(), positions=atoms2.positions.tolist()
        ),
    ]


def test_geom_reporter_initialization():
    reporter = GeomReporter()
    assert reporter.name == "geometry_reporter"
    assert "Analyzes the geometric changes" in reporter.description


def test_geom_reporter_run(sample_structures):
    reporter = GeomReporter()
    result = reporter._run(structures=sample_structures)

    assert "bond_changes" in result
    assert "angle_changes" in result
    assert "dihedral_changes" in result

def test_analyze_bonds():
    reporter = GeomReporter()
    atoms1 = Atoms("H2", positions=[(0, 0, 0), (0, 0, 0.74)])
    atoms2 = Atoms("H2", positions=[(0, 0, 0), (0, 0, 0.78)])

    changes = reporter._analyze_bonds([atoms1, atoms2])

    assert len(changes) == 1
    assert changes[0]["type"] == "bond"
    assert changes[0]["atoms"] == [0, 1]
    assert pytest.approx(changes[0]["change"]) == 0.04


def test_analyze_angles():
    expected_angle_changes = [2.72631099390626, 0.0, 0.0]
    reporter = GeomReporter()
    atoms1 = Atoms("H2O", positions=[(0, 0, 0), (0, 0, 1), (1, 0, 0)])
    atoms2 = Atoms("H2O", positions=[(0, 0, 0), (0, 0, 1), (1.1, 0, 0)])

    changes = reporter._analyze_angles([atoms1, atoms2])
    angle_changes = [change["change"] for change in changes]

    assert len(changes) > 0
    assert all(change["type"] == "angle" for change in changes)
    assert all(len(change["atoms"]) == 3 for change in changes)
    assert pytest.approx(angle_changes) == expected_angle_changes


def test_analyze_dihedrals():
    expected_dihedral_changes = [
        360.0,
        -5.710593137499643,
        -5.710593137499643,
        5.710593137499643,
        5.710593137499643,
    ]
    reporter = GeomReporter()
    atoms1 = Atoms(
        "C2H4",
        positions=[
            (0, 0, 0),
            (0, 0, 1.54),
            (-1, 0, -0.5),
            (1, 0, -0.5),
            (0, 1, 2),
            (0, -1, 2),
        ],
    )
    atoms2 = Atoms(
        "C2H4",
        positions=[
            (0, 0, 0),
            (0, 0, 1.54),
            (-1, 0, -0.5),
            (1, 0, -0.5),
            (0.1, 1, 2),
            (-0.1, -1, 2),
        ],
    )

    changes = reporter._analyze_dihedrals([atoms1, atoms2])
    dihedral_changes = [change["change"] for change in changes]

    assert len(changes) > 0
    assert all(change["type"] == "dihedral" for change in changes)
    assert all(len(change["atoms"]) == 4 for change in changes)
    assert pytest.approx(dihedral_changes) == expected_dihedral_changes
