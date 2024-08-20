from ase.atoms import Atoms
from ase.build import bulk
from ase.constraints import UnitCellFilter
from ase.eos import calculate_eos, plot
from ase.optimize import LBFGS
from ase.units import kJ
from langchain.agents import tool

from toddgpt.tools.datatypes import AtomsDict


# Example tool
@tool
def get_atom_dict_bulk_structure(chemical_symbol: str) -> AtomsDict:
    """
    Returns bulk structure atoms dictionary for a given chemical symbol

    Args:
        chemical_symbol (str): chemical symbol as string

    Returns:
        AtomsDict: DataClass representing the atomic structure
    """
    atoms = bulk(name=chemical_symbol)
    return AtomsDict(**{k: v.tolist() for k, v in atoms.todict().items()})
