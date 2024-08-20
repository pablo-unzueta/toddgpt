from ase.atoms import Atoms
from ase.build import bulk
from ase.constraints import UnitCellFilter
from ase.eos import calculate_eos, plot
from ase.optimize import LBFGS
from ase.units import kJ
from langchain.agents import tool
import numpy as np
from toddgpt.tools.datatypes import AtomsDict


# Example tool
@tool
def get_distances(atoms_dict: AtomsDict) -> np.ndarray:
    atoms = Atoms(**atoms_dict.dict())
    return atoms.get_all_distances()
