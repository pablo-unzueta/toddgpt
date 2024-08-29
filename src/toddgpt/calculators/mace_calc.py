import os
import sys
import warnings
from typing import Dict, List, Type, Union

from ase import Atoms
from ase.optimize import LBFGS
from langchain.pydantic_v1 import BaseModel
from langchain.tools import BaseTool

from toddgpt.tools.datatypes import AtomsDict

# Suppress all warnings
warnings.filterwarnings("ignore")

# Suppress the specific PyTorch warning
os.environ["PYTHONWARNINGS"] = "ignore::FutureWarning"

# Redirect stdout and stderr to devnull
original_stdout = sys.stdout
original_stderr = sys.stderr
sys.stdout = open(os.devnull, "w")
sys.stderr = open(os.devnull, "w")

try:
    from mace.calculators import mace_off
finally:
    # Restore stdout and stderr
    sys.stdout.close()
    sys.stderr.close()
    sys.stdout = original_stdout
    sys.stderr = original_stderr


class MaceCalculatorInput(BaseModel):
    atoms_dict: AtomsDict
    run_type: str


class MaceCalculator(BaseTool):
    name = "mace_calculator"
    description = (
        "Use this tool after grabbing a geometry from PubChem to minimize the geometry. If the user asks for a geometry optimization, set the run_type to minimize_positions. "
    )
    args_schema: Type[BaseModel] = MaceCalculatorInput

    def __init__(self):
        super().__init__()

    def init_calc(self):
        return mace_off(model="small", device="cpu")

    def _run(
        self, input: MaceCalculatorInput
    ) -> Dict[str, Union[float, List[List[float]], List[float]]]:
        atoms = Atoms(symbols=input.atoms_dict.symbols, positions=input.atoms_dict.positions)
        atoms.calc = self.init_calc()

        if input.run_type == "sp_energy":
            return {"sp_energy": atoms.get_potential_energy()}
        elif input.run_type == "forces":
            return {"forces": atoms.get_forces().tolist()}
        elif input.run_type == "minimize_positions":
            initial_positions = atoms.get_positions()
            opt = LBFGS(atoms)
            opt.run(fmax=0.001)
            return {
                "sp_energy": atoms.get_potential_energy(),
                "initial_positions": initial_positions.tolist(),
                "minimize_positions": atoms.get_positions().tolist(),
            }
        else:
            raise ValueError(f"Invalid run type: {input.run_type}")
