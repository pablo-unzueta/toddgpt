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
os.environ['PYTHONWARNINGS'] = 'ignore::FutureWarning'

# Redirect stdout and stderr to devnull
original_stdout = sys.stdout
original_stderr = sys.stderr
sys.stdout = open(os.devnull, 'w')
sys.stderr = open(os.devnull, 'w')

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
    description = "Use this tool to calculate the potential energy of a molecule using the MACE calculator."
    args_schema: Type[BaseModel] = MaceCalculatorInput

    def __init__(self):
        super().__init__()

    def init_calc(self):
        return mace_off(model="small", device='cpu')

    def _run(self, atoms_dict: AtomsDict, run_type: str) -> Dict[str, Union[float, List[List[float]], List[float]]]:
        atoms = Atoms(symbols=atoms_dict.symbols, positions=atoms_dict.positions)
        atoms.calc = self.init_calc()
        
        if run_type == "sp_energy":
            return {"sp_energy": atoms.get_potential_energy()}
        elif run_type == "forces":
            return {"forces": atoms.get_forces().tolist()}
        elif run_type == "minimize_positions":
            opt = LBFGS(atoms)
            opt.run(fmax=0.001)
            return {
                "sp_energy": atoms.get_potential_energy(),
                "minimize_positions": atoms.get_positions().tolist()
            }
        else:
            raise ValueError(f"Invalid run type: {run_type}")