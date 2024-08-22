import os

import numpy as np
# print(f"Python path: {sys.path}")
from ase import units
# from langchain.agents import tool
# Import things that are needed generically
from langchain.pydantic_v1 import BaseModel, PrivateAttr
from langchain.tools import BaseTool
from qcio import FileInput, ProgramOutput, Structure

from toddgpt.tools.datatypes import AtomsDict

try:
    from chemcloud import CCClient

    print("Successfully imported CCClient")
except ImportError as e:
    print(f"Import error: {e}")

from typing import Optional, Type


class TerachemInput(BaseModel):
    tc_input: str
    atoms_dict: AtomsDict


class RunTerachem(BaseTool):
    name = "run_terachem"
    description = "Use this tool to run a terachem calculation using the ChemCloud API."
    args_schema: Type[BaseModel] = TerachemInput

    _chemcloud_client: Optional[CCClient] = PrivateAttr(default=None)

    def __init__(self):
        super().__init__()

    def _run(self, tc_input: str, atoms_dict: AtomsDict) -> str:
        return self.run_terachem(tc_input, atoms_dict)

    def setup_terachem_qcio(self, tc_input: str, atoms_dict: AtomsDict) -> FileInput:
        """
        Useful for creating a FileInput object for run_terachem.
        """
        structure = Structure(
            symbols=atoms_dict.symbols,
            geometry=np.array(atoms_dict.positions) / units.Bohr,
        )
        xyz_str = structure.to_xyz()
        file_inp = FileInput(
            files={"tc.in": tc_input, "coords.xyz": xyz_str}, cmdline_args=["tc.in"]
        )
        return file_inp

    def run_terachem(self, tc_input: str, atoms_dict: AtomsDict) -> str:
        """
        Useful for running a TeraChem calculation.
        """
        if self._chemcloud_client is None:
            self._chemcloud_client = self.initialize_chemcloud_client()

        file_inp = self.setup_terachem_qcio(tc_input, atoms_dict)
        future_result = self._chemcloud_client.compute(
            "terachem", file_inp, collect_files=True
        )
        prog_output: ProgramOutput = future_result.get()
        return prog_output.stdout

    def initialize_chemcloud_client(self) -> CCClient:
        """
        Useful for initializing the ChemCloud client.
        """
        if self._chemcloud_client is None:
            self._chemcloud_client = CCClient()

            if "CHEMCLOUD_USER" not in os.environ:
                print(
                    "CHEMCLOUD_USER environment variable not set. Please set it for future use."
                )
                self._chemcloud_client.configure()

        return self._chemcloud_client


# Example usage
# if __name__ == "__main__":
#     client = initialize_chemcloud_client()
#     print(client.hello_world("Pablo"))
#     # ethylene excited state test
#     geom = np.array(
#         [
#             [0.8529965186, -0.5975442252, 0.9074387146],
#             [0.2785897792, 0.1371645653, 0.2978951136],
#             [0.9930687044, 0.5757821515, -0.4385703161],
#             [-0.8461571629, -0.4412009103, -0.3904973045],
#             [-0.0044685103, 0.9531329855, 1.0002446570],
#             [-1.2741293290, -0.6273345668, -1.3769108644],
#         ]
#     )
#     geom_bohr = [x / units.Bohr for x in geom]
#     ethylene = Structure(symbols=["H", "C", "H", "C", "H", "H"], geometry=geom_bohr)

#     tc_input = Path("/Users/pablo/software/toddgpt/tc-inp/minimize.in").read_text()
#     structure = Structure.open("/Users/pablo/software/toddgpt/tc-inp/starting-geom.xyz")
#     xyz_str = structure.to_xyz()  # type: ignore
#     file_inp = FileInput(
#         files={"tc.in": tc_input, "coords.xyz": xyz_str}, cmdline_args=["tc.in"]
#     )
#     future_result = client.compute("terachem", file_inp, collect_files=True)
#     prog_output: ProgramOutput = future_result.get()
#     print(prog_output.stdout)
