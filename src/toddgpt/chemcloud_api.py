import os
from pathlib import Path

import numpy as np
# print(f"Python path: {sys.path}")
from ase import units
from qcio import FileInput, ProgramOutput, Structure

from toddgpt.tools.datatypes import AtomsDict

try:
    from chemcloud import CCClient

    print("Successfully imported CCClient")
except ImportError as e:
    print(f"Import error: {e}")


def setup_terachem_qcio(tc_input: str, atoms_dict: AtomsDict) -> FileInput:
    # Create Structure object from atoms_dict
    structure = Structure(
        symbols=atoms_dict.symbols,
        geometry=np.array(atoms_dict.positions) / units.Bohr,
    )
    tc_input = Path(tc_input).read_text()
    xyz_str = structure.to_xyz()  # type: ignore
    # Create a FileInput object for TeraChem
    file_inp = FileInput(
        files={"tc.in": tc_input, "coords.xyz": xyz_str}, cmdline_args=["tc.in"]
    )
    return file_inp


def run_terachem(file_inp: FileInput) -> str:
    # This will write the files to disk in a temporary directory and then run
    # "terachem tc.in" in that directory.
    future_result = client.compute("terachem", file_inp, collect_files=True)
    prog_output: ProgramOutput = future_result.get()
    # ProgramOutput object containing all returned data
    return prog_output.stdout


# Create a function to initialize and configure the client
def initialize_chemcloud_client():
    client = CCClient()

    if "CHEMCLOUD_USER" not in os.environ:
        print(
            "CHEMCLOUD_USER environment variable not set. Please set it for future use."
        )
        client.configure()  # only run this if CHEMCLOUD_USER is not set
    return client


# Example usage
if __name__ == "__main__":
    client = initialize_chemcloud_client()
    print(client.hello_world("Pablo"))
    # ethylene excited state test
    geom = np.array(
        [
            [0.8529965186, -0.5975442252, 0.9074387146],
            [0.2785897792, 0.1371645653, 0.2978951136],
            [0.9930687044, 0.5757821515, -0.4385703161],
            [-0.8461571629, -0.4412009103, -0.3904973045],
            [-0.0044685103, 0.9531329855, 1.0002446570],
            [-1.2741293290, -0.6273345668, -1.3769108644],
        ]
    )
    geom_bohr = [x / units.Bohr for x in geom]
    ethylene = Structure(symbols=["H", "C", "H", "C", "H", "H"], geometry=geom_bohr)

    tc_input = Path("/Users/pablo/software/toddgpt/tc-inp/minimize.in").read_text()
    structure = Structure.open("/Users/pablo/software/toddgpt/tc-inp/starting-geom.xyz")
    xyz_str = structure.to_xyz()  # type: ignore
    file_inp = FileInput(
        files={"tc.in": tc_input, "coords.xyz": xyz_str}, cmdline_args=["tc.in"]
    )
    future_result = client.compute("terachem", file_inp, collect_files=True)
    prog_output: ProgramOutput = future_result.get()
    print(prog_output.stdout)
