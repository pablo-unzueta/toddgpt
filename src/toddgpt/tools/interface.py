from ase.atoms import Atoms
from ase.build import bulk
from ase.constraints import UnitCellFilter
from ase.eos import calculate_eos, plot
from ase.optimize import LBFGS
from ase.units import kJ
from langchain.agents import tool
import numpy as np
from toddgpt.tools.datatypes import AtomsDict
import requests
import json
from typing import Dict, Any
from pathlib import Path
from ase.io import read
from typing import Optional, Type

from langchain.pydantic_v1 import BaseModel
from langchain_core.callbacks import (
    AsyncCallbackManagerForToolRun,
    CallbackManagerForToolRun,
)
from langchain_core.tools import BaseTool


        
@tool
def get_distances(atoms_dict: AtomsDict) -> np.ndarray:
    """Get all distances between atoms in the system.

    Args:
        atoms_dict (AtomsDict): AtomsDict object

    Returns:
        np.ndarray: 2D array of distances between atoms
    """
    atoms = Atoms(**atoms_dict.dict())
    return atoms.get_all_distances()

@tool
def read_geometry_from_file(path: str) -> AtomsDict:
    """Read the geometry from a specified file path if provided.

    Args:
        path (Path): Path to the file containing the geometry.

    Returns:
        AtomsDict: AtomsDict object representing the read geometry.

    Raises:
        FileNotFoundError: If the specified file does not exist.
        ValueError: If the file format is not supported or the file is invalid.
    """
    try:
        path = Path(path)
        if not path.exists():
            raise FileNotFoundError(f"The file {path} does not exist.")

        atoms = read(path)
        return AtomsDict(
            numbers=atoms.get_atomic_numbers().tolist(),
            positions=atoms.positions.tolist(),
        )
    except Exception as e:
        raise ValueError(f"Error reading geometry from {path}: {str(e)}")

@tool
def extract_molecule_from_pubchem(compound_name: str) -> Dict[str, Any]:
    """
    Extract molecule information from PubChem based on the compound name using PUG REST API if no path is provided.

    Args:
        compound_name (str): The name of the compound to search for.

    Returns:
        Dict[str, Any]: A dictionary containing the following information:
            - 'atoms': List of atomic symbols
            - 'positions': List of 3D coordinates for each atom
            - 'formula': Molecular formula
            - 'molecular_weight': Molecular weight
            - 'iupac_name': IUPAC name of the compound
            - 'smiles': SMILES representation of the molecule
    """

    try:
        # Search for the compound using PUG REST API
        search_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound_name}/JSON"
        response = requests.get(search_url)
        response.raise_for_status()
        data = response.json()

        if "PC_Compounds" not in data:
            return {"error": f"No compound found for '{compound_name}'"}

        compound = data["PC_Compounds"][0]
        cid = compound["id"]["id"]["cid"]

        # Fetch 3D coordinates
        coord_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/JSON?record_type=3d"
        coord_response = requests.get(coord_url)
        coord_response.raise_for_status()
        coord_data = coord_response.json()

        atoms = coord_data["PC_Compounds"][0]["atoms"]
        conformers = coord_data["PC_Compounds"][0]["coords"][0]["conformers"][0]

        # Fetch additional properties
        prop_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/MolecularFormula,MolecularWeight,IUPACName,IsomericSMILES/JSON"
        prop_response = requests.get(prop_url)
        prop_response.raise_for_status()
        prop_data = prop_response.json()

        properties = prop_data["PropertyTable"]["Properties"][0]

        # Prepare the return dictionary
        result = {
            "atoms": [atoms["element"][i] for i in range(len(atoms["element"]))],
            "positions": list(
                zip(conformers["x"], conformers["y"], conformers["z"])
            ),
            "formula": properties["MolecularFormula"],
            "molecular_weight": properties["MolecularWeight"],
            "iupac_name": properties["IUPACName"],
            "smiles": properties["IsomericSMILES"],
        }

        return result

    except requests.exceptions.RequestException as e:
        return {"error": f"An error occurred while fetching data: {str(e)}"}
    except KeyError as e:
        return {"error": f"An error occurred while parsing data: {str(e)}"}
    except Exception as e:
        return {"error": f"An unexpected error occurred: {str(e)}"}

# class Interface(BaseTool):
#     name = "Interface"
#     description = "Run this tool to get the distances between atoms in a system or to read a geometry from a file."

#     def _run(self, query: str) -> str:
#         ...


#     def get_distances(self, atoms_dict: AtomsDict) -> np.ndarray:
#         """Get all distances between atoms in the system.

#         Args:
#             atoms_dict (AtomsDict): AtomsDict object

#         Returns:
#             np.ndarray: 2D array of distances between atoms
#         """
#         atoms = Atoms(**atoms_dict.dict())
#         return atoms.get_all_distances()

#     def read_geometry_from_file(self, path: str) -> AtomsDict:
#         """Read the geometry from a specified file path if provided.

#         Args:
#             path (Path): Path to the file containing the geometry.

#         Returns:
#             AtomsDict: AtomsDict object representing the read geometry.

#         Raises:
#             FileNotFoundError: If the specified file does not exist.
#             ValueError: If the file format is not supported or the file is invalid.
#         """
#         try:
#             path = Path(path)
#             if not path.exists():
#                 raise FileNotFoundError(f"The file {path} does not exist.")

#             atoms = read(path)
#             return AtomsDict(
#                 numbers=atoms.get_atomic_numbers().tolist(),
#                 positions=atoms.positions.tolist(),
#             )
#         except Exception as e:
#             raise ValueError(f"Error reading geometry from {path}: {str(e)}")

#     def extract_molecule_from_pubchem(self, compound_name: str) -> Dict[str, Any]:
#         """
#         Extract molecule information from PubChem based on the compound name using PUG REST API if no path is provided.

#         Args:
#             compound_name (str): The name of the compound to search for.

#         Returns:
#             Dict[str, Any]: A dictionary containing the following information:
#                 - 'atoms': List of atomic symbols
#                 - 'positions': List of 3D coordinates for each atom
#                 - 'formula': Molecular formula
#                 - 'molecular_weight': Molecular weight
#                 - 'iupac_name': IUPAC name of the compound
#                 - 'smiles': SMILES representation of the molecule
#         """

#         try:
#             # Search for the compound using PUG REST API
#             search_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound_name}/JSON"
#             response = requests.get(search_url)
#             response.raise_for_status()
#             data = response.json()

#             if "PC_Compounds" not in data:
#                 return {"error": f"No compound found for '{compound_name}'"}

#             compound = data["PC_Compounds"][0]
#             cid = compound["id"]["id"]["cid"]

#             # Fetch 3D coordinates
#             coord_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/JSON?record_type=3d"
#             coord_response = requests.get(coord_url)
#             coord_response.raise_for_status()
#             coord_data = coord_response.json()

#             atoms = coord_data["PC_Compounds"][0]["atoms"]
#             conformers = coord_data["PC_Compounds"][0]["coords"][0]["conformers"][0]

#             # Fetch additional properties
#             prop_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/MolecularFormula,MolecularWeight,IUPACName,IsomericSMILES/JSON"
#             prop_response = requests.get(prop_url)
#             prop_response.raise_for_status()
#             prop_data = prop_response.json()

#             properties = prop_data["PropertyTable"]["Properties"][0]

#             # Prepare the return dictionary
#             result = {
#                 "atoms": [atoms["element"][i] for i in range(len(atoms["element"]))],
#                 "positions": list(
#                     zip(conformers["x"], conformers["y"], conformers["z"])
#                 ),
#                 "formula": properties["MolecularFormula"],
#                 "molecular_weight": properties["MolecularWeight"],
#                 "iupac_name": properties["IUPACName"],
#                 "smiles": properties["IsomericSMILES"],
#             }

#             return result

#         except requests.exceptions.RequestException as e:
#             return {"error": f"An error occurred while fetching data: {str(e)}"}
#         except KeyError as e:
#             return {"error": f"An error occurred while parsing data: {str(e)}"}
#         except Exception as e:
#             return {"error": f"An unexpected error occurred: {str(e)}"}
