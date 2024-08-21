import logging
from pathlib import Path
from typing import Any, Dict

import numpy as np
import requests
from ase.atoms import Atoms
from ase.io import read
from langchain.agents import tool

from toddgpt.tools.datatypes import AtomsDict
from toddgpt.tools.input_examples import TERACHEM_INPUT_EXAMPLES

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

PAGE_SIZE = 10  # Define a constant for the page size


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
def read_geometry_from_file(file_path: str) -> AtomsDict:
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
        path = Path(file_path)
        if not path.exists():
            raise FileNotFoundError(f"The file {file_path} does not exist.")

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
            "positions": list(zip(conformers["x"], conformers["y"], conformers["z"])),
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


@tool
def return_terachem_input(atoms_dict: AtomsDict) -> str:
    """
    Return a terachem input file for the given description from the user.
    """
    raise NotImplementedError("This tool is not implemented yet.")


@tool
def list_terachem_input_examples(page: int = 1) -> Dict[str, Any]:
    """
    List TeraChem input file example names with pagination.

    Args:
        page (int): The page number to retrieve (default is 1).

    Returns:
        Dict[str, Any]: A dictionary containing:
            - 'examples': List of example names for the requested page.
            - 'total_pages': Total number of pages.
            - 'current_page': Current page number.
            - 'total_examples': Total number of examples.
    """
    try:
        logger.info(f"Listing TeraChem input examples - Page {page}")
        all_examples = list(TERACHEM_INPUT_EXAMPLES.keys())
        total_examples = len(all_examples)
        total_pages = (total_examples + PAGE_SIZE - 1) // PAGE_SIZE

        if page < 1 or page > total_pages:
            return {
                "error": f"Invalid page number. Please choose a page between 1 and {total_pages}."
            }

        start_index = (page - 1) * PAGE_SIZE
        end_index = min(start_index + PAGE_SIZE, total_examples)
        page_examples = all_examples[start_index:end_index]

        return {
            "examples": page_examples,
            "total_pages": total_pages,
            "current_page": page,
            "total_examples": total_examples,
        }
    except Exception as e:
        logger.error(f"Error listing TeraChem input examples: {str(e)}")
        return {"error": str(e)}


@tool
def get_terachem_input_example(example_name: str) -> Dict[str, Any]:
    """
    Retrieve a TeraChem input file example by name.

    Args:
        example_name (str): The name of the example input file to retrieve.

    Returns:
        Dict[str, Any]: A dictionary containing:
            - 'content': The content of the requested TeraChem input file example.
            - 'name': The name of the example.
    """
    try:
        content = TERACHEM_INPUT_EXAMPLES.get(example_name)
        if content is None:
            return {"error": f"Example '{example_name}' not found."}
        return {"name": example_name, "content": content}
    except Exception as e:
        logger.error(f"Error retrieving TeraChem input example: {str(e)}")
        return {"error": str(e)}


@tool
def search_terachem_input_examples(keyword: str, page: int = 1) -> Dict[str, Any]:
    """
    Search for TeraChem input file examples containing a specific keyword, with pagination.

    Args:
        keyword (str): The keyword to search for in example names and content.
        page (int): The page number to retrieve (default is 1).

    Returns:
        Dict[str, Any]: A dictionary containing:
            - 'examples': List of example names for the requested page that match the keyword.
            - 'total_pages': Total number of pages of search results.
            - 'current_page': Current page number.
            - 'total_matches': Total number of examples matching the keyword.
    """
    try:
        keyword = keyword.lower()
        matching_examples = [
            name
            for name, content in TERACHEM_INPUT_EXAMPLES.items()
            if keyword in name.lower() or keyword in content.lower()
        ]
        total_matches = len(matching_examples)
        total_pages = (total_matches + PAGE_SIZE - 1) // PAGE_SIZE

        if page < 1 or page > total_pages:
            return {
                "error": f"Invalid page number. Please choose a page between 1 and {total_pages}."
            }

        start_index = (page - 1) * PAGE_SIZE
        end_index = min(start_index + PAGE_SIZE, total_matches)
        page_examples = matching_examples[start_index:end_index]

        return {
            "examples": page_examples,
            "total_pages": total_pages,
            "current_page": page,
            "total_matches": total_matches,
        }
    except Exception as e:
        logger.error(f"Error searching TeraChem input examples: {str(e)}")
        return {"error": str(e)}
