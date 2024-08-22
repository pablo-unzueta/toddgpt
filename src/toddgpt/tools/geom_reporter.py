# Class for analyzing geometric changes from a trajectory or optimization 
# add reporting to the user what the changes are in understandable terms

from typing import Dict, List, Type, Union

from ase import Atoms
from langchain.pydantic_v1 import BaseModel
from langchain.tools import BaseTool

from toddgpt.tools.datatypes import AtomsDict


class GeomReporterInput(BaseModel):
    structures: List[AtomsDict]


class GeomReporter(BaseTool):
    name = "geometry_reporter"
    description = "Analyzes the geometric changes from a trajectory or optimization"
    args_schema: Type[BaseModel] = GeomReporterInput

    def _run(self, structures: List[AtomsDict]) -> Dict[str, List[Dict[str, Union[str, List[int], float]]]]:
        atoms_list = [Atoms(symbols=structure.symbols, positions=structure.positions) for structure in structures]
        
        bond_changes = self._analyze_bonds(atoms_list)
        angle_changes = self._analyze_angles(atoms_list)
        dihedral_changes = self._analyze_dihedrals(atoms_list)
        
        return {
            "bond_changes": bond_changes,
            "angle_changes": angle_changes,
            "dihedral_changes": dihedral_changes
        }

    def _analyze_bonds(self, atoms_list: List[Atoms]) -> List[Dict[str, Union[str, List[int], float]]]:
        initial_atoms = atoms_list[0]
        final_atoms = atoms_list[-1]
        
        changes = []
        for first_index in range(len(initial_atoms)):
            for second_index in range(first_index+1, len(initial_atoms)):
                initial_distance = initial_atoms.get_distance(first_index, second_index)
                final_distance = final_atoms.get_distance(first_index, second_index)
                change = final_distance - initial_distance
                changes.append({"type": "bond", "atoms": [first_index, second_index], "change": change})
        
        sorted_changes = sorted(changes, key=lambda x: abs(x["change"]), reverse=True)
        return sorted_changes[:5]  # Report top 5 changes

    def _analyze_angles(self, atoms_list: List[Atoms]) -> List[Dict[str, Union[str, List[int], float]]]:
        initial_atoms = atoms_list[0]
        final_atoms = atoms_list[-1]
        
        changes = []
        for first_index in range(len(initial_atoms)):
            for second_index in range(len(initial_atoms)):
                for third_index in range(second_index+1, len(initial_atoms)):
                    if first_index != second_index and first_index != third_index:
                        initial_angle = initial_atoms.get_angle(first_index, second_index, third_index)
                        final_angle = final_atoms.get_angle(first_index, second_index, third_index)
                        change = final_angle - initial_angle
                        changes.append({"type": "angle", "atoms": [first_index, second_index, third_index], "change": change})
        
        sorted_changes = sorted(changes, key=lambda x: abs(x["change"]), reverse=True)
        return sorted_changes[:5]  # Report top 5 changes

    def _analyze_dihedrals(self, atoms_list: List[Atoms]) -> List[Dict[str, Union[str, List[int], float]]]:
        initial_atoms = atoms_list[0]
        final_atoms = atoms_list[-1]
        
        changes = []
        for first_index in range(len(initial_atoms)):
            for second_index in range(len(initial_atoms)):
                for third_index in range(len(initial_atoms)):
                    for fourth_index in range(third_index+1, len(initial_atoms)):
                        if len(set([first_index, second_index, third_index, fourth_index])) == 4:
                            initial_dihedral = initial_atoms.get_dihedral(first_index, second_index, third_index, fourth_index)
                            final_dihedral = final_atoms.get_dihedral(first_index, second_index, third_index, fourth_index)
                            change = final_dihedral - initial_dihedral
                            changes.append({"type": "dihedral", "atoms": [first_index, second_index, third_index, fourth_index], "change": change})
        
        sorted_changes = sorted(changes, key=lambda x: abs(x["change"]), reverse=True)
        return sorted_changes[:5]  # Report top 5 changes

