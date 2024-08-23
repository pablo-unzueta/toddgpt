import difflib
from typing import Dict, Type

from langchain.pydantic_v1 import BaseModel, Field
from langchain.tools import BaseTool


class TerachemInputDescription(BaseModel):
    calculation_description: str = Field(..., description="Description of the calculation type")

class SuggestTerachemRunType(BaseTool):
    name = "suggest_terachem_run_type"
    description = "Suggests a TeraChem run mode based on the description of the calculation type."
    args_schema: Type[BaseModel] = TerachemInputDescription

    calculation_type_map: Dict[str, str] = {
        "energy": "energy",
        "single point": "energy",
        "single-point": "energy",
        "gradient": "gradient",
        "force": "gradient",
        "optimization": "minimize",
        "dynamics": "md",
        "molecular dynamics": "md",
        "MD simulation": "md",
        "MD trajectory": "md",
        "conical": "conical",
        "minimal energy conical intersection": "conical",
        "meci": "conical",
        "ts": "ts",
        "transition state": "ts",
        "nudged elastic band": "neb",
        "neb": "neb",
        "frequencies": "frequencies",
        "frequency": "frequencies",
        "hessian": "frequencies",
        "initial conditions": "initcond",
        "starting conditions": "initcond",
        "coupling": "coupling",
        "nonadiabatic coupling": "coupling",
        "non-adiabatic coupling": "coupling",
    }

    def _run(self, calculation_description: str) -> str:
        return self.suggest_input(calculation_description)

    def suggest_input(self, calculation_description: str) -> str:
        """
        Suggests a TeraChem run type based on the description of the calculation.
        
        Args:
            calculation_description (str): A description of the desired calculation.
        
        Returns:
            str: The suggested TeraChem run type.
        """
        calculation_description = calculation_description.lower()
        calculation_description = calculation_description.replace("_", " ").replace("-", " ")
        calculation_description = self.spell_check_description(calculation_description)
        
        matched_keys = [key for key in self.calculation_type_map.keys() if key in calculation_description]
        
        if matched_keys:
            def calculate_match_probability(search, key):
                search_words = set(search.split())
                key_words = set(key.split())
                overlap = len(search_words.intersection(key_words))
                total = len(search_words.union(key_words))
                return overlap / total if total > 0 else 0

            matches_with_probabilities = [
                (key, self.calculation_type_map[key], calculate_match_probability(calculation_description, key))
                for key in matched_keys
            ]
            
            # Sort by probability, highest first
            sorted_matches = sorted(matches_with_probabilities, key=lambda x: x[2], reverse=True)
            print(sorted_matches)
            # Return list of tuples: (matched_key, run_type, probability)
            return sorted_matches[0][1] if sorted_matches else "Unknown calculation type. Please provide more specific information."
        else:
            return "Unknown calculation type. Please provide more specific information."

    def spell_check_description(self, description: str) -> str:
        """
        Performs an automated spell check on the input description using the difflib library.
        
        Args:
            description (str): The input description to spell check.
        
        Returns:
            str: The spell-checked description.
        """
        # Define a list of known correct words related to TeraChem calculations
        known_words = [
            "single", "point", "energy", "calculation", "optimize", "optimization",
            "geometry", "molecular", "dynamics", "simulation", "conical", "intersection",
            "transition", "state", "nudged", "elastic", "band", "frequencies", "frequency",
            "hessian", "initial", "conditions", "coupling", "nonadiabatic", "adiabatic",
            "gradient", "forces"
        ]
        
        # Split the description into words
        words = description.lower().split()
        
        # Correct each word
        corrected_words = []
        for word in words:
            if word not in known_words:
                # Find the closest match in known_words
                matches = difflib.get_close_matches(word, known_words, n=1, cutoff=0.8)
                if matches:
                    corrected_words.append(matches[0])
                else:
                    corrected_words.append(word)
            else:
                corrected_words.append(word)
        
        # Join the corrected words back into a string
        corrected_description = " ".join(corrected_words)
        print(corrected_description)
        
        # Check for "minimal energy conical intersection" and replace with "meci"
        if "conical" in corrected_description:
            corrected_description = "meci"

        if "md" in corrected_description:
            corrected_description = "molecular dynamics"
        
        return corrected_description