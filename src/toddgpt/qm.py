import os
import subprocess
from string import Template
from typing import Any, Dict


class QMCalculator:
    def __init__(self, model: str, template_file: str):
        self.model = model
        self.input_file = "input.inp"
        self.output_file = "output.out"
        self.template_file = template_file

    def write_input(self, molecule: str, parameters: Dict[str, Any]) -> None:
        """Write input file for QM calculation using a template."""
        if not os.path.exists(self.template_file):
            raise FileNotFoundError(f"Template file {self.template_file} not found.")

        with open(self.template_file, "r") as template_file:
            template_content = template_file.read()

        template = Template(template_content)

        # Combine molecule and parameters into a single dictionary
        input_data = {"molecule": molecule, **parameters}

        # Substitute placeholders in the template with actual values
        input_content = template.safe_substitute(input_data)

        with open(self.input_file, "w") as f:
            f.write(input_content)

    def run_calculation(self) -> None:
        """Run QM calculation via command line."""
        try:
            subprocess.run(
                [self.model, self.input_file],
                check=True,
                capture_output=True,
                text=True,
            )
        except subprocess.CalledProcessError as e:
            print(f"Error running calculation: {e}")
            print(f"Stdout: {e.stdout}")
            print(f"Stderr: {e.stderr}")

    def read_output(self) -> Dict[str, Any]:
        """Read and parse output file from QM calculation."""
        results: Dict[str, Any] = {}
        # TODO: Implement output parsing logic
        with open(self.output_file, "r") as f:
            for line in f:
                # Parse relevant information from output file
                # Example: results['energy'] = float(line.split()[-1])
                pass
        return results

    def perform_calculation(
        self, molecule: str, parameters: Dict[str, Any]
    ) -> Dict[str, Any]:
        """Perform full QM calculation workflow."""
        self.write_input(molecule, parameters)
        self.run_calculation()
        return self.read_output()
