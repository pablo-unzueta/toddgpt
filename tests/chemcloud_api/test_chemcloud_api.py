import os
import tempfile
import unittest
from unittest.mock import MagicMock, patch

from qcio import FileInput

from toddgpt.chemcloud_api import setup_terachem_qcio
from toddgpt.tools.datatypes import AtomsDict


class TestSetupTerachemQcio(unittest.TestCase):
    def setUp(self):
        self.atoms_dict = AtomsDict(
            numbers=["1", "6", "1", "6", "1", "1"],
            positions=[
                [0.8529965186, -0.5975442252, 0.9074387146],
                [0.2785897792, 0.1371645653, 0.2978951136],
                [0.9930687044, 0.5757821515, -0.4385703161],
                [-0.8461571629, -0.4412009103, -0.3904973045],
                [-0.0044685103, 0.9531329855, 1.0002446570],
                [-1.2741293290, -0.6273345668, -1.3769108644],
            ],
        )

        self.tc_content = "# TeraChem input file\nmethod xyzabc\nbasis 6-31g"
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".in"
        ) as temp_file:
            temp_file.write(self.tc_content)
            self.tc_input = temp_file.name

    def tearDown(self):
        os.unlink(self.tc_input)

    @patch("toddgpt.chemcloud_api.CCClient")
    @patch("toddgpt.chemcloud_api.os.environ")
    def test_initialize_chemcloud_client(self, mock_environ, mock_CCClient):
        from toddgpt.chemcloud_api import initialize_chemcloud_client

        # Test when CHEMCLOUD_USER is set
        mock_environ.get.return_value = "test_user"
        mock_client = MagicMock()
        mock_CCClient.return_value = mock_client

        client = initialize_chemcloud_client()

        mock_CCClient.assert_called_once()
        self.assertEqual(client, mock_client)
        mock_client.configure.assert_not_called()

        # Test when CHEMCLOUD_USER is not set
        mock_environ.get.return_value = None
        mock_CCClient.reset_mock()
        mock_client.reset_mock()

        client = initialize_chemcloud_client()

        mock_CCClient.assert_called_once()
        self.assertEqual(client, mock_client)
        mock_client.configure.assert_called_once()

    @patch("pathlib.Path.read_text")
    def test_setup_terachem_qcio(self, mock_read_text):
        mock_read_text.return_value = self.tc_content

        file_inp = setup_terachem_qcio(self.tc_input, self.atoms_dict)

        self.assertIsInstance(file_inp, FileInput)
        self.assertIn("tc.in", file_inp.files)
        self.assertIn("coords.xyz", file_inp.files)
        self.assertEqual(file_inp.files["tc.in"], self.tc_content)

    #     # Check if the XYZ content is correct
    #     xyz_content = file_inp.files["coords.xyz"]
    #     self.assertIn("6", xyz_content)  # Number of atoms
    #     self.assertIn("H", xyz_content)
    #     self.assertIn("C", xyz_content)

    #     # Check if the positions are correctly converted to Bohr
    #     structure = Structure.from_xyz(xyz_content)
    #     np.testing.assert_allclose(
    #         structure.geometry * units.Bohr, self.atoms_dict.positions, rtol=1e-5
    #     )

    #     self.assertEqual(file_inp.cmdline_args, ["tc.in"])

    # @patch('toddgpt.chemcloud_api.CCClient')
    # @patch('toddgpt.chemcloud_api.run_terachem')
    # def test_run_terachem(self, mock_run_terachem, mock_CCClient):
    #     # Create a mock ProgramOutput
    #     mock_output = MagicMock(spec=ProgramOutput)
    #     mock_output.stdout = "TeraChem calculation completed successfully"

    #     # Set the return value for the mocked run_terachem function
    #     mock_run_terachem.return_value = mock_output

    #     # Create a sample FileInput
    #     file_input = FileInput(
    #         files={"tc.in": "# TeraChem input", "coords.xyz": "XYZ coordinates"},
    #         cmdline_args=["tc.in"]
    #     )

    #     # Run the function
    #     result = run_terachem(file_input)

    #     # Assert that the mocked run_terachem function was called with the correct argument
    #     mock_run_terachem.assert_called_once_with(file_input)

    #     # Check the result
    #     self.assertEqual(result.stdout, "TeraChem calculation completed successfully")

    #     # Ensure that CCClient was not actually instantiated
    #     mock_CCClient.assert_not_called()

    #     # # Assert that the future's get method was called
    #     # mock_future.get.assert_called_once()

    #     # # Check the result
    #     # self.assertEqual(result, "TeraChem calculation completed successfully")


if __name__ == "__main__":
    unittest.main()
