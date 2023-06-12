# Standard
import argparse
from typing import Optional

# Third-party
from openbabel import pybel
import pymol


class Handler:
    """
    Handles the protein and ligand preparation for docking.
    """

    def __init__(self):
        """Constructor
        """
        self.protein: Optional[str] = None
        self.protein_output: Optional[str] = None
        self.ligand_from_structure: bool = False
        self.ligand_output: Optional[str] = None
        self.to_pdbqt: Optional[str] = None

    def parse_arguments(self):
        """
        Pareses all cmd arguments.
        """
        # Creates a parser
        parser = argparse.ArgumentParser()

        # Add arguments to the parser
        parser.add_argument("-p", "--protein")
        parser.add_argument("-op", "--output_protein", default="protein.pdb")
        parser.add_argument("-ls", "--ligand_from_structure", action="store_true")
        parser.add_argument("-ol", "--output_ligand", default="ligand.pdb")
        parser.add_argument("-pdbqt", "--to_pdbqt")

        # Save arguments to the class variables
        args = parser.parse_args()
        self.protein = args.protein
        self.protein_output = args.output_protein
        self.ligand_from_structure = args.ligand_from_structure
        self.ligand_output = args.output_ligand
        self.to_pdbqt = args.to_pdbqt

    def prepare_protein(self, *, dimer: bool = False):
        """
        Prepares the protein for the docking process.
        Saves the result into a file.

        :param: dimer: Indicates whether protein is dimer
        """
        # Load the protein into pymol
        pymol.cmd.load(self.protein, "receptor")

        # Select needed parts
        if dimer is True:
            selection = "chain A and polymer"
            selection = pymol.cmd.select("receptor", selection)

            # Save the result as into a file
            pymol.cmd.save(self.protein_output, selection)
        else:
            print("No selection has been configured yet.")

    def extract_ligand(self):
        """
        Prepares the ligand for the docking process.
        Saves the result into a file.
        """
        # Load the ligand into pymol
        pymol.cmd.load(self.ligand_from_structure, "import")

        # Select needed parts
        selection_name = "ligand"
        selection = "chain A and resname PNM"
        pymol.cmd.select(selection_name, selection)

        print(f"Ran selection: {selection}")

        # Save the result as into a file
        pymol.cmd.save(self.ligand_output, selection_name)

    @staticmethod
    def pdb_to_pdbqt(input_file: str, *, file_type: str, output_file: str = "output.pdbqt"):
        """
        Converts pdb file to pdbqt file.
        """
        molecule = next(pybel.readfile(file_type, input_file))
        molecule.write("pdbqt", output_file, overwrite=True)


def main():
    handler = Handler()
    handler.parse_arguments()

    if handler.protein is not None:
        handler.prepare_protein(dimer=True)

        if handler.ligand_from_structure is True:
            handler.extract_ligand()

    if handler.to_pdbqt is not None:
        filename_split = handler.to_pdbqt.rsplit(".", 1)
        output_file = "".join([filename_split[0], ".pdbqt"])
        handler.pdb_to_pdbqt(input_file=handler.to_pdbqt, file_type=filename_split[1], output_file=output_file)


if __name__ == '__main__':
    main()
