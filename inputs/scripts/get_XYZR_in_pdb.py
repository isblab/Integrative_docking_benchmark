from Bio.PDB import PDBParser, PDBIO
import sys, os

# Atom, ca, chain, residue, x, y, z, r

def extract_calpha_with_coordinates(input_pdb_file, output_pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", input_pdb_file)

    calpha_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.has_id('CA'):
                    calpha_atoms.append(residue['CA'])

    with open(output_pdb_file, 'w') as f:
        count = 1
        for atom in calpha_atoms:
            coord = atom.get_coord()
            line = 'ATOM '
            line += str(count) + ' '
            line += str(atom.get_id()) + ' '
            line += str(atom.get_parent().get_parent().get_id()) + ' '
            line += str(atom.get_parent().get_resname()) + ' '
            line += str(atom.get_full_id()[3][1]) + ' '
            line += str(atom.get_coord()[0]) + ' '
            line += str(atom.get_coord()[1]) + ' '
            line += str(atom.get_coord()[2]) + ' '
            line += '2.8 \n'
            f.write(line)
            count += 1
            # res_num = atom.get_full_id()[3][1]
            # print(res_num[3][1])



input_file = sys.argv[1]
output_file = sys.argv[2]
extract_calpha_with_coordinates(input_file, output_file)
