from Bio.PDB import MMCIFParser, PDBIO
import os
import glob

cif_folder = "/home/muskaan/EASAL/scripts/pdbs/"
cif_files = glob.glob(os.path.join(cif_folder, '*.cif'))

for cif_filename in cif_files:
    base_filename = os.path.splitext(os.path.basename(cif_filename))[0]
    pdb_filename = os.path.join(cif_folder, f"{base_filename.split('_')[1]}.pdb")

    cif_parser = MMCIFParser()
    structure = cif_parser.get_structure("structure", cif_filename)

    pdb_io = PDBIO()
    pdb_io.set_structure(structure[0])
    pdb_io.save(pdb_filename)
