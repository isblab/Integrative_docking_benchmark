import os
import sys
from Bio.PDB import PDBParser
import csv

def calculate_distance(pdb_file, chain_A, chain_B, residue1, residue2):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file) 
    #TODO same problem here as in IMP script. parsing the structure for each xlink obtained 
    # TODO write down the algo on paper first, spot inefficiencies and then program it
    model_A = structure[0][chain_A]
    model_B = structure[0][chain_B]
    ca_atom1 = model_A[residue1]['CA']
    ca_atom2 = model_B[residue2]['CA']
    distance = ca_atom1 - ca_atom2
    return distance

def process_csv(pdb_file, chain_A, chain_B, xlink_file, output_file):
    distances = []
    with open(xlink_file, 'r') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)
        for row in reader:
            residue1, residue2 = map(int, (row[0], row[2]))
            distance = calculate_distance(pdb_file, chain_A, chain_B, residue1, residue2)
            distances.append(f'{pdb_file}_Distance between residues {residue1} and {residue2}: {distance}')
            #TODO unnecessary words stored in output file 

    with open(output_file, 'a') as output:
        output.write('\n'.join(distances) + '\n')

def main():
    chain_A, chain_B, xlink_file = sys.argv[1:4]
    output_directory = os.path.expanduser('~/easal_output_muskaan/crosslink_distances/')
    os.makedirs(output_directory, exist_ok=True)
    xlink_filename = os.path.splitext(os.path.basename(xlink_file))[0]
    output_file = os.path.join(output_directory, f'{xlink_filename}_distances.txt')

    for pdb_file in os.listdir(os.getcwd()):
        if pdb_file.endswith(".pdb"):
            process_csv(pdb_file, chain_A, chain_B, xlink_file, output_file)
            #TODO change file name to something more meaningful such as calc xlink distances for pdb

if __name__ == "__main__":
    main()
