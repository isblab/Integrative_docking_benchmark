import os
import sys
from Bio.PDB import PDBParser
import csv

def calculate_distance(structure, model_A, model_B, residue1, residue2):

    ca_atom1 = model_A[residue1]['CA']
    ca_atom2 = model_B[residue2]['CA']
    distance = ca_atom1 - ca_atom2
    return distance

def process_csv(pdb_file, chain_A, chain_B, xlink_file, output_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    model_A = structure[0][chain_A]
    model_B = structure[0][chain_B]
    distances = []
    with open(xlink_file, 'r') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)
        for row in reader:
            residue1, residue2 = map(int, (row[0], row[2]))
            distance = calculate_distance(structure, model_A, model_B, residue1, residue2)
            distances.append(f'{pdb_file} {distance}')
            #TODO unnecessary words stored in output file

    with open(output_file, 'a') as output:
        output.write('\n'.join(distances) + '\n')

def main():
    chain_A, chain_B, xlink_file = sys.argv[1:4]
    xlink_filename = os.path.splitext(os.path.basename(xlink_file))[0]
    output_file = os.path.join('/home/muskaan/easal_output_muskaan/crosslink_distances/', f'{xlink_filename}_distances.txt')
    xl_satisfaction = os.path.join('/home/muskaan/easal_output_muskaan/xl_satisfaction/', f'{xlink_filename}_perc_satisfied.txt')

    for pdb_file in os.listdir(os.getcwd()):
        if pdb_file.endswith(".pdb"):
            xl_satisfied = pdb_file.split('_')[2]
            process_csv(pdb_file, chain_A, chain_B, xlink_file, output_file)

            perc = 0
            for xl in xl_satisfied:
                if xl == '1':
                    perc += 1

            with open(xl_satisfaction, 'a') as perc_satisfied:
                perc_satisfied.write(f'{pdb_file} {perc/len(xl_satisfied) *100}\n')

if __name__ == "__main__":
    main()
