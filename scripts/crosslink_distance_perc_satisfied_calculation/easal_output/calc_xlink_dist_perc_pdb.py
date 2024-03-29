import os
import sys
from Bio.PDB import PDBParser
import csv
import numpy as np

def calculate_distance(structure, model_A, model_B, residue1, residue2):
    #TODO structure argument is never used 

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
    #TODO cleanup avg_dist everywhere if not using
    avg_dist = []
    with open(xlink_file, 'r') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)
        for row in reader:
            residue1, residue2 = map(int, (row[0], row[2]))
            distance = calculate_distance(structure, model_A, model_B, residue1, residue2)
            #TODO can calculate in-line in this function instead of passing structure objects to another function 
            #TODO this is all we need distance = model_A['CA']-model_B['CA']
            
            avg_dist.append(distance)
            distances.append(f'{pdb_file} {distance}')

    with open(output_file, 'a') as output:
        output.write('\n'.join(distances) + '\n')

    return avg_dist

def main():
    chain_A, chain_B, xlink_file = sys.argv[1:4]
    xlink_filename = os.path.splitext(os.path.basename(xlink_file))[0]
    output_file = os.path.join('/home/muskaan/easal/easal_output/crosslink_distances/', f'{xlink_filename}_distances.txt')
    #TODO rename as output_xl_satisfaction 
    xl_satisfaction = os.path.join('/home/muskaan/easal/easal_output/xl_satisfaction/', f'{xlink_filename}_perc_satisfied.txt')

    for pdb_file in os.listdir(os.getcwd()):
        avg_across_xlinks = []
        #TODO avg_across is not used
        
        if pdb_file.endswith(".pdb"):
            xl_satisfied = pdb_file.split('_')[2] #TODO comment that EASAL outputs a binary string denoting "1" for satisfied and "0" for unsatisfied xlinks. 
            #TODO non-trivial things like this need to be commented, or else a new person will get confused. 
            
            avg_dist = process_csv(pdb_file, chain_A, chain_B, xlink_file, output_file)
            #TODO where is the avg_dist used here? 

            perc = 0
            for xl in xl_satisfied:
                if xl == '1':
                    perc += 1

            with open(xl_satisfaction, 'a') as perc_satisfied:
                perc_satisfied.write(f'{pdb_file} {(perc/len(xl_satisfied)) *100}\n')

if __name__ == "__main__":
    main()
