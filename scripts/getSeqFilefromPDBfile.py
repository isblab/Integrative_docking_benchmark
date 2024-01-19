#!/usr/bin/python
# -*- coding: utf-8 -*-

# This program uses the Biopython distribution. (http://www.biopython.org)

# 01/28/2010
# @Author Shruthi Viswanath

import os,sys,math,numpy,string
from subprocess import *
from optparse import OptionParser
import glob
from Bio.PDB import *
import warnings

def getAA_Alphabet(resnam):
	
	if resnam=='ALA':
		return 'A'
	elif resnam=='CYS':
		return 'C'
	elif resnam=='ASP':
		return 'D'
	elif resnam=='GLU':
		return 'E'
	elif resnam=='PHE':
		return 'F'
	elif resnam=='GLY':
		return 'G'
	elif resnam=='HIS':
		return 'H'
	elif resnam=='ILE':
		return 'I'
	elif resnam=='LYS':
		return 'K'
	elif resnam=='LEU':
		return 'L'
	elif resnam=='MET':
		return 'M'
	elif resnam=='ASN':
		return 'N'
	elif resnam=='PRO':
		return 'P'
	elif resnam=='GLN':
		return 'Q'
	elif resnam=='ARG':
		return 'R'
	elif resnam=='SER':
		return 'S'
	elif resnam=='THR':
		return 'T'
	elif resnam=='VAL':
		return 'V'
	elif resnam=='TRP':
		return 'W'
	elif resnam=='TYR':
		return 'Y'
	else:
		return 'G'
	# default is GLYcine	
	
def makeSeqFile():

	warnings.filterwarnings('ignore')


	# making the sequence file with all residues in lower case.
	inputFile=sys.argv[1]
	chainName=sys.argv[2]
	outputFile=sys.argv[3]

	# Get the structure
	parser=PDBParser()
	structure=parser.get_structure('pdbStruct', inputFile)
	
	# Get the list of residues in current structure
	op=open(outputFile,'w')
	
	for res in structure[0][chainName]:
		nam=res.get_resname()
		alphabet=getAA_Alphabet(nam)
				
		op.write(alphabet)
	
	op.close()
	
	
    
# MAIN
if __name__ == '__main__':
    makeSeqFile()
    

