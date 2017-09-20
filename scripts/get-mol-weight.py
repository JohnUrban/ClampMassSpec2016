#!/usr/bin/env python
import sys
import argparse
from Bio import SeqIO
import re

from msfxns import *

parser = argparse.ArgumentParser(description="""

DESCRIPTION
    Takes in protein fasta... returns tab-sep NAME + MOLWEIGHT (kDa) + AA_LENGTH
    
    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('fasta', metavar='fasta', nargs='+',
                   type= str, 
                   help='''Path(s) to fasta file(s). ''')

args = parser.parse_args()


monoisotopicMassTable = {"U":168.964203, "A":71.03711,"C":103.00919,"D":115.02694,"E":129.04259,"F":147.06841,"G":57.02146,"H":137.05891,"I":113.08406,"K":128.09496,"L":113.08406,"M":131.04049,"N":114.04293,"P":97.05276,"Q":128.05858,"R":156.10111,"S":87.03203,"T":101.04768,"V":99.06841,"W":186.07931,"Y":163.06333}
mean = sum(monoisotopicMassTable.values())/float(len(monoisotopicMassTable.values()))
monoisotopicMassTable['X'] = mean


intMassToAA = {71:["A"],103:["C"],115:["D"],129:["E"],147:["F"],57:["G"],137:["H"],113:["I","L"],128:["K","Q"],131:["M"],114:["N"],97:["P"],156:["R"],87:["S"],101:["T"],99:["V"],186:["W"],163:["Y"]}


##def peptideMass(peptide, integerMass=False):
##    mass = 0
##    if integerMass:
##        for aa in peptide:
##            mass += int(monoisotopicMassTable[aa])
##    else:
##        for aa in peptide:
##            mass += monoisotopicMassTable[aa]
##    return mass
##
##def countPeps(seq, regex='[RK]'):
##    regex = re.compile(regex)
##    length = len(seq)
##    count = 0
##    # Trypsin cleaves peptide chains mainly at the carboxyl side of the amino acids lysine or arginine, except when either is followed by proline.
##    # count number of eligible R and K (those not next to proline, P)
##    # shortcoming is also counting really short peptides that are likely not able to be assigned to a protein
##    for m in re.finditer(regex, seq):
####        print seq[m.start():m.start()+2]  #DEBUGGING CODE
##        if m.start() != length-1 and seq[m.start()+1] != 'P':
##            count += 1
##    # Add 1 to turn numKR into numPeps produced from digestion
##    count += 1
##    return count, length

## DEBUGGING CODE
##for f in args.fasta:
##    for fa in SeqIO.parse(f, 'fasta'):
##        print countPeps(str(fa.seq).upper())    


for f in args.fasta:
    for fa in SeqIO.parse(f, 'fasta'):
        kDa = peptideMass(str(fa.seq).upper())/1000.0
        count, length = countPeps(str(fa.seq).upper())       
        print ("\t").join([str(e) for e in [fa.description, kDa, length, count]])
