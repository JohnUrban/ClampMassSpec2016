#!/usr/bin/env python
import sys
from os import system
import numpy as np
import re
import argparse
from collections import defaultdict

from msfxns import *

########################
'''Usage'''
########################
if len(sys.argv) == 1:
    print
    print """
Usage: mergeRepInfoOnNumUniqPeps.py br1-t1-pro  br1-t1-pep  br1-t2-pro  br1-t2-pep  br2-pro  br2-pep

mergeRepInfoOnNumUniqPeps -- i.e. merge Info about Uniq Peps from tech and bio Replicates....
This is only for procesing the non-XLIP table files (output by process_jens_mass_spec.py)
The XLIP files do not have technical replicates nor biological replicates, so no merging of information needs to take place.
Instead awk and other cmdline tools can be used on the table files for XLIP...
"""
    print
    quit()



########################
'''Files'''
########################
## Protein (a) and peptide (b) files for BioRep1 Tech Reps 1
f1a = sys.argv[1]
f1b = sys.argv[2]

## Protein (a) and peptide (b) files for BioRep1 Tech Reps 2
f2a = sys.argv[3]
f2b = sys.argv[4]

## Protein (a) and peptide (b) files for BioRep2
bf2a = sys.argv[5]
bf2b = sys.argv[6]



########################
'''Functions'''
########################
##def process_protein(pf, acc2ints, acc2name):
##    with open(pf) as f:
##        for line in f:
##            line = line.strip().split("\t")
##            acc = line[1]
##            peps = set([int(e) for e in line[5].split(" ")])
##            acc2ints[acc] = acc2ints[acc].union(peps)
##            name = line[2]
##            if acc2name[acc]:
##                if acc2name[acc] != name:
##                    sys.stderr.write("accession: " + str(acc) + "had more than one name associated with it.\n")
##                    acc2name[acc] += ","+name
##            else:
##                acc2name[acc] = name
##    return acc2ints, acc2name
##
##def process_peptide(pf, ints2peps):
##    with open(pf) as f:
##        for line in f:
##            line = line.strip().split("\t")
##            i = int(line[0])
##            pep = line[1]
##            ints2peps[i] = pep
##    return ints2peps



########################
'''Dictionaries'''
########################
acc2name = defaultdict(str)
acc2ints1 = defaultdict(set)
acc2ints2 = defaultdict(set)
acc2intsBR2 = defaultdict(set)
ints2peps1 = defaultdict(set)
ints2peps2 = defaultdict(set)
ints2pepsBR2 = defaultdict(set)
acc2peps1 = defaultdict(set)
acc2peps2 = defaultdict(set)
acc2peps = defaultdict(set) ## both br1 tech reps
acc2pepsBR2 = defaultdict(set)
acc2pepsBioReps = defaultdict(set)



########################
'''Process protein files'''
########################

#### BIOREP1 TECH2
acc2ints1, acc2name = process_protein(f1a, acc2ints1, acc2name)


#### BIOREP1 TECH2
acc2ints2, acc2name = process_protein(f2a, acc2ints2, acc2name)



#### BIOREP2 
acc2intsBR2, acc2name = process_protein(bf2a, acc2intsBR2, acc2name)
            
        
########################
'''Process peptide files'''
########################

#### BIOREP1 TECH1
ints2peps1 = process_peptide(f1b, ints2peps1)

#### BIOREP1 TECH2
ints2peps2 = process_peptide(f2b, ints2peps2)        

#### BIOREP2 
ints2pepsBR2 = process_peptide(bf2b, ints2pepsBR2)


####################
''' Map Accessions to Peptides.'''
####################

#### BIOREP1 TECH1
for acc in acc2ints1.keys():
    for pepint in acc2ints1[acc]:
        acc2peps[acc].add( ints2peps1[pepint] )
        acc2peps1[acc].add( ints2peps1[pepint] )
        acc2pepsBioReps[acc].add( ints2peps1[pepint] )



#### BIOREP1 TECH2
for acc in acc2ints2.keys():
    for pepint in acc2ints2[acc]:
        acc2peps[acc].add( ints2peps2[pepint] )
        acc2peps2[acc].add( ints2peps2[pepint] )
        acc2pepsBioReps[acc].add( ints2peps2[pepint] )

#### BIOREP2 
for acc in acc2intsBR2.keys():
    for pepint in acc2intsBR2[acc]:
        acc2pepsBR2[acc].add( ints2pepsBR2[pepint] )
        acc2pepsBioReps[acc].add( ints2pepsBR2[pepint] )




########################
'''Output'''
########################
for acc in acc2pepsBioReps.keys():
    biorep1_n1 = len(acc2peps1[acc])
    biorep1_n2 = len(acc2peps2[acc])
    biorep1_n_mu = (biorep1_n1+biorep1_n2)/2.0
    biorep1_union = len(acc2peps[acc])
    biorep2_n = len(acc2pepsBR2[acc])
    bioreps_union = len(acc2pepsBioReps[acc])
    bioreps_n_mu = (biorep1_n_mu + biorep2_n)/2.0
    print ("\t").join([ acc, str(biorep1_union), str(biorep1_n1),  str(biorep1_n2),  str(biorep1_n_mu), str(biorep2_n), str(bioreps_union), str(bioreps_n_mu), acc2name[acc] ])






