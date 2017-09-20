#!/usr/bin/env python
import sys
import argparse
from Bio import SeqIO
from msfxns import *

parser = argparse.ArgumentParser(description="""

DESCRIPTION
    Takes in 
    
    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('-s', '--scores', type=str, required=True,
                    help=''' S2_CLAMP_IP.txt or similar file.''')
parser.add_argument('key', metavar='key', nargs='+', 
                    help=''' unprot-proteins.molweightAndLengths or similar tab-sep file.
                            The file has protein names, MWs, lengths  in first 3 cols.
                            Can provide multiple such files (this is a positional arg).''')
normtype = parser.add_mutually_exclusive_group(required=True)
normtype.add_argument('-m', '--molweight', action='store_true', default=False,
                    help=''' normalize by molecular weight ''')

normtype.add_argument('-l', '--length', action='store_true', default=False,
                    help=''' normalize by AA length instead of MW ''')

normtype.add_argument('-n', '--numpeps', action='store_true', default=False,
                    help=''' normalize by number possible peptides created by Trypsin instead of MW. ''')

parser.add_argument('-c', '--constant', type=int, default=0,
                   help=''' Multiply normalized values by this constant.
                            Defaults:
                                If normalizing by MW: 1000
                                If normalizing by length: 1000
                                If normalizing by numpeps: 100
                            To get scores not multipled by a constant, set to 1: -c 1''')

parser.add_argument('-x', '--xlip', action='store_true', default=False,
                    help=''' The scores file is from a XLIP table.  S2_CLAMP_IP.XLIP.txt or similar file.
                            Also can be used for any downstream 3-column files set up like XLIP.''')


args = parser.parse_args()

## shorter var name and define norm type
## norm=0 for MW normalization, norm=1 for length normalization, norm=2 for numpeps normalization
if args.constant == 0:
    if args.molweight:
        c = 1000
        norm = 0
    elif args.length:
        c = 1000
        norm = 1
    elif args.numpeps:
        c = 100
        norm = 2




## opening files and dumping all their lines into lists (fine for small files)
keyfile = []
for kf in args.key:        
    with open(kf) as f:
        keyfile += f.readlines()

with open(args.scores) as f:
    scoresfile = f.readlines()

## defining function to turn info from key file into a dictionary
##def parse_key(keyfile):
##    ''' keyfile is list of lines from file from readlines()'''
##    key = {}
##    for line in keyfile:
####        print line 
##        line = line.strip().split('\t')
##        k = line[0].split('|')
####        print k
##        key[k[1]] = [float(line[1]), float(line[2]), float(line[3])]
##        if k[1].endswith('-1'):
##            key[k[1][:-2]] = [float(line[1]), float(line[2]), float(line[3])]
##        elif k[1][-2] != '-':
##            key[k[1]+'-1'] = [float(line[1]), float(line[2]), float(line[3])]
##    return key
##
##
##
##
#### going through scores file to normalize based on key and accessions found
##def normalize_scores(scoresfile, key, c):
##    for line in scoresfile:
##        line = line.strip().split("\t")
##        try:
##            d = key[line[0]][norm]
##            l = []
##            l.append(line[0])
##            for i in range(1,8):
##                l.append(c*float(line[i])/d)
##            l.append(line[8])
##            print ("\t").join([str(e) for e in l])
##        except ValueError as e:
##            print ("\t").join([line[0]] + ['NA']*7 + [line[8]])
##
##def normalize_scores_xlip(scoresfile, key, c):
##    for line in scoresfile:
##        line = line.strip().split("\t")
##        try:
##            d = key[line[0]][norm]
##            l = []
##            l.append(line[0])
##            l.append(c*float(line[1])/d)
##            l.append(line[2])
##            print ("\t").join([str(e) for e in l])
##        except ValueError as e:
##            print ("\t").join([line[0]] + ['NA']+ [line[2]])

key = parse_key(keyfile)

if args.xlip:
    normalize_scores_xlip(scoresfile, key, c, norm)
else:
    normalize_scores(scoresfile, key, c, norm)
