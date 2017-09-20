#!/usr/bin/env python
import sys
import argparse
from Bio import SeqIO
import re

from msfxns import *

parser = argparse.ArgumentParser(description="""

DESCRIPTION
    Takes in protein fasta... returns tab-sep:
    NAME
    LENGTH
    NUM specific AA requested (default Q)
    pct specific AA
    max specific AA density found in window sze
    comma-sep densities in windows of given size (slides by 1)
    
    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('fasta', metavar='fasta', nargs='+',
                   type= str, 
                   help='''Path(s) to fasta file(s). ''')

parser.add_argument('-a', '--aminoacid', type=str, default='Q',
                    help=''' Single letter amino acid symbol. Default: Q. Case-insensitive.''')
parser.add_argument('-w', '--window', type=int, default=20,
                    help=''' Window size. Default=20.''')

args = parser.parse_args()


re = '[' + args.aminoacid.upper() + args.aminoacid.lower() + ']'

for f in args.fasta:
    for fa in SeqIO.parse(f, 'fasta'):
        name = ("_").join(str(fa.description).split())
        name = ("_").join(name.split("|"))
        name = ("").join(name.split(","))
        name = ("_").join(name.split("/"))
        length = len(fa)
        starts = positionsOfAA(str(fa.seq), regex=re)
        num = len(starts)
        pct = 100.0*num/length
        counts, maxcount = densityOfAA(starts, length, window=args.window)
        csvcounts = (",").join([str(e) for e in counts])
        maxdens = float(maxcount)/args.window
        densities = [float(e)/args.window for e in counts]
        csvdens = (",").join([str(e) for e in densities])
        count, length = countPeps(str(fa.seq).upper())       
        print ("\t").join([str(e) for e in [name, length, num, pct, maxcount, maxdens, csvcounts, csvdens]])
