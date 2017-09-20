#!/usr/bin/env python
import sys
import argparse
from Bio import SeqIO
from collections import defaultdict

from msfxns import *

parser = argparse.ArgumentParser(description="""

DESCRIPTION
    Takes in CLAMP IP and IgG IP
    If a protein in CLAMP shows up in IgG,
    the IgG is subtracted from CLAMP:
        CLAMP - IgG
    

    Alternatively, you can remove the protein from the CLAMP set.
    Subtraction almost certainly will result in the protein being lowly ranked
    with the advantage of keeping it in the set to observe.
    However, as this is not a quantitative MS experiment and the
    samples are not normalized to each other ... meaning subtraction is
    only a way to lower the rank of the given protein.

    
    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('-s', '--scores', type=str, required=True,
                    help=''' S2_CLAMP_IP.txt or similar file.''')

parser.add_argument('-g', '--igg', type=str, required=True,
                    help=''' Corresponding IgG file -- e.g. S2_IgG_IP.txt or similar file.''')


parser.add_argument('-x', '--xlip', action='store_true', default=False,
                    help=''' The scores file is from a XLIP table.  S2_CLAMP_IP.XLIP.txt or similar file. ''')

parser.add_argument('-r', '--remove', action='store_true', default=False,
                    help='''Alternative: If CLAMP interacting protein in IgG, remove it instead of subtracting from score...''')

args = parser.parse_args()


## DEAL WITH ARGS
if args.xlip:
    i,j = 1,2
else:
    i,j = 7,8


## DEFINE FUNCTIONS
##def get_lines(s):
##    with open(s) as f:
##        lines = f.readlines()
##    return lines
##
##def get_clamp_info(lines,i=7,j=8):
##    info = {}
##    for line in lines:
##        line = line.strip().split()
##        info[line[0]] = [float(line[i]),line[j]]
##    return info
##
##def get_igg_counts(lines,i=7):
##    info = defaultdict(int)
##    for line in lines:
##        line = line.strip().split()
##        info[line[0]] = float(line[i])
##    return info    



## EXECUTE
clamp = get_lines(args.scores)
igg = get_lines(args.igg)
clamp_info = get_clamp_info(clamp, i=i,  j=j)
igg_counts = get_igg_counts(igg, i=i)

## IF REMOVING
if args.remove:
    igg_set = set(igg.keys())

## OUTPUT
for acc in clamp_info.keys():
    if args.remove:
        if acc not in igg_set:
            print ("\t").join([str(e) for e in [acc, clamp_info[acc][0], clamp_info[acc][1]]])
    else:
##        print clamp_info[acc][0], igg_counts[acc]
        # note if acc not in igg dict, it will default to sbtracting 0
        new_clamp_count = clamp_info[acc][0] - igg_counts[acc]
        print ("\t").join([str(e) for e in [acc, new_clamp_count, clamp_info[acc][1], clamp_info[acc][0], igg_counts[acc]]])




