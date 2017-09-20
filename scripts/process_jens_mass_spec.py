#!/usr/bin/env python
import sys
from os import system
import numpy as np
import re
import argparse
from collections import defaultdict
parser = argparse.ArgumentParser(description="""

DESCRIPTION
    
    Takes in Jen's Mass Spec files either as a text file or as an RTF file.
    --> Text files are those converted from original RTFs using "textutil -convert txt file.rtf"
    --> This python script converts rtf to text automatically, but depends on textutil (Mac OS).
        --> If not on a Mac OS, convert to txt some other way.
        --> I have not yet seen the conversion by other methods produce a different result.

    Outputs:
        1. text version of RTF
        2. Table 1 (identified proteins)
        3. Table 2 (peptide information)


    No guarantees this script will work on anything but the files it was designed for...
    """, formatter_class= argparse.RawTextHelpFormatter)


parser.add_argument('-i', "--input",
                   type=str,
                   help='''Point this at the input text file (file.txt) or RTF file (file.rtf). If RTF, specify -r or --rtf.''')

parser.add_argument('-o','--outprefix', type=str, default=False,
                    help=''' Give the output prefix for the table filesresulting from this...
                            If not specified, will retain input filename prefix (recommended).''')

parser.add_argument('-r','--rtf', action='store_true', default=False, help='''If providing as RTF file, specify -r or --rtf''')

args = parser.parse_args()


if args.rtf:
    cmd = 'textutil -convert txt ' + args.input
    system(cmd)
    args.input = args.input[:-3] + 'txt'

if not args.outprefix:
    args.outprefix = args.input[:-3] + "table"



with open(args.input) as f:
    i=0
    d = defaultdict(int)
    lines = []
    while "Grp Nr." not in f.next():
        continue
    for line in f:
        if "By the rule of parsimony" in line or "Number of proteins per group" in line:
            break
        headercount = 0
        for e in ["Grp Nr.", "Accession Number", "Protein Name","Protein Score","Unique PSMs","PSM Serial Nrs.","Other Grp.","Score (other)", "UniProt"]:
            if e in line.strip():
                headercount += 1
        if headercount == 0 and line.strip() != '':
            try:
                int(line.strip())
                lastint = line
                lines.append(line)
                for i in range(7):
                    lines.append(f.next())
            except:
                lines.append(lastint)
                lines.append(line)
                for i in range(6):
                    lines.append(f.next())

a = open(args.outprefix+".1.txt",'w')        
for i in range(0,len(lines),8):
    out = []
    for line in lines[i:i+8]:
        if line.strip() == '':
            out.append("NA")
        else:
            out.append(line.strip())
    a.write( ("\t").join(out) + "\n")
a.close()

lines = []
with open(args.input) as f:
    while "Scan Nr." not in f.next():
        continue
    for line in f:
        if "Serial (sequential)" in line or "Modification codes" in line:
            break
        headercount=0
        for e in ["Assigned Peptides", "#", "Sequence", "PTM Site", "Nr. Scans", "Mascot Score", "Expectation", "Isolated Mass", "Delta Mass", "Charge State","Matched Ions", "Scan Nr."]:
            if e in line.strip():
                headercount += 1
        if headercount == 0 and line.strip() != '':
            try:
                int(line.strip())
                lastint = line
                lines.append(line)
                for i in range(10):
                    lines.append(f.next())
            except:
                lines.append(lastint)
                lines.append(line)
                for i in range(9):
                    lines.append(f.next())
a = open(args.outprefix+".2.txt",'w')        
for i in range(0,len(lines),11):
    out = []
    for line in lines[i:i+11]:
        if line.strip() == '':
            out.append("NA")
        else:
            out.append(line.strip())
    a.write( ("\t").join(out) + "\n")
a.close()
