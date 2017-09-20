#!/usr/bin/env python
import sys
from os import system
import numpy as np
import re
import argparse
from collections import defaultdict

f1 = sys.argv[1] ## Kc
f2 = sys.argv[2] ## S2



with open(f1) as f:
    lines1 = f.readlines()
with open(f2) as f:
    lines2 = f.readlines()


def make_acc2line(lines):
    acc2line = {}
    for line in lines:
        line = line.strip().split("\t")
        acc2line[line[0]] = line[1:]
    return acc2line

acc2line1 = make_acc2line(lines1)
acc2line2 = make_acc2line(lines2)

allkeys = list(set(acc2line1.keys() + acc2line2.keys()))


for key in allkeys:
    try:
        part1 = acc2line1[key][:-1]
        name1 = acc2line1[key][-1]
    except:
        part1 = ['0']*7
        name1 = '-'
    try:
        part2 = acc2line2[key][:-1]
        name2 = acc2line2[key][-1]
    except:
        part2 = ['0']*7
        name2 = '-'
    print ("\t").join([str(e) for e in [key, "Kc"]+part1+["S2"]+part2+[name1,name2] ])
 
