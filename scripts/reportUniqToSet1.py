#!/usr/bin/env python
import sys

f1 = sys.argv[1]
f2 = sys.argv[2]



s1 = set([])
s2 = set([])

with open(f1) as f:
    lines1 = f.readlines()

with open(f2) as f:
    lines2 = f.readlines()

for line in lines2:
    line = line.strip().split("\t")
    s2.add(line[0])


for line in lines1:
    line = line.strip().split("\t")
    s1.add(line[0])
    if line[0] not in s2:
        print ("\t").join(line)


##print s1.difference(s2)
##print len(s1)
##print len(s2)
