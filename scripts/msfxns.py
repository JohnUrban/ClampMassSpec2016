import sys
import argparse
from Bio import SeqIO
import re
from os import system
import numpy as np
from collections import defaultdict, deque


## FROM MERGE REP INFO ON NUM UNIQ PEPS
def process_protein(pf, acc2ints, acc2name):
    with open(pf) as f:
        for line in f:
            line = line.strip().split("\t")
            acc = line[1]
            peps = set([int(e) for e in line[5].split(" ")])
            acc2ints[acc] = acc2ints[acc].union(peps)
            name = line[2]
            if acc2name[acc]:
                if acc2name[acc] != name:
                    sys.stderr.write("accession: " + str(acc) + "had more than one name associated with it.\n")
                    acc2name[acc] += ","+name
            else:
                acc2name[acc] = name
    return acc2ints, acc2name

def process_peptide(pf, ints2peps):
    with open(pf) as f:
        for line in f:
            line = line.strip().split("\t")
            i = int(line[0])
            pep = line[1]
            ints2peps[i] = pep
    return ints2peps

## FROM GET_MOL_WEIGHT
monoisotopicMassTable = {"U":168.964203, "A":71.03711,"C":103.00919,"D":115.02694,"E":129.04259,"F":147.06841,"G":57.02146,"H":137.05891,"I":113.08406,"K":128.09496,"L":113.08406,"M":131.04049,"N":114.04293,"P":97.05276,"Q":128.05858,"R":156.10111,"S":87.03203,"T":101.04768,"V":99.06841,"W":186.07931,"Y":163.06333}
mean = sum(monoisotopicMassTable.values())/float(len(monoisotopicMassTable.values()))
monoisotopicMassTable['X'] = mean
intMassToAA = {71:["A"],103:["C"],115:["D"],129:["E"],147:["F"],57:["G"],137:["H"],113:["I","L"],128:["K","Q"],131:["M"],114:["N"],97:["P"],156:["R"],87:["S"],101:["T"],99:["V"],186:["W"],163:["Y"]}

def peptideMass(peptide, integerMass=False):
    mass = 0
    if integerMass:
        for aa in peptide:
            mass += int(monoisotopicMassTable[aa])
    else:
        for aa in peptide:
            mass += monoisotopicMassTable[aa]
    return mass

def countPeps(seq, regex='[RK]'):
    regex = re.compile(regex)
    length = len(seq)
    count = 0
    # Trypsin cleaves peptide chains mainly at the carboxyl side of the amino acids lysine or arginine, except when either is followed by proline.
    # count number of eligible R and K (those not next to proline, P)
    # shortcoming is also counting really short peptides that are likely not able to be assigned to a protein
    for m in re.finditer(regex, seq):
##        print seq[m.start():m.start()+2]  #DEBUGGING CODE
        if m.start() != length-1 and seq[m.start()+1] != 'P':
            count += 1
    # Add 1 to turn numKR into numPeps produced from digestion
    count += 1
    return count, length

def positionsOfAA(seq, regex='[Qq]'):
    regex = re.compile(regex)
    starts = []
    for m in re.finditer(regex, seq):
        starts.append(m.start())
    return starts

def densityOfAA(positions, length, window=20):
    densities = []
    positions = deque(positions)
    binary = []
    for i in range(length):
        if len(positions) > 0:
            if i == positions[0]:
                positions.popleft()
                binary.append(1)
            else:
                binary.append(0)
        else:
            binary.append(0)
    curr = sum(binary[:window])
    maxden = curr
    densities.append(curr)
    for i in range(1, length-window+1):
        curr -= binary[i]
        curr += binary[i+window-1]
        densities.append(curr)
        if curr > maxden:
            maxden = curr
    return densities, maxden

    

### FROM NORMALIZE_SCORES
def parse_key(keyfile):
    ''' keyfile is list of lines from file from readlines()'''
    key = {}
    for line in keyfile:
##        print line 
        line = line.strip().split('\t')
        k = line[0].split('|')
##        print k
        key[k[1]] = [float(line[1]), float(line[2]), float(line[3])]
        if k[1].endswith('-1'):
            key[k[1][:-2]] = [float(line[1]), float(line[2]), float(line[3])]
        elif k[1][-2] != '-':
            key[k[1]+'-1'] = [float(line[1]), float(line[2]), float(line[3])]
    return key




## going through scores file to normalize based on key and accessions found
def normalize_scores(scoresfile, key, c, norm):
    for line in scoresfile:
        line = line.strip().split("\t")
        try:
            d = key[line[0]][norm]
            l = []
            l.append(line[0])
            for i in range(1,8):
                l.append(c*float(line[i])/d)
            l.append(line[8])
            print ("\t").join([str(e) for e in l])
        except ValueError as e:
            print ("\t").join([line[0]] + ['NA']*7 + [line[8]])

def normalize_scores_xlip(scoresfile, key, c, norm):
    for line in scoresfile:
        line = line.strip().split("\t")
        try:
            d = key[line[0]][norm]
            l = []
            l.append(line[0])
            l.append(c*float(line[1])/d)
            l.append(line[2])
            print ("\t").join([str(e) for e in l])
        except ValueError as e:
            print ("\t").join([line[0]] + ['NA']+ [line[2]])


## FROM SUBTRACTIGG
def get_lines(s):
    with open(s) as f:
        lines = f.readlines()
    return lines

def get_clamp_info(lines,i=7,j=8):
    info = {}
    for line in lines:
        line = line.strip().split()
        info[line[0]] = [float(line[i]),line[j]]
    return info

def get_igg_counts(lines,i=7):
    info = defaultdict(int)
    for line in lines:
        line = line.strip().split()
        info[line[0]] = float(line[i])
    return info


## FINAL TABLES
def get_file(fname):
    with open(fname) as f:
        out = f.readlines()
    return out

def parse_score_file(scoresfile, xlip=False):
    scores = {}
    for line in scoresfile:
        line = line.strip().split("\t")
        l = []
        if xlip:
            l = [float(line[1]), line[2]]
        else:
            for i in range(1,8):
                l.append(float(line[i]))
            l.append(line[8])
        scores[line[0]] = l
    return scores
            
def get_key_files(keys):
    keyfiles = []
    for kf in keys:        
        with open(kf) as f:
            keyfiles += f.readlines()
    return keyfiles

def normalize(e, c, d):
    return float(c)*e/float(d)

def normalize_list(l, c, d):
    return [normalize(e,c,d) for e in l]

def get_union(s1, s2):
    ''' provide 2 score dicts'''
    return set(s1.keys()).union(set(s2.keys()))

def get_intersect(s1, s2):
    ''' provide 2 score dicts'''
    return set(s1.keys()).intersection(set(s2.keys()))

def get_uniq(s1, s2):
    ''' provide 2 score dicts'''
    return set(s1.keys()).difference(set(s2.keys()))

def get_score(scores, acc, i, c , d):
    try:
        return normalize(scores[acc][i], c, d)
    except KeyError: ## not in current dict
        return 0.0

def get_names(l, xl):
    ''' l is a list of scores dicts... '''
    ''' xl is a list of scores dicts from xlip'''
    names = {}
    for d in l:
        for acc in d.keys():
            names[acc] = d[acc][7]
    for d in xl:
        for acc in d.keys():
            names[acc] = d[acc][1]
    return names




def table1(outfile, acclist, names, key, c, d, s2, kc, s2igg=None, kcigg=None, s2xlip=None, s2xlipigg=None):
    with open(outfile, 'w') as out:
        outlist = ["Name", "S2_CLAMP"]
        if s2igg:
            outlist.append("S2_IgG")
            outlist.append("S2_CLAMP-IgG")
        outlist.append("Kc_CLAMP")
        if kcigg:
            outlist.append("Kc_IgG")
            outlist.append("Kc_CLAMP-IgG")
        if s2xlip:
            outlist.append("S2XLIP_CLAMP")
        if s2xlipigg:
            outlist.append("S2XLIP_IgG")
            outlist.append("S2XLIP_CLAMP-IgG")
        line = (",").join([ str(e) for e in outlist ])
        out.write(line+"\n")
        #LOOP
        for acc in acclist:
            outlist = []
            name = names[acc]
            outlist.append(name)
            s2score = get_score(s2, acc, 6, c, key[acc][d])
            outlist.append(s2score)
            if s2igg:
                s2iggscore = get_score(s2igg, acc, 6, c, key[acc][d])
                s2sub = s2score-s2iggscore
                outlist.append(s2iggscore)
                outlist.append(s2sub)
            kcscore = get_score(kc, acc, 6, c, key[acc][d])
            outlist.append(kcscore)
            if kcigg:
                kciggscore = get_score(kcigg, acc, 6, c, key[acc][d])
                kcsub = kcscore - kciggscore
                outlist.append(kciggscore)
                outlist.append(kcsub)
            if s2xlip:
                s2xlipscore = get_score(s2xlip, acc, 0, c, key[acc][d])
                outlist.append(s2xlipscore)
                if s2xlipigg:
                    s2xlipiggscore = get_score(s2xlipigg, acc, 0, c, key[acc][d])
                    s2xlipsub = s2xlipscore - s2xlipiggscore
                    outlist.append(s2xlipiggscore)
                    outlist.append(s2xlipsub)
            line = (",").join([ str(e) for e in outlist ])
            out.write(line+"\n")

def table2(outfile, names, key, scores, xlip=False):
    with open(outfile, 'w') as out:
        if xlip:
            outlist = ["Accession","Name","BioRep1","kDA","AA_Length","num_peptides"]
        else:
            outlist = ["Accession","Name","BioRep1_techrep1", "BioRep1_techrep2","BioRep1_union","BioRep1_mean","BioRep2","Union", "Mean", "kDA","AA_Length","num_peptides"]
        line = (",").join([ str(e) for e in outlist ])
        out.write(line+"\n")
        for acc in scores.keys():
            outlist = []
            outlist.append(acc)
            outlist.append(names[acc])
            if xlip:
                outlist.append(scores[acc][0])
            else:
                outlist.append(scores[acc][1])
                outlist.append(scores[acc][2])
                outlist.append(scores[acc][0])
                outlist += scores[acc][3:7]
            outlist += key[acc]
            line = (",").join([ str(e) for e in outlist ])
            out.write(line+"\n")
