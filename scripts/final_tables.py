#!/usr/bin/env python
import sys
import argparse
from Bio import SeqIO

from msfxns import *

parser = argparse.ArgumentParser(description="""

DESCRIPTION
    Takes in stuff

    Outputs all Jen's requested/final table structures...
    """, formatter_class= argparse.RawTextHelpFormatter)

parser.add_argument('-s2', '--s2', type=str, required=True,
                    help=''' S2_CLAMP_IP.txt or similar file.''')

parser.add_argument('-s2igg', '--s2igg', type=str, required=True,
                    help=''' S2_IgG_IP.txt or similar file.''')

parser.add_argument('-kc', '--kc', type=str, required=True,
                    help=''' ''')

parser.add_argument('-kcigg', '--kcigg', type=str, required=True,
                    help=''' ''')


parser.add_argument('-s2xlip', '--s2xlip', type=str, required=True,
                    help=''' ''')

parser.add_argument('-s2xlipigg', '--s2xlipigg', type=str, required=True,
                    help=''' ''')

##parser.add_argument('-x', '--xlip', action='store_true', default=False,
##                    help=''' The scores file is from a XLIP table.  S2_CLAMP_IP.XLIP.txt or similar file.
##                            Also can be used for any downstream 3-column files set up like XLIP.''')

parser.add_argument('key', metavar='key', nargs='+', 
                    help=''' unprot-proteins.molweightAndLengths or similar tab-sep file.
                            The file has protein names, MWs, lengths  in first 3 cols.
                            Can provide multiple such files (this is a positional arg).''')

##parser.add_argument('-c', '--constant', type=int, default=0,
##                   help=''' Multiply normalized values by this constant.
##                            Defaults:
##                                If normalizing by MW: 1000
##                                If normalizing by length: 1000
##                                If normalizing by numpeps: 100
##                            To get scores not multipled by a constant, set to 1: -c 1''')

parser.add_argument('-o', '--outdir', type=str, default='',
                    help='''Provide output directory. Default = current working directory.''')


args = parser.parse_args()

if not args.outdir:
    args.outdir = "./"
elif args.outdir[-1] != "/":
    args.outdir += "/"

## ESTABLISH CONSTANTS AND INDEXES FOR DIFFERENT NORMS IN KEY
c_mw = 1000
mw = 0
c_len = 1000
length = 1
c_pep = 100
numpep = 2

## OPEN and PARSE SCORES FILES
s2 = parse_score_file(get_file(args.s2))
s2igg = parse_score_file(get_file(args.s2igg))
kc = parse_score_file(get_file(args.kc))
kcigg = parse_score_file(get_file(args.kcigg))
s2xlip = parse_score_file(get_file(args.s2xlip), xlip=True)
s2xlipigg = parse_score_file(get_file(args.s2xlipigg), xlip=True)

## Get all protein names
names = get_names(l=[s2,s2igg,kc,kcigg], xl=[s2xlip, s2xlipigg])

## OPEN and PARSE KEY FILES
key = parse_key(get_key_files(args.key))


## TABLE1 - Intersection S2 and KC ##########################
## 1=ProteinName, 2=S2NormNumUniqPep, 3 =KCNormNumUniqPep (norm by length)
table1(outfile=args.outdir+'table1.lengthnorm.intersect.csv', acclist=list(get_intersect(s2,kc)), names=names, key=key, c=c_len, d=length, s2=s2, kc=kc)

## TABLE1 - UNION S2 or KC
## 1=ProteinName, 2=S2NormNumUniqPep, 3 =KCNormNumUniqPep (norm by length)
table1(outfile=args.outdir+'table1.lengthnorm.union.csv', acclist=list(get_union(s2,kc)), names=names, key=key, c=c_len, d=length, s2=s2, kc=kc)

## TABLE1 -  S2 NOT KC
## 1=ProteinName, 2=S2NormNumUniqPep, 3 =KCNormNumUniqPep (norm by length)
table1(outfile=args.outdir+'table1.lengthnorm.S2notKC.csv', acclist=list(get_uniq(s2,kc)), names=names, key=key, c=c_len, d=length, s2=s2, kc=kc)

## TABLE1 - KC NOT S2
## 1=ProteinName, 2=S2NormNumUniqPep, 3 =KCNormNumUniqPep (norm by length)
table1(outfile=args.outdir+'table1.lengthnorm.KCnotS2.csv', acclist=list(get_uniq(kc,s2)), names=names, key=key, c=c_len, d=length, s2=s2, kc=kc)





## WITH XLIP ##########################
## TABLE1 - Intersection S2 and KC
## 1=ProteinName, 2=S2NormNumUniqPep, 3 =KCNormNumUniqPep (norm by length)
table1(outfile=args.outdir+'table1.lengthnorm.withXLIP.intersect.csv', acclist=list(get_intersect(s2,kc)), names=names, key=key, c=c_len, d=length, s2=s2, kc=kc, s2xlip=s2xlip)

## TABLE1 - UNION S2 or KC
## 1=ProteinName, 2=S2NormNumUniqPep, 3 =KCNormNumUniqPep (norm by length)
table1(outfile=args.outdir+'table1.lengthnorm.withXLIP.union.csv', acclist=list(get_union(s2,kc)), names=names, key=key, c=c_len, d=length, s2=s2, kc=kc, s2xlip=s2xlip)

## TABLE1 -  S2 NOT KC
## 1=ProteinName, 2=S2NormNumUniqPep, 3 =KCNormNumUniqPep (norm by length)
table1(outfile=args.outdir+'table1.lengthnorm.withXLIP.S2notKC.csv', acclist=list(get_uniq(s2,kc)), names=names, key=key, c=c_len, d=length, s2=s2, kc=kc, s2xlip=s2xlip)

## TABLE1 - KC NOT S2
## 1=ProteinName, 2=S2NormNumUniqPep, 3 =KCNormNumUniqPep (norm by length)
table1(outfile=args.outdir+'table1.lengthnorm.withXLIP.KCnotS2.csv', acclist=list(get_uniq(kc,s2)), names=names, key=key, c=c_len, d=length, s2=s2, kc=kc, s2xlip=s2xlip)




## WITH IGG ##########################

## TABLE1 WITH IgG - Intersection S2 and KC 
## 1=ProteinName, 2=S2NormNumUniqPep, 3=S2IgG.. 4=S2clamp-igg 5,6,7=same as 2,3,4 for KCNormNumUniqPep 
table1(outfile=args.outdir+'table1.lengthnorm.withIgG.intersect.csv', acclist=list(get_intersect(s2,kc)), names=names, key=key, c=c_len, d=length, s2=s2, kc=kc, s2igg=s2igg, kcigg=kcigg)


## TABLE1 WITH IgG - UNION S2 or KC 
## 1=ProteinName, 2=S2NormNumUniqPep, 3=S2IgG.. 4=S2clamp-igg 5,6,7=same as 2,3,4 for KCNormNumUniqPep 
table1(outfile=args.outdir+'table1.lengthnorm.withIgG.union.csv', acclist=list(get_union(s2,kc)), names=names, key=key, c=c_len, d=length, s2=s2, kc=kc, s2igg=s2igg, kcigg=kcigg)

## TABLE1 -  S2 NOT KC
## 1=ProteinName, 2=S2NormNumUniqPep, 3 =KCNormNumUniqPep (norm by length)
table1(outfile=args.outdir+'table1.lengthnorm.withIgG.S2notKC.csv', acclist=list(get_uniq(s2,kc)), names=names, key=key, c=c_len, d=length, s2=s2, kc=kc, s2igg=s2igg, kcigg=kcigg)

## TABLE1 - KC NOT S2
## 1=ProteinName, 2=S2NormNumUniqPep, 3 =KCNormNumUniqPep (norm by length)
table1(outfile=args.outdir+'table1.lengthnorm.withIgG.KCnotS2.csv', acclist=list(get_uniq(kc,s2)), names=names, key=key, c=c_len, d=length, s2=s2, kc=kc, s2igg=s2igg, kcigg=kcigg)






## WITH IGG AND WITH XLIP  ##########################
## TABLE1 - Intersection S2 and KC
## 1=ProteinName, 2=S2NormNumUniqPep, 3 =KCNormNumUniqPep (norm by length)
table1(outfile=args.outdir+'table1.lengthnorm.withIgG.withXLIP.intersect.csv', acclist=list(get_intersect(s2,kc)), names=names, key=key, c=c_len, d=length, s2=s2, kc=kc, s2igg=s2igg, kcigg=kcigg, s2xlip=s2xlip, s2xlipigg=s2xlipigg)

## TABLE1 - UNION S2 or KC
## 1=ProteinName, 2=S2NormNumUniqPep, 3 =KCNormNumUniqPep (norm by length)
table1(outfile=args.outdir+'table1.lengthnorm.withIgG.withXLIP.union.csv', acclist=list(get_union(s2,kc)), names=names, key=key, c=c_len, d=length, s2=s2, kc=kc, s2igg=s2igg, kcigg=kcigg, s2xlip=s2xlip, s2xlipigg=s2xlipigg)

## TABLE1 -  S2 NOT KC
## 1=ProteinName, 2=S2NormNumUniqPep, 3 =KCNormNumUniqPep (norm by length)
table1(outfile=args.outdir+'table1.lengthnorm.withIgG.withXLIP.S2notKC.csv', acclist=list(get_uniq(s2,kc)), names=names, key=key, c=c_len, d=length, s2=s2, kc=kc, s2igg=s2igg, kcigg=kcigg, s2xlip=s2xlip, s2xlipigg=s2xlipigg)

## TABLE1 - KC NOT S2
## 1=ProteinName, 2=S2NormNumUniqPep, 3 =KCNormNumUniqPep (norm by length)
table1(outfile=args.outdir+'table1.lengthnorm.withIgG.withXLIP.KCnotS2.csv', acclist=list(get_uniq(kc,s2)), names=names, key=key, c=c_len, d=length, s2=s2, kc=kc, s2igg=s2igg, kcigg=kcigg, s2xlip=s2xlip, s2xlipigg=s2xlipigg)



##with open(args.outdir+'table2.S2.csv', 'w'):
##    outlist = ["Accession","Name","BioRep1_techrep1", "BioRep1_techrep2","BioRep1_union","BioRep1_mean","BioRep2","Union","kDA","AA_Length","num_peptides"]
##    print (",").join([ str(e) for e in outlist ])
##    for acc in s2.keys():
##        outlist = []
##        outlist.append(acc)
##        outlist.append(names[acc])
##        outlist.append(s2[acc][1])
##        outlist.append(s2[acc][2])
##        outlist.append(s2[acc][0])
##        outlist += s2[acc][3:7]
##        outlist += key[acc]
##        print (",").join([ str(e) for e in outlist ])

table2(args.outdir+'table2.S2.csv', names, key, s2)
table2(args.outdir+'table2.Kc.csv', names, key, kc)
table2(args.outdir+'table2.S2xlip.csv', names, key, s2xlip, xlip=True)
