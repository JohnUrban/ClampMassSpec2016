#!/bin/bash
export PATH=../../analysisofmassspecconvertingtotableformat/:$PATH

bash cmd.2.sh 

cat *txt | cut -f 1 | sort | uniq > accessions.toget

## you dont have this file -- but I made one for in this dir called : uniprot-all.names
##cut -f 1 ~/data/sciara/genomeAssembly/for-maker/uniprotswissprot/uniprot-all.sizes.txt > uniprot-all.names

grep -f accessions.toget uniprot-all.names > names.to.get


# you dont have this function nor the superset fasta file -- I am making a fasta full of the proteins in all these txt files... proteins.fasta
#extractFastxEntries.py --fa -f ~/data/sciara/genomeAssembly/for-maker/uniprotswissprot/uniprot-all.fasta -n names.to.get > proteins.fasta
## only found 168 of 374 in the uniprot database....

## so instead used accessions on uniprot download server -- got all:
## 373 out of 373 UniProtKB AC/ID identifiers were successfully mapped to 329 UniProtKB IDs in the table below.
## had to change O16797-4 to O16797 to get the protein -- both are ribosomal...
#unprot-proteins.fasta


./get-mol-weight.py unprot-proteins.fasta > unprot-proteins.molweightAndLengths


## also got this: full-dmel-proteome-uniprot.fasta 
## http://www.uniprot.org/proteomes/?query=drosophila+melanogaster&sort=score

./get-mol-weight.py full-dmel-proteome-uniprot.fasta > full-dmel-proteome-uniprot.molweightAndLengths


## finally got this:
## http://www.uniprot.org/uniprot/?query=yourlist:M20161021C2335653E4FA1B8AECF5153189FA788F2F10038&sort=yourlist:M20161021C2335653E4FA1B8AECF5153189FA788F2F10038&columns=yourlist(M20161021C2335653E4FA1B8AECF5153189FA788F2F10038),id,entry%20name,reviewed,protein%20names,genes,organism,length
## unprot-proteins.andisoforms.fasta
## same as first above -- but also inclduded isoforms in download

./get-mol-weight.py unprot-proteins.andisoforms.fasta > unprot-proteins.andisoforms.molweightAndLengths


Some accessions go changed by uniprot so have to change back and include:

## DELETED FROM UNIPROT...
Q9VHL3


grep P78385 unprot-proteins.andisoforms.molweightAndLengths | awk '{sub(/P78385/,"Q6NT21"); print}' > subtituted.accessions.molweightAndLengths
grep Q61765 unprot-proteins.andisoforms.molweightAndLengths | awk '{sub(/Q61765/,"A2A5Y0"); print}' >> subtituted.accessions.molweightAndLengths
grep Q5KR48 unprot-proteins.andisoforms.molweightAndLengths | awk '{sub(/Q5KR48/,"Q3SX28"); print}' >> subtituted.accessions.molweightAndLengths
...
...and a ton more!

Norm example:
./normalize-scores.py -s S2_CLAMP_IP.2.txt unprot-proteins.andisoforms.molweightAndLengths subtituted.accessions.molweightAndLengths full-dmel-proteome-uniprot.molweightAndLengths | sort -k8,8nr | grep -v -E -i 'keratin|ribon|bovine|trypsin|adp|ribosom|Elongation|Heat|actin|histone|serum|translation' | less
./normalize-scores.py --length -s S2_CLAMP_IP.2.txt unprot-proteins.andisoforms.molweightAndLengths subtituted.accessions.molweightAndLengths full-dmel-proteome-uniprot.molweightAndLengths | sort -k8,8nr | grep -v -E -i 'keratin|ribon|bovine|trypsin|adp|ribosom|Elongation|Heat|actin|histone|serum|translation' | less



  848  grep P78385 unprot-proteins.andisoforms.molweightAndLengths | awk '{sub(/P78385/,"Q6NT21"); print}' > subtituted.accessions.molweightAndLengths
  854  grep Q61765 unprot-proteins.andisoforms.molweightAndLengths | awk '{sub(/Q61765/,"A2A5Y0"); print}' >> subtituted.accessions.molweightAndLengths
  860  grep Q5KR48 unprot-proteins.andisoforms.molweightAndLengths | awk '{sub(/Q5KR48/,"Q3SX28"); print}' >> subtituted.accessions.molweightAndLengths
  866  grep Q9W4M7 unprot-proteins.andisoforms.molweightAndLengths | awk '{sub(/Q9W4M7/,"E1JJD4"); print}' >> subtituted.accessions.molweightAndLengths
  869  grep E1JJD6 unprot-proteins.andisoforms.molweightAndLengths | awk '{sub(/E1JJD6/,"E1JJD7"); print}' >> subtituted.accessions.molweightAndLengths
  888  #grep  unprot-proteins.andisoforms.molweightAndLengths | awk '{sub(/E1JJD6/,"E1JJD7"); print}' >> subtituted.accessions.molweightAndLengths
  892  grep Q8WTI8 unprot-proteins.andisoforms.molweightAndLengths | awk '{sub(/Q8WTI8/,"Q9VHL3"); print}' >> subtituted.accessions.molweightAndLengths
  901  grep P55830 unprot-proteins.andisoforms.molweightAndLengths | awk '{sub(/P55830/,"H9XVL9"); print}' >> subtituted.accessions.molweightAndLengths
  907  grep Q7KNS5 unprot-proteins.andisoforms.molweightAndLengths | awk '{sub(/Q7KNS5/,"A1ZAT7"); print}' >> subtituted.accessions.molweightAndLengths
  909  grep Q15323 unprot-proteins.andisoforms.molweightAndLengths | awk '{sub(/Q15323/,"Q9UE12"); print}' >> subtituted.accessions.molweightAndLengths
  919  grep A0A023GQA5 unprot-proteins.andisoforms.molweightAndLengths | awk '{sub(/A0A023GQA5/,"E1JHL5"); print}' >> subtituted.accessions.molweightAndLengths
  923  grep Q9W392 unprot-proteins.andisoforms.molweightAndLengths | awk '{sub(/Q9W392/,"A4V464"); print}' >> subtituted.accessions.molweightAndLengths
  931  grep Q8N1N4 unprot-proteins.andisoforms.molweightAndLengths | awk '{sub(/Q8N1N4/,"Q7RTT2"); print}' >> subtituted.accessions.molweightAndLengths
  937  grep Q9NSB2 unprot-proteins.andisoforms.molweightAndLengths | awk '{sub(/Q9NSB2/,"Q6ISB0"); print}' >> subtituted.accessions.molweightAndLengths
  942  grep Q29443 unprot-proteins.andisoforms.molweightAndLengths | awk '{sub(/Q29443/,"Q0IIK2"); print}' >> subtituted.accessions.molweightAndLengths
  944  grep Q7PL76 unprot-proteins.andisoforms.molweightAndLengths | awk '{sub(/Q7PL76/,"Q7PL77"); print}' >> subtituted.accessions.molweightAndLengths
  949  grep Q9VN21 unprot-proteins.andisoforms.molweightAndLengths | awk '{sub(/Q9VN21/,"Q8IPM3"); print}' #>> subtituted.accessions.molweightAndLengths
  950  grep Q9VN21 unprot-proteins.andisoforms.molweightAndLengths | awk '{sub(/Q9VN21/,"Q8IPM3"); print}' >> subtituted.accessions.molweightAndLengths


Had to add some proteins manually to protein_info/unprot-proteins.andisoforms.fasta 
>tr|A0A023GQA5|A0A023GQA5_DROME Lethal (2) 37Cc, isoform D OS=Drosophila melanogaster GN=l(2)37Cc PE=1 SV=1 ADDED_THIS_MANUALLY_AS_REPLACEMENT_FOR_E1JHL5
>tr|Q8WTI8|Q8WTI8_DROME PFTAIRE-interacting factor 1A OS=Drosophila melanogaster GN=Pif1A PE=1 SV=1 ADDED_THIS_MANUALLY_AS_REPLACEMENT_FOR_Q9VHL3
