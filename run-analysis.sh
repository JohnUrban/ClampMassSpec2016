#!/bin/bash

MAIN=$PWD

##if False; then ## DELTE ME

# 1. Put scripts in PATH
echo STEP 1: Putting scripts in PATH
export PATH=${PWD}/scripts/:$PATH


# 2. Get uname -- if Darwin run as on Mac OS starting with RTF files; else start from txt files
echo STEP 2: Ascertaining if this is a Mac or other. 
U=`uname`

if [ $U == "Darwin" ]; then
   echo ... Darwin/MacOS detected... Running RTF pipeline
   DATADIR=data-as-rtf
else
   echo ... Darwin/MacOS NOT detected... Running TXT pipeline
   DATADIR=data-as-text
fi


# 3. get data from tarball
echo STEP 3: Getting data from appropriate tarball
tar -xzf ${DATADIR}.tar.gz

# 4. enter data directory and process RTF files in each folder (only in Darwin MODE).
#       Covert to Text format
#       Text output is messy and contains information for 2 distinct tables that needed to be separated.
#       Parse text file to make 2 clean table files.
cd ${DATADIR}

if [ $DATADIR == "data-as-rtf" ]; then TYPE="--rtf"; EXT="rtf"; else TYPE=""; EXT="txt"; fi

echo "STEP 4: process mass spec data into clean tables"
for dir in RTF*/; do
    if [ -d $dir ]; then
        echo Entering $dir...
        cd $dir
        for f in *.${EXT}; do
            echo ... Working on $f ...
            process_jens_mass_spec.py -i $f ${TYPE}
        done
        cd ../
    fi
done


# 5. Merge information from replicates
FILTER='keratin|taurus|trypsin|albumin|ribosom|tubulin'

## FILTER FOR THINGS TO BE TREATED AS WORDS
FILTER2='actin'

echo STEP 5: Merging information on the number of unique peptides from technical and biological replicates...
echo ... Also filtering out any entries that contain any of the following case-insensitive strings :
echo ... $FILTER
d1=RTF_files_R1_1; 
d2=RTF_files_R1_2; 
d3=RTF_files_R2;
TABLE1=.table.1.txt
TABLE2=.table.2.txt
while read line; do 
    b=`basename $line .rtf`; 
    echo ... Working on $b ... ; 
    mergeRepInfoOnNumUniqPeps.py $d1/${b}*${TABLE1} $d1/${b}*${TABLE2} $d2/${b}*${TABLE1} $d2/${b}*${TABLE2} $d3/${b}*${TABLE1} $d3/${b}*${TABLE2} | \
      awk 'OFS="\t" {gsub(/,/,""); gsub(/\ /,"_"); gsub(/OS=/,"\t"); print $1,$2,$3,$4,$5,$6,$7,$8,$9}' | grep -v -E -i $FILTER | grep -v -E -w -i $FILTER2 | sort -k8,8nr | \
      sed 's/\(.*\)_/\1 /' > ${b}.txt 
done < input.fofn


# 6. Clean up XLIP tables
echo STEP 6: Clean up XLIP tables ...
echo ... Also filtering out any entries that contain any of the following case-insensitive strings :
echo ... $FILTER ...as strings found anywhere...
echo ... AND $FILTER2 ...as string found as its own word...
for f in RTF_files_S2_CLAMP_XLIP/*${TABLE1}; do
    b=`basename $f ${TABLE1}`
    awk 'OFS="\t" {gsub(/,/,""); gsub(/\ /,"_"); gsub(/OS=/,"\t"); print $2, $5, $3}' $f  | grep -v -E -i $FILTER | grep -v -E -w -i $FILTER2 |sort -k2,2nr | \
    sed 's/\(.*\)_/\1 /' > RTF_files_S2_CLAMP_XLIP/${b}.txt
done
mv RTF_files_S2_CLAMP_XLIP/S2_CLAMP_IP.XLIP.txt .
mv RTF_files_S2_CLAMP_XLIP/S2_IgG_IP.XLIP.txt .

#7. Getting list of accessions - 
ACC=accessions.txt
echo STEP 7: Generating list of accessions that can be used to retreive proteins from UniProt...
echo "	You can upload ${ACC} here:"
echo "	http://www.uniprot.org/uploadlists/"
echo "We downloaded proteins on October 22, 2016 for this analysis and will use those to reproduce our analyses here."
echo "Please read README.md for more information on how we dealt with proteins that were deleted from UniProt or had uncommon accessions."

cat *txt | cut -f 1 | sort | uniq > $ACC
echo There are `cat accessions.txt | wc -l` $ACC
#echo "Note: This differs from README.md because various things were filtered out before making this accession list.
#echo "....$FILTER $FILTER2"

##fi ; cd data-as-rtf/ ##DELETE ME

#8. Get information about proteins
echo STEP 8: Getting Moelcular Weight, Amino Acid Length, and Number of Trypson Peptides for each protein.....
echo ... Putting in table file called: uniprot-proteins.andisoforms.info.txt

get-mol-weight.py ${MAIN}/protein_info/uniprot-proteins.andisoforms.fasta > uniprot-proteins.andisoforms.info.txt

echo ... Need to update table file with the uncommon accessions found in our mass spec data by assigning them the same info given to the more common UniProt accession.
echo ... Please read README.md for more information.

F1=uniprot-proteins.andisoforms.info.txt
F2=substituted.accessions.info.txt
while read line; do 
    A=`echo $line | awk '{print $1}'`; 
    B=`echo $line | awk '{print $2}'`; 
    grep $A $F1 | sed s/$A/$B/; 
done < ${MAIN}/protein_info/subs.txt > $F2
cat $F2 >> $F1
rm $F2

echo ... Some accessions were discontinued from UniProt database.
echo ... We used most recent versions in the history of the UniProt database before they were deleted.
echo ... Adding deleted protein info to table file...
get-mol-weight.py ${MAIN}/protein_info/deleted-proteins.fasta >> $F1

#9. Normalize protein scores.... and ranking...
mkdir tables
mkdir ranks
echo STEP 9: Normalizing scores based on number of unique peptides AND Getting Ranks of Proteins for each method

NORMS="length molweight numpeps"

for NORM in $NORMS; do
    echo ... Normalizing by ${NORM}
    for f in *_IP.txt; do
        b=`basename $f .txt`
        echo ... ... Working on $b ...
        normalize-scores.py --${NORM} -s $f uniprot-proteins.andisoforms.info.txt > tables/$b.${NORM}.txt
        sort -k8,8nr tables/$b.${NORM}.txt | awk 'OFS="\t" {print NR,$1,$8,$9}' | sort -k2,2 > ranks/$b.${NORM}.txt
    done
    for f in *XLIP.txt; do
        b=`basename $f .txt`
        echo ... ... Working on $b ...
        normalize-scores.py --xlip --${NORM} -s $f uniprot-proteins.andisoforms.info.txt > tables/$b.${NORM}.txt
        sort -k2,2nr tables/$b.${NORM}.txt | awk 'OFS="\t" {print NR,$1,$2,$3}' | sort -k2,2 > ranks/$b.${NORM}.txt
    done
done

#rank non-norm ones and mv to tables dir
for f in *_IP.txt; do
    sort -k8,8nr $f | awk 'OFS="\t" {print NR,$1,$8,$9}' | sort -k2,2 > ranks/$f
done
for f in *XLIP.txt; do
    sort -k2,2nr $f | awk 'OFS="\t" {print NR,$1,$2,$3}' | sort -k2,2 > ranks/$f
done

## Also do this with subtracting IgG from CLAMP
echo Also looking at results after associated IgGs are subtracted from CLAMP IPs...and ranking
mkdir minusIgG
subtractIgG.py -s S2_CLAMP_IP.txt -g S2_IgG_IP.txt > S2.IgGsubtracted.txt
subtractIgG.py -s Kc_CLAMP_IP.txt -g Kc_IgG_IP.txt > Kc.IgGsubtracted.txt
subtractIgG.py --xlip -s S2_CLAMP_IP.XLIP.txt -g S2_IgG_IP.XLIP.txt > S2.XLIP.IgGsubtracted.txt
for NORM in $NORMS; do
    for f in *IgGsubtracted.txt; do
        b=`basename $f .txt`
        echo ... ... Working on $b ...
        normalize-scores.py --xlip --${NORM} -s $f uniprot-proteins.andisoforms.info.txt > tables/$b.${NORM}.txt
        sort -k2,2nr tables/$b.${NORM}.txt | awk 'OFS="\t" {print NR,$1,$2,$3}' | sort -k2,2 > ranks/$b.${NORM}.txt
    done
done
#rank non-norms
for f in *IgGsubtracted.txt; do
    sort -k2,2nr $f | awk 'OFS="\t" {print NR,$1,$2,$3}' | sort -k2,2 > ranks/$f
done

mv *IgGsubtracted.txt tables/
mv *IP.txt tables/

# 10 correlations
#combine all ranks into single file

echo STEP 10: Getting correlations of different ways to process and score: 
echo "......no normalization"
echo "......length normalization" 
echo "......molecular weight normalization"
echo "......number of trypsin peptides normalization"
echo "......Subtracting IgG counts vs. Not doing so"
echo Combining ranks into single file for R analysis
for f in ranks/* ; do 
    b=`basename $f .txt`; 
    awk -v "file=$b" 'OFS="\t" {print $0,file}' $f; 
done > ranks.combined.txt
 
mkdir correlations
Rscript ${MAIN}/scripts/mscor.R ranks.combined.txt correlations/
echo "CSV files are in correlations/"
echo These can be directly imported into Excel...
echo

###fi; cd data-as-rtf/ ##DELETE ME

## 11. Make final csv output files Jen requested
echo STEP 11: Making CSV files for tables for paper
mkdir final_tables
T=tables
final_tables.py -s2 $T/S2_CLAMP_IP.txt -s2igg $T/S2_IgG_IP.txt -kc $T/Kc_CLAMP_IP.txt -kcigg $T/Kc_IgG_IP.txt \
   -s2xlip $T/S2_CLAMP_IP.XLIP.txt -s2xlipigg $T/S2_IgG_IP.XLIP.txt -o final_tables uniprot-proteins.andisoforms.info.txt

##sorting
mkdir sorted_tmp
for f in final_tables/table1*csv; do
    b=`basename $f`
    #sort -k2,2nr -t , $f > sorted_tmp/$b
    ( head -n 1 $f && tail -n +2 $f | sort -k2,2nr -t , ) > sorted_tmp/$b
done

files="final_tables/table1.lengthnorm.KCnotS2.csv final_tables/table1.lengthnorm.withXLIP.KCnotS2.csv"
for f in $files; do
    b=`basename $f`
    ( head -n 1 $f && tail -n +2 $f | sort -k3,3nr -t , ) > sorted_tmp/$b
done

files="final_tables/table1.lengthnorm.withIgG.KCnotS2.csv final_tables/table1.lengthnorm.withIgG.withXLIP.KCnotS2.csv"
for f in $files; do
    b=`basename $f`
    ( head -n 1 $f && tail -n +2 $f | sort -k7,7nr -t , ) > sorted_tmp/$b
done

for f in table2.S2.csv table2.Kc.csv; do
    #sort -k9,9 -t , final_tables/$f > sorted_tmp/$f
    ( head -n 1 final_tables/$f && tail -n +2 final_tables/$f | sort -k9,9nr -t , ) > sorted_tmp/$f
done

f=table2.S2xlip.csv
#sort -k3,3nr -t , final_tables/$f > sorted_tmp/$f 
( head -n 1 final_tables/$f && tail -n +2 final_tables/$f | sort -k3,3nr -t , ) > sorted_tmp/$f

mv sorted_tmp/* final_tables/
rmdir sorted_tmp

echo "CSV files are in final_tables/"
echo These can be directly imported into Excel...
echo

echo STEP..N... Looking at glutamine densities of selected proteins...
echo
mkdir glu
cd glu
aminoprofiler.py $MAIN/protein_info/uniprot-proteins.andisoforms.fasta $MAIN/protein_info/deleted-proteins.fasta > glutamine-profiles.txt
grep -i -E 'clamp|nocte|nelf|argonaute|ataxin|A1ZAB5|Q24572|P25822|E1JH02|E1JIF4|A8MPH9|A1ZB80|P13469|Q59E30|O61602|A1ZB83|A4V2Z1|A4V364|CG1832|Elongation_factor_1-alpha|CG17838|Uridine_kinase|lingerer|Calreticulin|Polyubiquitin|GH17761p|Polyadenylate-binding_protein|Voltage-dependent_anion-selective_channel|ADPATP_carrier_protein|Isoform_A_of_ADPATP_carrier_protein|Elongation_factor_1-alpha_1|AAA_family_protein_Bor|Histone_H2A.v' glutamine-profiles.txt > glutamine-profiles-of-interest.txt
while read line; do
 pre=`echo $line | awk '{print $1}' | tr "|" "_"`
 name=${pre}.txt
 echo $line | awk '{print $8}' | tr "," "\n" > "$name"
 Rscript $MAIN/scripts/glu.R "$name" "$pre"
done < glutamine-profiles-of-interest.txt
cd ../


cd ../
mv $DATADIR finished_analysis
echo Results in finished_analysis
echo DONE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

echo
echo If errors were thrown or analysis did not finish successfully for any reason,
echo  then pre-completed analysis can be found in completed_analysis.tar.gz

