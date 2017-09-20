d1=RTF_files_R1_1; 
d2=RTF_files_R1_2; 
d3=RTF_files_R2;
while read line; do b=`basename $line .rtf`; 
echo $b; 
#f1a = sys.argv[1]
#f1b = sys.argv[2]
#f2a = sys.argv[3]
#f2b = sys.argv[4]
#bf2a = sys.argv[5]
#bf2b = sys.argv[6]
echo $d1/${b}*.bullshit.1.txt $d1/${b}*.bullshit.2.txt $d2/${b}*.bullshit.1.txt $d2/${b}*.bullshit.2.txt $d3/${b}*.bullshit.1.txt $d3/${b}*.bullshit.2.txt
./numUniqPeps.FIX.py $d1/${b}*.bullshit.1.txt $d1/${b}*.bullshit.2.txt $d2/${b}*.bullshit.1.txt $d2/${b}*.bullshit.2.txt $d3/${b}*.bullshit.1.txt $d3/${b}*.bullshit.2.txt > ${b}.2.txt; 
done < input.fofn
