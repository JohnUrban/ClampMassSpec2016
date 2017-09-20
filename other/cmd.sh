MAIN=/Users/johnurban/searchPaths/github/ClampMassSpec2016
grep -i -E 'clamp|nocte|nelf|argonaute|ataxin' glutamine-profiles.txt > glutamine-profiles-of-interest.txt
while read line; do
 pre=`echo $line | awk '{print $1}' | tr "|" "_"`
 name=${pre}.txt
 echo $line | awk '{print $8}' | tr "," "\n" > "$name"
 Rscript $MAIN/scripts/glu.R "$name" "$pre"
done < glutamine-profiles-of-interest.txt
