# ClampMassSpec2016
Scripts and data needed to reproduce Mass Spec results in CLAMP paper

Regenerate Analysis by running:
$ bash run_analysis.sh

Depends on Python and R being installed.
Uses lattice package in R, which should be part of normal installation.
If not, try the following at R prompt:
  install.packages("lattice")

Depends on biopython package. If you do not have it, try:
pip install biopython

If you do not have pip, first try:
sudo easy_install pip

This runs on Mac OS well. It should run on Linux too.
If it fails early, try editing line 14 to: U=anything
And run again:
$ bash run_analysis.sh

Ultimately, the finished analysis directory will either be called data-as-rtf/ or data-as-text/.

If that fails or you see error messages while running the analysis, then you can view the code to see the analysis and see completed analysis results by:
tar -xzf completed_analysis.tar.gz



DESCRIPTION
1. Data was provided as Rich Text Format (RTF)
2. Convert to text. Text output is messy and contains information for 2 distinct tables that needed to be separated.
3. Parse text file to make 2 clean table files.
4. BioRep1 has 2 technical replicates (a second injection of BioRep1 was performed simply because the first did not yield a lot of peptide data).
5. The information from 2 technical replicates needed to be integrated, and the integrated BioRep1 information needed to be integrated with BioRep2.
6. We focused on number of unique peptides. One way to merge information is to create the union (set) of unique peptides. 
	Another way to merge information is to take the mean number of unique peptides between the two.
	If one set has 5 unique peptides and the other has 0, the mean is 2.5.
	The union approach is biased against short proteins. If the number of unique peptides is already saturated in set1, set2 will not contribute any information.
	Conversely, long proteins that are not saturated in set 1 will be more likeley to get new information from set 2.
	The mean approach offsets this bias a little bit. However, in both, long proteins are arbitrarily able to have a higher number of unique peptides than short proteins.
	Thus, one needs to normalize the number of unique peptides observed for each protein by the number of possible unique peptides it could have.
	Two proxies for the number of possible unique peptides are molecular weight and amino acid length of the protein.
	Molecular weight applies different normalization weights dependent on amino acid composition, which may or may not be appropriate.
	Amino acid length applies a uniform normalization weight to all proteins independent of amino acid composition.
	We tried both ways.
7. For weights we obtained monoisotopic masses for each amino acid including U. For X we used the mean monoisotopic mass.	
8. We obtained protein sequences from UniProt with 696 unique accessions in our datasets (from CLAMP_IP and IgG_IP files, NOT Input files) 
	URL: http://www.uniprot.org/uploadlists/
   	Downloaded on October 22, 2016.
   	Options: Select From: UniProtKB AC/ID To: UniProtKB   
   	Downloaded all "canonical & isoform" protein sequence data as Fasta.
   UniProt reported: 694 out of 696 UniProtKB AC/ID identifiers were successfully mapped to 603 UniProtKB IDs in the table below. Not mapped: O16797-4 and Q8MSV2-3.
   The un-mapped ones (O16797-4 and Q8MSV2-3) were substituted with another intermediate-sized isoform (see below).
   However, 7 were deleted from most recent UniProt database:
	A4V2J2
	A4V464
	E1JHL5
	Q1RMK2
	Q8IPM3
	Q9H552
	Q9VHL3
   These were obtained from this UniProt history manually - the most recent version in the history for each was used.
   A handful of accessions were real but the uncommon version and the downloaded UniProt sequences came with the more common accessions.
   These are handled below.

9. Ran get-mol-weight.py on protein fasta files to get molecular weights and amino acid lengths associated with each name.

10. At this point, because some of the accessions in our data are uncommon, we had to manually add the uncommon accessions to the name/molwt/length tables we generated.
	For all below, new entries were made where column#1 (the common accession that UniProt returned) is substituted with value in column#2 (the uncommon accession for same UniProt entry that was in our dataset).
	Thus, the final database of names/accessions with MW/Len information had matches for all accessions in our mass spec data.
	P78385 ---> Q6NT21
	Q61765 ---> A2A5Y0
	Q5KR48 ---> Q3SX28
	Q9W4M7 ---> E1JJD4
	E1JJD6 ---> E1JJD7
	Q8N1N4 ---> Q7RTT2
	P55830 ---> H9XVL9
	Q7KNS5 ---> A1ZAT7
	Q15323 ---> Q9UE12
	Q9NSB2 ---> Q6ISB0
	Q29443 ---> Q0IIK2
	Q7PL76 ---> Q7PL77
	Q9V463 ---> Q9VKL5
	Q8MSV2 ---> E1JID7
	Q8MSV2-2 ---> Q8MSV2-3 ## one of the unmapped ones (see #8 for more info)
	D1Z3A1 ---> Q9VZS4
	O16797-2 ---> O16797-4 ## one of the unmapped ones (see #8 for more info)


