# Code Book

mkdir ~/data/blast

#creates directory for BLAST database

gunzip proteomes/*.gz

#unzip proteomes in seperate directory

cat  proteomes/* > ~/data/blast/allprotein.fas

#puts all protein sequences into a single file

makeblastdb -in ~/data/blast/allprotein.fas -parse_seqids -dbtype prot

#creates BLAST database

ncbi-acc-download -F fasta -m protein XP_032219055.1

#dowloads amino acid sequence for Neamto Stella CAMSAP

blastp -db ~/data/blast/allprotein.fas -query XP_032219055.1.fa -outfmt 0 -max_hsps 1 > XP_032219055.blastp.typical.out\

#preforms  a blast search with the query protein

blastp -db ~/data/blast/allprotein.fas -query XP_032219055.1.fa -outfmt "6 sseqid pident length mismatch gapopen evalue bitscore pident stitle" -max_hsps 1 -out XP_032219055.blastp.detail.out

#creates a more detailed output of the same analysis

awk '{if ($6<0.00000000000001)print $1 }' XP_032219055.blastp.detail.out > XP_032219055.blastp.detail.filtered.out

#Filters Blast search for homlogues such that the search will only contain homolougues with an e value less than 1e-14

wc -l XP_032219055.blastp.detail.filtered.out

#counts the number of blast results

seqkit grep --pattern-file XP_032219055.blastp.detail.filtered.out ~/data/blast/allprotein.fas > XP_032219055.blastp.detail.filtered.fas

#aligns the gene family sequences.  At this point the homolougues are renamed to make them easily identifiable

muscle -in XP_032219055.blastp.detail.filtered.fas -out XP_032219055.blastp.detail.filtered.aligned.fas

#preforms a global multi sequence alignement 

t_coffee -other_pg seq_reformat -in XP_032219055.blastp.detail.filtered.aligned.fas -output sim

#provides statistics about the alignment

t_coffee -other_pg seq_reformat -in XP_032219055.blastp.detail.filtered.aligned.fas -action +rm_gap 50 -out allhomologs.aligned.r50.9055.f

#removes all columns from allignment that contain greater than 50% gapped residues

alv -kli --majority allhomologs.aligned.r50.f | less -RS

#view final alignment

iqtree -s allhomologs.aligned.r50.9055.f -nt 2

#creates phylogenetic tree from alignment

nw_display allhomologs.aligned.r50.9055.f.treefile

#view tree with arbitrary root

gotree draw png -w 1000 -i allhomologs.aligned.r50.9055.f.treefile -r -o allhomologs.aligned.r50.9055.f.treefile.png

#view unrooted tree. Creates png of unrooted tree

gotree reroot midpoint -i allhomologs.aligned.r50.9055.f.treefile -o allhomologs.aligned.r50.9055.f.midpoint.treefile

#roots tree at the midpoint

nw_display allhomologs.aligned.r50.9055.f.midpoint.treefile

#view tree in console

nw_display -s allhomologs.aligned.r50.9055.f.midpoint.treefile -w 1000 -b 'opacity:0' > allhomologs.aligned.r50.9055.f.midpoint.svg

#create svg of midpoint rooted tree
