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

iqtree -s allhomologs.aligned.r50.9055.f -bb 1000 -nt 2 --prefix XP_032219055.r50.ufboot

#generates bootstrap support values

gotree reroot midpoint -i XP_032219055.r50.ufboot.treefile -o XP_032219055.r50.ufboot.midpoint.treefile

#creates midpoint rooted tree with bootstrap support

nw_display XP_032219055.r50.ufboot.midpoint.treefile

#view midpoint rooted tree with bootstrap support in console

nw_display -s XP_032219055.r50.ufboot.midpoint.treefile -w 1000 -b 'opacity:0' >XP_032219055.r50.ufboot.midpoint.treefile.svg

#creates svg image of midpoint rooted tree with bootstrap support

nano species.tre

#creare a species phylogeny for the 5 species being exaimined. paste "(((Homo_sapiens,Strongylocentrotus_purpuratus)Deuterostomia,Drosophila_melanogaster)Bilateria,(Nematostella_vectensis,Pocillopora_damicornis)Cnidaria)Eumetazoa;" as the text in the file

nano camsaptree.txt

#paste species.tre and XP_032219055.r50.ufboot.midpoint.treefile into this file

java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar -b camsaptree.txt --reconcile --speciestag prefix --savepng --treestats --events --phylogenomics

#reconciles gene tree with species tree

python ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g XP_032219055.r50.ufboot.midpoint.treefile.reconciled --include.species

#Generates a RecPhyloXML object to view the gene-within-species tree. push this file to a repository and examine at http://phylariane.univ-lyon1.fr/recphyloxml/recphylovisu.

less XP_032219055.r50.ufboot.midpoint.treefile.reconciled.events.txt

#view total dupications and losses for the reconcilled tree

java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar -b camsaptree.txt --root --speciestag prefix --savepng --treestats --events --phylogenomics

#reroot tree to minimize duplications and losses. Compare to reconcilled tree

nano camsaptree.txt

#replace gene tree with XP_032219055.r50.ufboot.midpoint.treefile.rooting.0.  This is the rerooted tree

java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar -b camsaptree.txt --rearrange --speciestag prefix  --savepng --treestats --events  --outputdir CAMSAPtreeroot_rearrange --edgeweights name --threshold 90 

#rearranges rerooted tree using a threshold of 90 to rearrange poorly supported branches. compare the re rooted tree to the rearranged tree in the new "CAMSAPtreeroot_rearrange" file.  If they are identical, (as they were when I ran this code) stop and continue after the command to run the topology test.  If not, note changes between trees and continue

python ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g XP_032219055.r50.ufboot.midpoint.treefile.rooting.0.rearrange.0 --include.species

#creates XML to view rearranged tree in species tree.  View file at http://phylariane.univ-lyon1.fr/recphyloxml/recphylovisu.

 java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar -g CAMSAPtreeroot_rearrange/XP_032219055.r50.ufboot.midpoint.treefile.rooting.0.rearrange.0   -s species.tre --reconcile --speciestag prefix  --treeoutput newick --nolosses
 
 #converts rearranged tree from Notung to Newick.  Make sure this is done in the main directory
 
gotree unroot -i XP_032219055.r50.ufboot.midpoint.treefile.rooting.0.rearrange.0.reconciled -o XP_032219055.r50.ufboot.unrooted.treefile.rearrange

#unroots rearranged tree

cat XP_032219055.r50.ufboot.treefile XP_032219055.r50.ufboot.unrooted.treefile.rearrange > XP_032219055.r50.alternativetrees

#concentrates unrooted optimal tree and unrooted rearranged tree into a single file

iqtree -s allhomologs.aligned.r50.9055.f -z XP_032219055.r50.alternativetrees -au -zb 10000 --prefix CAMSAP_altTrees -m LG+F+R5 -nt 2 -te XP_032219055.r50.ufboot.treefile

#runs topology test.  If the p value for tree two (the rearranged tree) is less than 0.05 it is signifigantly worse than tree one and should be rejected.  If it is greater than 0.05, we cannot reject the rearranged tree and it must be used as part of the confidence set.
