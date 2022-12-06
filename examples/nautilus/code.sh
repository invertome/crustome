
#1 - Dowload from OrthoDB, orthologs for Drosophila melanogaster's nautilus gene (Dmel + Bmori + crustaceans) - Myogenic-determination protein

## FIRST ITERATION

#2 - Align reference sequences
mafft --dash --originalseqonly --genafpair --maxiterate 10000 --thread 14 nau_reference.fasta > reference.aligned.fasta

#3 - Build an HMM profile using previously aligned reference sequences
hmmbuild -n nautilus -o nau.1.log --cpu 14 nau.1.hmm reference.aligned.fasta

#4 - Index crustome aminoacid database 
#esl-sfetch --index crustome.pep 

#5 - Use HMM profile to search vs. CrusTome's transcriptomes
hmmsearch --cpu 14 --tblout nau.1.hits.tbl -E 1e-18 nau.1.hmm crustome.pep

#6 - Extract hits as full sequences from CrusTome
grep -v "^#" nau.1.hits.tbl | awk '{print $1}' | esl-sfetch -f crustome.pep - > nau.1.hits.fa

#7 - Clean fasta file 
python ~/bin/filter_faa.py nau.1.hits.fa > nau.1.hits.clean.fa

#8 - Remove sequences shorter than 120 aa 
seqtk seq -C -L 120 nau.1.hits.clean.fa > nau.1.hits.filtered.fa

## SECOND ITERATION

#9 - Align reference sequences
mafft --dash --originalseqonly --genafpair --maxiterate 10000 --thread 14 nau.1.hits.filtered.fa > nau.1.hits.filtered.aligned.fasta

#10 - Build an HMM profile using previously aligned reference sequences
hmmbuild -n nautilus -o nau.2.log --cpu 14 nau.2.hmm nau.1.hits.filtered.aligned.fasta

#11 - Use HMM profile to search vs. CrusTome's transcriptomes
hmmsearch --cpu 14 --tblout nau.2.hits.tbl -E 1e-18 nau.2.hmm crustome.pep

#12 - Extract hits as full sequences from CrusTome
grep -v "^#" nau.2.hits.tbl | awk '{print $1}' | esl-sfetch -f crustome.pep - > nau.2.hits.fa

#13 - Clean fasta file 
python ~/bin/filter_faa.py nau.2.hits.fa > nau.2.hits.clean.fa

#14 - Remove sequences shorter than 120 aa 
seqtk seq -C -L 120 nau.2.hits.clean.fa > nau.2.hits.filtered.fa

## THIRD ITERATION

#15 - Align reference sequences
mafft --dash --originalseqonly --genafpair --maxiterate 10000 --thread 14 nau.2.hits.filtered.fa > nau.2.hits.filtered.aligned.fasta

#16 - Build an HMM profile using previously aligned reference sequences
hmmbuild -n nautilus -o nau.3.log --cpu 14 nau.3.hmm nau.2.hits.filtered.aligned.fasta

#17 - Use HMM profile to search vs. CrusTome's transcriptomes
hmmsearch --cpu 14 --tblout nau.3.hits.tbl -E 1e-18 nau.3.hmm crustome.pep

#18 - Extract hits as full sequences from CrusTome
grep -v "^#" nau.3.hits.tbl | awk '{print $1}' | esl-sfetch -f crustome.pep - > nau.3.hits.fa

#19 - Clean fasta file 
python ~/bin/filter_faa.py nau.3.hits.fa > nau.3.hits.clean.fa

#20 - Remove sequences shorter than 120 aa 
seqtk seq -C -L 120 nau.3.hits.clean.fa > nau.3.hits.filtered.fa

## FOURTH ITERATION

#21 - Align reference sequences
mafft --dash --originalseqonly --genafpair --maxiterate 10000 --thread 14 nau.3.hits.filtered.fa > nau.3.hits.filtered.aligned.fasta

#22 - Build an HMM profile using previously aligned reference sequences
hmmbuild -n nautilus -o nau.4.log --cpu 14 nau.4.hmm nau.3.hits.filtered.aligned.fasta

#23 - Use HMM profile to search vs. CrusTome's transcriptomes
hmmsearch --cpu 14 --tblout nau.4.hits.tbl -E 1e-18 nau.4.hmm crustome.pep

#24 - Extract hits as full sequences from CrusTome
grep -v "^#" nau.4.hits.tbl | awk '{print $1}' | esl-sfetch -f crustome.pep - > nau.4.hits.fa

#25 - Clean fasta file 
python ~/bin/filter_faa.py nau.4.hits.fa > nau.4.hits.clean.fa

#26 - Remove sequences shorter than 120 aa 
seqtk seq -C -L 120 nau.4.hits.clean.fa > nau.4.hits.filtered.fa

## FIFTH ITERATION

#27 - Align reference sequences
mafft --dash --originalseqonly --genafpair --maxiterate 10000 --thread 14 nau.4.hits.filtered.fa > nau.4.hits.filtered.aligned.fasta

#28 - Build an HMM profile using previously aligned reference sequences
hmmbuild -n nautilus -o nau.5.log --cpu 14 nau.5.hmm nau.4.hits.filtered.aligned.fasta

#29 - Use HMM profile to search vs. CrusTome's transcriptomes
hmmsearch --cpu 14 --tblout nau.5.hits.tbl -E 1e-18 nau.5.hmm crustome.pep

#30 - Extract hits as full sequences from CrusTome
grep -v "^#" nau.5.hits.tbl | awk '{print $1}' | esl-sfetch -f crustome.pep - > nau.5.hits.fa

#31 - Clean fasta file 
python ~/bin/filter_faa.py nau.5.hits.fa > nau.5.hits.clean.fa

#32 - Remove sequences shorter than 120 aa 
seqtk seq -C -L 120 nau.5.hits.clean.fa > nau.5.hits.filtered.fa

## PHYLOGENETIC INFERENCE

#32 - Align final hit list and reference sequences
cat nau.5.hits.filtered.fa nau_reference.fasta > nau.all.filtered.fa
mafft --dash --originalseqonly --genafpair --maxiterate 10000 --thread 14 nau.all.filtered.fa > nau.all.aligned.fa

#33 - Trim alignment gaps with clipkit

cut -d ' ' -f1 nau.all.aligned.fa > nau.all.aligned.clean.fa
clipkit nau.all.aligned.clean.fa -m smart-gap


#34 - Infer phylogeny with iqtree2 for trimmed alignment

iqtree -s nau.all.aligned.clean.fa.clipkit -pre nau.trim -nt 14 -msub nuclear -bb 1000 -bnni 

#35 - Infer phylogeny with iqtree2 for untrimmed alignment

iqtree -s nau.all.aligned.clean.fa -pre nau.notrim -nt 14 -msub nuclear -bb 1000 -bnni

#36 - Prune long-branches from phylogenetic trees using treeshrink 
mkdir trees
mkdir trees/trim
cp nau.trim.contree trees/trim/input.tree
cp nau.all.aligned.clean.fa.clipkit trees/trim/input.fasta

mkdir trees
mkdir trees/notrim
cp nau.notrim.contree trees/notrim/input.tree
cp nau.all.aligned.clean.fa trees/notrim/input.fasta

run_treeshrink.py -i trees -q 0.05 --force

#37 - Re-align pruned trimmed alignment 

mafft --dash --originalseqonly --genafpair --maxiterate 10000 --thread 14 trees/trim/output.fasta > trees/trim/realigned.fasta


#38 - Infer trimmed phylogeny 

iqtree -s trees/trim/realigned.fasta -pre nau.trim.pruned -nt 14 -m JTT+F+R4 -nstop 250 -bb 10000 -bnni -abayes

#39 - Replace CrusTome IDs with species names

sed -f <(sed -E 's_(.+)\t(.+)_s/\1/\2/g_' dict) nau.trim.pruned.treefile > nau.trim.pruned.final.tree

#40 - Re-align pruned untrimmed alignment

mafft --dash --originalseqonly --genafpair --maxiterate 10000 --thread 14 trees/notrim/output.fasta > trees/notrim/realigned.fasta

#41 - Infer untrimmed phylogeny 

iqtree -s trees/notrim/realigned.fasta -pre nau.notrim.pruned -nt 14 -m JTT+F+R4 -nstop 250 -bb 10000 -bnni -abayes

#42 - Replace CrusTome IDs with species names

sed -f <(sed -E 's_(.+)\t(.+)_s/\1/\2/g_' dict) nau.notrim.pruned.treefile >  nau.notrim.pruned.final.tree

#END
