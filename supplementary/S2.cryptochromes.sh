## Code for the example analysis of cryptochromes and photolyases using the CrusTome database.

## CRYPTOCHROMES
#Initial search with Insect reference sequences
blastp -query cryptochromes.fasta -db /crustome/crustome_aa -num_threads 14 -max_target_seqs 500 -evalue 1e-120 -outfmt "6 qseqid sseqid evalue bitscore pident nident qlen slen qstart qend sstart send length mismatch gapopen" -out hits.tsv

#List of hit IDS
cut -f2 hits.tsv > prehits.list
sed '$!N; /^\(.*\)\n\1$/!P; D' prehits.list > hits.list

#Extract sequences from BLAST db
blastdbcmd -db /crustome/crustome_aa -dbtype prot -entry_batch hits.list -out hits.fasta -outfmt %f

#Second search using 1st hits as input
blastp -query hits.fasta -db /crustome/crustome_aa -num_threads 14 -max_target_seqs 100 -evalue 1e-120 -outfmt "6 qseqid sseqid evalue bitscore pident nident qlen slen qstart qend sstart send length mismatch gapopen" -out hits2.tsv

#Populate List of hit IDs from blast results, then remove duplicates with awk
cut -f2 hits2.tsv > prehits2.list
awk '!seen[$0]++' prehits2.list > hits2.list

#Extract sequences from BLAST db
blastdbcmd -db /crustome/crustome_aa -dbtype prot -entry_batch hits2.list -out hits2.fasta -outfmt %f

#Concatenate with original insect reference seqs
cat hits2.fasta cryptochromes.fasta > input.fa

#Align sequences
mafft --dash --originalseqonly --genafpair --maxiterate 10000 --thread 14 input.fa > alignment.fa

#Trim alignment gaps with clipkit
cut -d ' ' -f1 alignment.fa > aligned.fa
clipkit aligned.fa -m smart-gap

#Infer phylogeny with iqtree2
iqtree -s aligned.fa.clipkit -pre cryptotrim -nt 14 -m TESTNEW -msub nuclear -nstop 250 -bb 10000 -bnni -abayes
iqtree -s aligned.fa -pre crypto -nt 14 -m TESTNEW -msub nuclear -nstop 250 -bb 10000 -bnni -abayes

#TreeShrink
mkdir trees
mkdir trees/crypto
cp crypto.contree trees/crypto/input.tree
cp aligned.fa trees/crypto/input.fasta
mkdir trees/cryptotrim
cp cryptotrim.contree trees/cryptotrim/input.tree
cp aligned.fa.clipkit trees/cryptotrim/input.fasta
run_treeshrink.py -i trees -q 0.05

#Infer phylogeny with iqtree2
mafft --dash --originalseqonly --genafpair --maxiterate 10000 --thread 14 trees/cryptotrim/output.fasta > trees/cryptotrim/realigned.fasta
iqtree -s trees/cryptotrim/realigned.fasta -pre cryptotrim3 -nt 14 -mset LG+R10 -nstop 250 -bb 10000 -bnni -abayes

#Replace IDs with names in final treefile
sed -f <(sed -E 's_(.+)\t(.+)_s/\1/\2/g_' dict) cryptotrim3.treefile > cryptotrim3.final.tree
