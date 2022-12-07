### REPLACE "project" with project name
### REPLACE "reference.fasta" with your input fasta file with reference sequences
### REPLACE "/path/to/crustome/" with the path to where you extracted CrusTome's db files 
### ADJUST evalues according to genes of interest
### ADJUST "-num_threads" "--thread" and "-nt" according to computational resources available
### ADJUST model for iqtree2 according to model-testing results
### EXECUTE commands either line by line (recommended for beginners), as a script, or via a job scheduling system if on an HPC (e.g., SLURM, LSF).

### 1 - Initial search with reference sequences
blastp -query reference.fasta -db /path/to/crustome/crustome_aa -num_threads 5 -max_target_seqs 1000 -evalue 1e-96 -outfmt "6 qseqid sseqid evalue bitscore pident nident qlen slen qstart qend sstart send length mismatch gapopen" -out hits.tsv

### 2 - List of hit IDS. Remove duplicates.
cut -f2 hits.tsv > prehits.list
awk '!seen[$0]++' prehits.list > hits.list

### 3 - Extract sequences from BLAST db
blastdbcmd -db /path/to/crustome/crustome_aa -dbtype prot -entry_batch hits.list -out hits.fasta -outfmt %f

### 4 - Second search using 1st hits as input
blastp -query hits.fasta -db /path/to/crustome/crustome_aa -num_threads 5 -max_target_seqs 1000 -evalue 1e-96 -outfmt "6 qseqid sseqid evalue bitscore pident nident qlen slen qstart qend sstart send length mismatch gapopen" -out hits2.tsv

### 5 - Populate List of hit IDs from blast results, then remove duplicates with awk
cut -f2 hits2.tsv > prehits2.list
awk '!seen[$0]++' prehits2.list > hits2.list

### 6 - Extract sequences from BLAST db
blastdbcmd -db /path/to/crustome/crustome_aa -dbtype prot -entry_batch hits2.list -out hits2.fasta -outfmt %f

### 7 - Third search using 2nd hits as input

blastp -query hits2.fasta -db /path/to/crustome/crustome_aa -num_threads 5 -max_target_seqs 1000 -evalue 1e-96 -outfmt "6 qseqid sseqid evalue bitscore pident nident qlen slen qstart qend sstart send length mismatch gapopen" -out hits3.tsv

### 8 - Populate List of hit IDs from blast results, then remove duplicates with awk

cut -f2 hits3.tsv > prehits3.list
awk '!seen[$0]++' prehits3.list > hits3.list

### 9 - Extract sequences from BLAST db

blastdbcmd -db /path/to/crustome/crustome_aa -dbtype prot -entry_batch hits3.list -out hits3.fasta -outfmt %f

### 10 - Fourth search using 3rd hits as input

blastp -query hits3.fasta -db /path/to/crustome/crustome_aa -num_threads 5 -max_target_seqs 1000 -evalue 1e-96 -outfmt "6 qseqid sseqid evalue bitscore pident nident qlen slen qstart qend sstart send length mismatch gapopen" -out hits4.tsv

### 11 - Populate List of hit IDs from blast results, then remove duplicates with awk

cut -f2 hits4.tsv > prehits4.list
cat prehits.list prehits2.list prehits3.list prehits4.list > prehits_all.list
awk '!seen[$0]++' prehits_all.list > hits_all.list

### 12 - Extract sequences from BLAST db

blastdbcmd -db /path/to/crustome/crustome_aa -dbtype prot -entry_batch hits_all.list -out hits_all.fasta -outfmt %f

### 13 - Concatenate with original reference seqs

cat hits_all.fasta reference.fasta > input.fa

### 14 - Align sequences
mafft --dash --originalseqonly --genafpair --maxiterate 10000 --thread 10 input.fa > alignment.fa

### 15 - Trim alignment gaps with clipkit
cut -d ' ' -f1 alignment.fa > aligned.fa
clipkit aligned.fa -m smart-gap

### 16 - Infer phylogeny with iqtree2
### Suggestion: Use different clipkit settings and infer a phylogeny for each trimmed alignment version
iqtree -s aligned.fa.clipkit -pre project.trimmed -nt 10 -m TESTNEW -msub nuclear -bb 1000 -bnni -abayes
iqtree -s aligned.fa -pre project -nt 10 -m TESTNEW -msub nuclear -bb 1000 -bnni -abayes

### 17 - TreeShrink
mkdir trees
mkdir trees/project
cp project.contree trees/project/input.tree
cp aligned.fa trees/project/input.fasta
mkdir trees/project.trimmed
cp project.trimmed.contree trees/project.trimmed/input.tree
cp aligned.fa.clipkit trees/project.trimmed/input.fasta
run_treeshrink.py -i trees -q 0.05

### 18 - Infer phylogeny with iqtree2
### REPLACE "MODEL" by model given on first tree-building run
mafft --dash --originalseqonly --genafpair --maxiterate 10000 --thread 10 trees/project.trimmed/output.fasta > trees/project.trimmed/realigned.fasta
iqtree -s trees/project.trimmed/realigned.fasta -pre project.trimmed.final -nt 10 -mset MODEL -nstop 250 -bb 10000 -bnni -abayes

### 19 - Replace IDs with names in final treefile
### "dict" file should be copied to the folder where this is run 
cp /path/to/crustome/dict ./
sed -f <(sed -E 's_(.+)\t(.+)_s/\1/\2/g_' dict) project.trimmed.final.treefile > project.trimmed.final.names.treefile
