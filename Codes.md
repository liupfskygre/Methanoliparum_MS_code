# Similarity of 47 MAGs belonging to ‘Ca. Methanoliparia’

## 16S rRNA gene sequence identity calculation

```
makeblastdb -in methanoliparum_16SrRNA.fa -dbtype nucl
blastn -db methanoliparum_16SrRNA.fa -query MAGs_16S.fasta -out MAGs_16S.txt -evalue 1e-5 -outfmt 6 -max_target_seqs 1 -num_threads 32
```

## AAI calculation

```
comparem aai_wf . Result --cpus 32 --file_ext fa
```

# Identification of genes coding for enzymes with alkane activation potentials

```
perl /scripts/SeqTools/limit2Length.pl -f assembled_scaffold.fasta -len 500 -o scaffold_500.fasta 
prodigal -f gff -i scaffold_500.fasta -o scaffold_500.fasta.gff -p meta -a scaffold_500.faa -d scaffold_500.fna
hmmsearch mcrA.hmm scaffold_500.faa > scaffold_500_mcrA.out
perl extract_fasta_from_list.pl -f scaffold_500.faa -l scaffold_500_mcrA.list -o mcrA.faa (used for phylogenetic analysis)
perl extract_fasta_from_list.pl -f scaffold_500.fna -l scaffold_500_mcrA.list -o mcrA.fna (used for relative abundance analysis)
/software/cdhit-4.6.2/cd-hit -i mcrA.fna -o ref_mcrA.fna -c 0.99 -d 20 -T 32
```

# Phylogenomic analyses

## Using concatenated ribosomal proteins (16 r-proteins)

```
checkm lineage_wf -x .fa /bins checkm_result -f checkm_result.txt -t 30 --pplacer_threads 30
perl /script/find_ribosome_faa_in_checkM2.0.pl -f checkm_result
cat Ribosomal*.txt > all_ribosomal.fas
blastp -query all_ribosomal.fas -db /new_tree_of_life/find_ribosomal/all_ribosomal.fas.clean -out all_ribosonmal_blast_new_tree.txt -outfmt 6 -evalue 1e-5 -max_target_seqs 1 -num_threads 30
mkdir ribosomal_for_tree/
perl /script/pick_ribosomal_with_blast_v2.0.pl
Input the fasta file: all_ribosomal.fas
Input the blast result:all_ribosonmal_blast_new_tree.txt
Input the out file: ribosomal_for_tree
cd ribosomal_for_tree/
mkdir aligned_ribosomals/
muscle -in ribosomal_for_tree_Ribosomal_L*(S*).fas -out aligned_ribosomals/ribosomal_for_tree_Ribosomal_L*(S*).muscle.fas
cd aligned_ribosomals/
trimal -in ribosomal_for_tree_Ribosomal_L*(S*).muscle.fas -out ribosomal_for_tree_Ribosomal_L*(S*).muscle.trim.fas -gt 0.95
perl /script/combine_aligned_ribosomal_v2.0.pl
Input alignment folder:./
/tools/iqtree/bin/iqtree-omp -s ref_ribosomal_combine.phy -nt 30 -m WAG -bb 1000
```

## Phylogenetic trees for AcrA/McrA, AssA related, BcrB, and BcrC protein sequence

```
makeblastdb -in mcrA_reference_protein_seqs.faa -dbtype prot
blastp -query bin.fa.faa -db mcrA_reference_protein_seqs.faa -out bin.fa.txt -outfmt 6 -evalue 1e-5 -max_target_seqs 1 -num_threads 30
muscle -in combined_mcrA.faa -out aligned_combined_mcrA.faa
/tools/iqtree/bin/iqtree-omp -s aligned_combined_mcrA.faa -nt 32 -m WAG -bb 1000
```

# Evaluation of the relative abundance and activity of ‘Ca. Methanoliparum’

## 16S rRNA gene sequences extraction and assignment from metagenomic data

```
kraken2-build --db /software/kraken2_db/silva_db --special silva
kraken2 --db /software/kraken2_db/silva_db/16S_SILVA138_k2db  --paired --classified-out cseqs.fq 1.clean.fq 2.clean.fq
kraken2 --db /software/kraken2_db/silva_db/16S_SILVA138_k2db  cseqs_1.fq cseqs_2.fq --output 16S_rRNA_copy.selection.fasta.out.txt --use-names --report sample_16S.txt --threads 32 
```

## Relative abundance and expression activity of dereplicated MAGs

```
dRep dereplicate outout_directory -g /bins/*.fa -p 32 -comp 50 -con 10 -sa 0.97 -genomeInfo ./checkM_gtdbtk_dRep.csv
/tools/bowtie2-2.3.4.1-linux-x86_64/bowtie2-build dRepMAGs.fa scaffold.fa
/tools/bowtie2-2.3.4.1-linux-x86_64/bowtie2 -x scaffold.fa -1 1.clean.fq -2 2.clean.fq -S all.sam -p 128
/tools/samtools-1.3.1/samtools view -bS all.sam > all.bam -@ 128
/tools/samtools-1.3.1/samtools sort all.bam -o all_sort.bam -@ 128
/tools/samtools-1.3.1/samtools index all_sort.bam
perl /scripts/SeqTools/length+GC.pl -f dRepMAGs.fa -len | sed 's/\s.*\t/\t/' | sed 's/\t/\t0\t/' > scaffold.bed
bedtools coverage -abam all_sort.bam -b scaffold.bed -counts > scaffold.bamstat
```

## Expression activity of all annoted genes of four ‘Ca. Methanoliparum’ representative MAGs

```
prodigal -f gff -i four_rep_MAGs.fasta -o four_rep_MAGs.fasta.gff -p meta -a four_rep_MAGs.faa -d four_rep_MAGs.fna
bwa index four_rep_MAGs.fna
bwa mem -t 32 four_rep_MAGs.fna 1.clean.fq 2.clean.fq > scaffold.sam 2> scaffold.log
/tools/samtools-1.3.1/samtools view -bS scaffold.sam > all.bam -@ 32
/tools/samtools-1.3.1/samtools sort all.bam -o all_sort.bam -@ 32
/tools/samtools-1.3.1/samtools index all_sort.bam
perl /scripts/SeqTools/length+GC.pl -f four_rep_MAGs.fna -len | sed 's/\s.*\t/\t/' | sed 's/\t/\t0\t/' > scaffold.bed
bedtools coverage -abam all_sort.bam -b scaffold.bed -counts > scaffold.bamstat
```

## Expression activity of all annoted genes of assembled scaffold

```
perl /scripts/SeqTools/limit2Length.pl -f assembled_scaffold.fasta -len 500 -o scaffold_500.fasta 
prodigal -f gff -i scaffold_500.fasta -o scaffold_500.fasta.gff -p meta -a scaffold_500.faa -d scaffold_500.fna
bwa index scaffold_500.fna
bwa mem -t 32 scaffold_500.fna 1.clean.fq 2.clean.fq > scaffold.sam 2> scaffold.log
/tools/samtools-1.3.1/samtools view -bS scaffold.sam > all.bam -@ 32
/tools/samtools-1.3.1/samtools sort all.bam -o all_sort.bam -@ 32
/tools/samtools-1.3.1/samtools index all_sort.bam
perl /scripts/SeqTools/length+GC.pl -f scaffold_500.fna -len | sed 's/\s.*\t/\t/' | sed 's/\t/\t0\t/' > scaffold.bed
bedtools coverage -abam all_sort.bam -b scaffold.bed -counts > scaffold.bamstat
```

## Expression activity of AcrA/McrA, AssA-related and BcrB/BcrC

```
bwa index ref_mcrA.fna
bwa mem -t 32 ref_mcrA.fna 1.clean.fq 2.clean.fq > scaffold.sam 2> scaffold.log
/tools/samtools-1.3.1/samtools view -bS scaffold.sam > all.bam -@ 32
/tools/samtools-1.3.1/samtools sort all.bam -o all_sort.bam -@ 32
/tools/samtools-1.3.1/samtools index all_sort.bam
perl /scripts/SeqTools/length+GC.pl -f ref_mcrA.fna -len | sed 's/\s.*\t/\t/' | sed 's/\t/\t0\t/' > scaffold.bed
bedtools coverage -abam all_sort.bam -b scaffold.bed -counts > scaffold.bamstat
```

*notes*
```
all in-house scripts used are avaiable upon request.
contacts: Cuijing Zhang, zhangcj@szu.edu.cn
```

