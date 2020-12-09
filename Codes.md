# 16S rRNA analysis

## 16S rRNA extraction and assignment from metagenomic data

```
kraken2-build --db /software/kraken2_db/silva_db --special silva
kraken2 --db /software/kraken2_db/silva_db/16S_SILVA138_k2db  --paired --classified-out cseqs#.fq 1.clean.fq 2.clean.fq
kraken2 --db /software/kraken2_db/silva_db/16S_SILVA138_k2db  cseqs_1.fq cseqs_2.fq --output 16S_rRNA_copy.selection.fasta.out.txt --use-names --report sample_16S.txt --threads 32
```
