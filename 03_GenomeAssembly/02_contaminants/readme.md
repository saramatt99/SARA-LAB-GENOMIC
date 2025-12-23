# Decontamination

La fase di **decontaminazione del genoma** consiste in una serie di procedure mirate a identificare ed eliminare sequenze geniche non appartenenti all'organismo di interesse (*Anopheles stephensi*). Di seguito sono riportati i comandi e i software utilizzati in questa fase.

## Ri-mappatura delle reads

Le short reads vengono rimappate sul genoma post-polishing per preparare il dataset alla fase di decontaminazione.
Il codice utilizzato si trova nel file `mapping_contaminants.sh`.

```bash
#|assembly|
minimap2 -ax sr --MD -t 8 Anoste_pol.fasta SRR11672503_1_paired_fastq SRR11672503_2_paired_fastq | samtools view -Sb - > Anoste_pol_sr.bam
samtools sort -@8 -o Anoste_pol_sr_sorted.bam Anoste_pol_sr.bam
samtools index Anoste_pol_sr_sorted.bam
rm Anoste_pol_sr.bam
```

## Annotazione tassonomica dei contigs

Ogni contig del genoma viene associato a un gruppo tassonomico grazie al confronto con la banca dati NCBI.

```bash
#|assembly|
blastn -query Anoste_pol.fasta -db <PATH/TO/nt/> \
-outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
-max_target_seqs 25 -max_hsps 1 -num_threads 25 -evalue 1e-25 -out Anoste_blast.tsv
```

## Individuazione dei contigs contaminati

Il programma **Blobtools** viene utilizzato per identificare contigs che non appartengono al genoma di *A. stephensi*.

```bash
#|sequence|
blobtools create -i Anoste_pol.fasta -b Anoste_pol_sr_sorted.bam -t Anoste_blast.tsv -o Anoste_blob
blobtools view -i Anoste_blob.blobDB.json -o Anoste
blobtools plot -i Anoste_blob.blobDB.json -o Anoste
```

## Salvataggio dei contigs non contaminati

Viene generato un file contenente solo i contigs associati a *A. stephensi* (identificati come "Arthropoda").

```bash
grep "Arthropoda" Anoste.Anoste_blob.blobDB.table.txt > contig_arthropoda.tsv
```

## Estrazione dei contigs decontaminati in formato fasta

Si associa il file dei contigs non contaminati con le sequenze genomiche per ottenere il genoma decontaminato finale.

```bash
awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' Anoste_pol.fasta \
| grep -w -Ff <(cut -f1 contig_arthropoda.tsv) - \
| tr "\t" "\n" > Anoste_decont
```
