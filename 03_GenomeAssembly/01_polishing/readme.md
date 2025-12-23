# Polishing

Questa sezione descrive i comandi e i software utilizzati per migliorare la qualità dell'assemblaggio del genoma.

## Mappatura delle reads

### Short reads

Le short reads vengono mappate sul genoma assemblato utilizzando minimap2. Il codice è contenuto nel file `mapping.sh`.

```bash
#|assembly|
minimap2 -ax --MD -t 6 Anoste_raw.fasta SRR11672503_1_paired_fastq SRR11672503_1_paired_fastq > Anoste_raw_sr.sam
samtools view -Sb Anoste_raw_sr.sam > Anoste_raw_sr.bam
rm Anoste_raw_sr.sam
samtools sort -@6 -o Anoste_raw_sr_sorted.bam Anoste_raw_sr.bam
samtools index Anoste_raw_sr_sorted.bam
rm Anoste_raw_sr.bam
```

### Long reads

Analogamente, le long reads vengono mappate sul genoma assemblato.

```bash
#|assembly|
minimap2 -ax --MD -t 6 Anoste_raw.fasta SRR11672503_1_paired_fastq SRR11672503_1_paired_fastq > Anoste_raw_lr.sam
samtools view -Sb Anoste_raw_lr.sam > Anoste_raw_lr.bam
rm Anoste_raw_lr.sam
samtools sort -@6 -o Anoste_raw_lr_sorted.bam Anoste_raw_lr.bam
samtools index Anoste_raw_lr_sorted.bam
rm Anoste_raw_lr.bam
```

## Pulizia dell'assemblaggio

La sequenza del genoma viene migliorata grazie al confronto con le reads mappate precedentemente.

```bash
#|assembly|
echo e- "$R1\n$R2" > Sr.path
hypo -d Anoste_raw.fasta -r @Sr.path -s 227054799 -c 136 -b Anoste_raw_sr_sorted.bam -B Anoste_raw_lr_sorted.bam -t 6
```

## Controllo qualità del genoma pulito

Verifica delle statistiche relative al genoma post-pulizia. I metodi utilizzati sono gli stessi dell'assemblaggio raw.

### N50

Metodo per verificare la contiguità relativa dell'assemblaggio.

```bash
#|assembly|
assembly-stats Anoste_pol.fasta > Anoste_pol.stats
```

### BUSCO

Confronta l'assemblaggio con un dataset di geni del livello tassonomico di interesse, restituendo la percentuale di geni ritrovati.

```bash
#|sequence|
busco -m geno -l $BUSCO/culicidae_odb12 -c 8 -o Anoste_pol_busco -i Anoste_pol.fasta
```

### Spectra-cn (KAT)

Metodo per confrontare le frequenze dei kmer nelle reads originali con quelle del genoma assemblato.

```bash
#|kat|
kat comp -t 8 -o Anoste_pol 'SRR11672503_1_paired.fastq SRR11672503_2_paired.fastq' Anoste_pol.fast
```
