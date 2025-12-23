## Kmer - Based Genome Survey

Il procedimento iniziale consente di ottenere informazioni sulla qualità delle short-reads ottenute mediante il sequenziamento dell'organismo *Anopheles stephensi*. Il workflow è suddiviso in tre fasi principali, riportate di seguito.

### 1. Controllo qualitativo

Analisi della qualità delle reads mediante FastQC.

```bash
#|assembly|
fastqc SRR11672503_1.fastq.gz SRR11672503_2.fastq.gz
```

### 2. Suddivisione delle sequenze paired e unpaired

Le short-reads vengono separate in reads paired e unpaired basandosi sulla complementarietà tra le letture sovrapposte.

```bash
#|assembly|
trimmomatic PE -threads 20 -phred33 \
  SRR11672503_1.fastq.gz SRR11672503_2.fastq.gz \
  SRR11672503_1_paired.fastq SRR11672503_1_unpaired.fastq \
  SRR11672503_2_paired.fastq SRR11672503_2_unpaired.fastq \
  ILLUMINACLIP:/opt/miniforge3/envs/assembly/share/trimmomatic-0.40-0/adapters/TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> stats_trimmomatic.log
```

### 3. Applicazione del metodo Kmer

Il metodo Kmer viene applicato alle reads paired risultanti dal passaggio precedente per valutare la distribuzione e la frequenza dei kmer.

```bash
#|kat|
kat hist -t 6 -m 27 -o Anoste_kmer27 SRR11672503_1_paired_fastqc.html SRR11672503_2_paired_fastqc.html
```
