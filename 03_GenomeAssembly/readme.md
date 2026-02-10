## Genome Assembly

Le tecnologie di sequenziamento attuali non permettono di leggere un intero cromosoma in un’unica volta. Al contrario, producono milioni di frammenti di DNA di lunghezza variabile che devono essere successivamente riassemblati.
L’**assembly genomico** è quindi il processo computazionale che ricostruisce la sequenza originale a partire da queste letture frammentate, sfruttando le regioni di sovrapposizione tra di esse.

### Livelli di assemblaggio

Durante il processo di ricostruzione genomica si distinguono diversi livelli strutturali:

* **Contig**: sequenza continua ottenuta assemblando reads sovrapposte, senza gap.
* **Scaffold**: insieme ordinato e orientato di contig, separati da gap (rappresentati da “N”), ricostruiti utilizzando informazioni di paired-end o long reads.
* **Cromosoma**: scaffold che rappresenta un cromosoma biologico completo; per raggiungere questo livello sono necessari dati addizionali come mappe genetiche o informazioni sulla struttura tridimensionale del DNA.

Nel presente lavoro l’assemblaggio è stato condotto fino al **livello di contig e scaffold**.

---

## 00 – Raw assembly

### Assemblaggio a livello di contig

L’assemblaggio iniziale è stato effettuato utilizzando **wtdbg2**, un assemblatore ottimizzato per dati long-read. È importante sottolineare che l’assemblaggio genomico è un processo fortemente empirico: non esiste un assemblatore universalmente “migliore”, ma la scelta dipende dal tipo di dati e dall’obiettivo dell’analisi.

```bash
ln -s /home/PERSONALE/mirko.martini3/Lab_CompGeno/00_practice/00_data/00_reads/SRR11672506.fastq.gz
conda activate assembly
wtdbg2 -x rs -g 227054799 -t 8 -i SRR11672506.fastq.gz -fo Anoste_raw
```

### Generazione della sequenza di consenso

L’output di wtdbg2 consiste in una struttura a grafo. Per ottenere una sequenza FASTA finale è necessario utilizzare **wtpoa-cns**, che applica un algoritmo di Partial Order Alignment (POA) per calcolare una sequenza di consenso a partire dal grafo di assemblaggio.

```bash
wtpoa-cns -t 8 -i Anoste_raw.ctg.lay.gz -fo Anoste_raw
```

Il file risultante (`Anoste_raw.fasta`) rappresenta l’assemblaggio grezzo a livello di contig.

---

## Valutazione della qualità dell’assemblaggio grezzo

La qualità dell’assemblaggio è stata valutata utilizzando tre metriche complementari:

### N50

Il valore **N50** fornisce una misura della continuità dell’assemblaggio. Rappresenta la lunghezza del contig più corto tra quelli che, sommati, coprono almeno il 50% della lunghezza totale del genoma. Valori elevati indicano un assemblaggio meno frammentato.

```bash
assembly-stats Anoste_raw.fasta > Anoste_raw.stats
```

### BUSCO

BUSCO valuta la **completezza biologica** dell’assemblaggio cercando geni altamente conservati appartenenti a un dataset tassonomico specifico.
Per questo lavoro è stato utilizzato il dataset `culicidae_odb12`, coerente con il gruppo tassonomico analizzato.

```bash
conda activate sequence
export NUMEXPR_MAX_THREADS=80
busco -m geno -l $BUSCO/culicidae_odb12 -c 6 -o Anoste_raw.busco -i Anoste_raw.fasta
```

### KAT

KAT confronta i k-mer presenti nelle reads originali con quelli dell’assemblaggio finale, permettendo di identificare regioni mancanti, duplicazioni o rumore di fondo.

```bash
conda activate kat
kat comp -t 8 -o Anoste_raw 'SRR11672503_1_paired.fastq SRR11672503_2_paired.fastq' Anoste_raw.fasta
```

---

## 01 – Polishing

### Allineamento delle reads

Il polishing richiede l’allineamento delle reads originali sull’assemblaggio grezzo. Questo passaggio consente di individuare errori di consenso e correggerli.
È stato utilizzato **minimap2** per mappare sia short reads che long reads.

```bash
minimap2 -ax sr --MD -t 6 Anoste_raw.fasta SRR11672503_1_paired.fastq SRR11672503_2_paired.fastq > Anoste_raw_sr.sam
samtools view -Sb Anoste_raw_sr.sam > Anoste_raw.bam
samtools sort -@6 -o Anoste_raw_sr_sorted.bam Anoste_raw.bam
samtools index Anoste_raw_sr_sorted.bam
```

```bash
minimap2 -ax map-pb --MD -t 6 Anoste_raw.fasta SRR11672506.fastq.gz > Anoste_raw_lr.sam
samtools view -Sb Anoste_raw_lr.sam > Anoste_raw_lr.bam
samtools sort -@6 -o Anoste_raw_lr_sorted.bam Anoste_raw_lr.bam
samtools index Anoste_raw_lr_sorted.bam
```

### Correzione con Hypo

Il polishing vero e proprio è stato effettuato con **Hypo**, che utilizza la copertura delle short reads per correggere errori residui nell’assemblaggio.

```bash
mosdepth -n --fast-mode --by 500 Anoste_raw_sr Anoste_raw_sr_sorted.bam
echo -e "$R1\n$R2" > Sr.Path
hypo -d Anoste_raw.fasta -r @Sr.Path -s 227m -c 136 \
     -b Anoste_raw_sr_sorted.bam -B Anoste_raw_lr_sorted.bam -t 6
mv hypo_Anoste_raw.fasta Anoste_pol.fasta
```

---

## Controllo qualità post-polishing

Dopo il polishing, l’intero controllo di qualità è stato ripetuto utilizzando le stesse metriche (N50, BUSCO e KAT), sostituendo l’input con `Anoste_pol.fasta`.

---

## 02 – Decontaminazione

Una volta ottenuto un assemblaggio polito, è fondamentale rimuovere eventuali sequenze contaminanti non appartenenti all’organismo target.
Le short reads sono state nuovamente mappate sul genoma polito e i contig sono stati annotati tassonomicamente tramite **blastn**.

```bash
blastn -query Anoste_pol.fasta -db <PATH/TO/nt/> \
-outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' \
-max_target_seqs 25 -max_hsps 1 -num_threads 25 -evalue 1e-25 \
-out Anoste_blast
```

### Identificazione dei contaminanti

Per l’identificazione visiva dei contaminanti è stato utilizzato **BlobTools**, che integra informazioni su GC content, copertura e tassonomia.

```bash
blobtools create -i Anoste_pol.fasta -b Anoste_pol_sorted.bam -t Anoste_blast.tsv -o Anoste_blob
blobtools plot -i Anoste_blob.blobDB.json -o Anoste
```

Sono stati mantenuti esclusivamente i contig assegnati al clade **Arthropoda**, ottenendo così un assemblaggio finale decontaminato a livello di contig.

---

## 03 – Scaffolding

Per migliorare ulteriormente la struttura dell’assemblaggio, è stato pianificato uno scaffolding basato su riferimento utilizzando **RagTag**, che sfrutta un genoma di una specie filogeneticamente vicina per ordinare e orientare i contig.

```bash
ragtag.py correct -t 20 <REFERENCE_GENOME> <DRAFT_GENOME>
ragtag.py scaffold -C -t 20 -o <OUTPUT_DIR> <REFERENCE_GENOME> <CORRECTED_DRAFTGENOME>
```

A causa di limitazioni temporali, questa fase non è stata eseguita direttamente e per le analisi successive è stato utilizzato l’assemblaggio fornito dal docente.
