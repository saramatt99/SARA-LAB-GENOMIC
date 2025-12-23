# Raw Assembly

In questa fase si procede con l'assemblaggio del genoma a livello di contig (contig-level). Un contig rappresenta una sequenza contigua di nucleotidi ottenuta mediante la sovrapposizione e l'allineamento di molteplici reads (short e long reads).

## Assemblaggio contig-level

Mediante l'allineamento e la sovrapposizione delle long reads si ottiene una sequenza univoca del genoma a livello di contig.

```bash
wtdbg2 -x rs -g 227054799 -t 8 -i SRR11672506.fastq.gz -fo Anoste_raw
```

## Creazione del consensus

I contig ottenuti vengono riordinati per creare la sequenza di consensus, ovvero l'insieme di contig correttamente disposti.

```bash
wtpoa-cns -t 7 -i Anoste_raw.ctg.lay.gz -fo Anoste_raw
```

---

## Controllo qualitativo

Controllo della qualità della sequenza di consensus ottenuta tramite tre metodi differenti: N50, BUSCO e Spectra-cn (KAT).

### N50

Metodo utilizzato per verificare la contiguità relativa dell'assemblaggio. Indica il numero di contig che rappresentano almeno il 50% del genoma assemblato.

```bash
Assembly-stats Anoste_raw > Anoste_raw.stats
```

### BUSCO

Metodo che confronta l'assemblaggio con un dataset di geni del livello tassonomico di interesse, restituendo una percentuale delle tipologie di geni ritrovati.

```bash
export NUMEXPR_MAX_THREADS=80
busco -m geno -l BUSCO/culicidae_odb12 -c 6 -o Anoste_raw_busco -i Anoste_raw.fasta
```

### Spectra-cn (KAT)

Metodo che confronta le frequenze dei kmer nelle reads originali con quelle nel genoma assemblato per stimare la coverage.

```bash
kat comp -t 8 -o Anoste_raw 'SRR11672503_1_paired.fastq SRR11672503_2_paired.fastq' A
```
