## Ambienti di calcolo utilizzati

Durante il progetto sono stati utilizzati diversi **ambienti Conda**, ciascuno dedicato a una specifica fase del workflow di genomica comparata. L’uso di ambienti separati ha permesso di evitare conflitti tra dipendenze e di mantenere le analisi riproducibili.

I comandi base per la gestione degli ambienti sono:

```bash
conda activate <nome_ambiente>
conda deactivate
```

Di seguito sono riportati gli ambienti effettivamente utilizzati e i principali software installati.

---

### Ambiente `assembly`

Usato per le fasi di **controllo qualità delle reads**, **assemblaggio genomico** e **annotazione strutturale**.

* python 3.13.5
* fastqc 0.12.1
* trimmomatic 0.40
* minimap2 2.28
* samtools 1.21
* mosdepth 0.3.10
* assembly-stats 1.0.1
* hypo 1.0.3
* augustus 3.1
* maker 3.01.04
* blast 2.16.0
* diamond 2.1.10

---

### Ambiente `kat`

Utilizzato nella fase di **genome survey pre-assemblaggio**, per l’analisi dei K-mer.

* kat 2.4.2

---

### Ambiente `sequence`

Ambiente dedicato alla **gestione e valutazione delle sequenze**, inclusa la preparazione dei dataset per le analisi filogenomiche.

* agat 1.4.1
* blobtools 1.1.1
* bmge 1.12
* BUSCO 6.0.0
* mafft 7.526
* ncbi-datasets 18.3.1

---

### Ambiente `tree`

Usato per le **analisi filogenetiche e filogenomiche** e per lo studio dell’evoluzione delle famiglie geniche.

* orthofinder 2.5.5
* disco 1.4.1
* iqtree 3.0.1
* cafe 5.1.0
* mafft 7.526

---

### Ambiente `base`

Utilizzato per le **analisi statistiche e funzionali in R**, in particolare per l’analisi di arricchimento GO.

* R 4.5.1
