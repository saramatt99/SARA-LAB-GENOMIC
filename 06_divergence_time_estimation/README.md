# Calibrazione Temporale dell'Albero

In questa fase, l'albero filogenetico generato precedentemente (nella directory `05_OG.Inference_Phylogenomic`) viene calibrato per ottenere un **Time Tree** che rappresenti i tempi di divergenza evolutiva.

### 1. Stima dei Tempi di Divergenza
Per ottenere i punti di calibrazione, è stata utilizzata la lista delle specie analizzate:

* *Anopheles gambiae*
* *Culex quinquefasciatus*
* *Aedes aegypti*
* *Aedes albopictus*
* *Sabethes cyaneus*
* *Anopheles stephensi*

> **Nota:** Durante il reperimento dei dati (es. tramite TimeTree/iTOL), per *Anopheles stephensi* non erano disponibili informazioni dirette (segnalato da un asterisco). Di conseguenza, è stata utilizzata come riferimento la specie affine *Anopheles karwari*. La divergenza interna al genere *Anopheles* non viene quindi considerata come punto di calibrazione in questo step.

### 2. Creazione del File di Calibrazione
Sulla base delle stime ottenute, è stato generato il file `calibration.txt`. Il formato prevede l'elenco dei taxa coinvolti nel nodo e il tempo di divergenza (in milioni di anni).

**Contenuto di `calibration.txt`:**
```text
Anogam,Anoara -1.5
Anogam,Anoara,Anofun -25
Anoste,Anoara,Anogam -20
Anogam,Anoara,Anofun,Anocul,Anoste -45

```

### 3. Esecuzione della Calibrazione (IQ-TREE)

Prima di procedere, è necessario linkare nella directory di lavoro corrente il file delle sequenze concatenate (generato nello step di trimming dei single copy orthologs).

```bash
# Link simbolico al file concatenato
ln -s ../05_OG.Inference_Phylogenomic/04_trimmed/00_single_complete/conc_species_tree .

```

Successivamente, si esegue **IQ-TREE** specificando il modello evolutivo (identificato nel file `.iqtree` precedente, qui `Q.INSECT+F+I+R3`) e passando il file di calibrazione.

```bash
iqtree -s conc_species_tree \
       --date calibration.txt \
       --date-tip 0 \
       -o Anoste,Anogam \
       -m Q.INSECT+F+I+R3 \
       -nt 13 \
       --prefix time_tree \
       --date-options "-u 1"

```

* `--date`: specifica il file con i vincoli temporali.
* `--date-tip 0`: assume che le sequenze attuali abbiano età 0 (presente).
* `-o`: definisce l'outgroup per il rooting.
* `--date-options "-u 1"`: imposta il limite superiore per la calibrazione.

```

```
