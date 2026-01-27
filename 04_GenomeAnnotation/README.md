# Genome Annotation

L'annotazione genomica rappresenta la fase cruciale successiva all'assemblaggio. In questo stadio, la sequenza grezza viene analizzata per identificare e descrivere le entità funzionali presenti, quali sequenze codificanti (CDS), pseudogeni, trasposoni, ncRNA e tRNA.

L'obiettivo operativo è trasformare una sequenza anonima di nucleotidi in una mappa biologica, associando coordinate precise a funzioni specifiche. Il workflow si avvale di software come **SNAP** e **AUGUSTUS** e viene eseguito attraverso **round iterativi** che permettono di affinare progressivamente i modelli statistici, migliorando la sensibilità e l'accuratezza della predizione genica.



---

## 1. Annotazione degli Elementi Ripetuti

Poiché il focus del progetto è rivolto alle regioni codificanti, è indispensabile identificare e "mascherare" preventivamente i **repetitive elements**. Queste strutture, se non gestite, possono interferire con gli algoritmi di predizione generando falsi positivi.

> **Tool utilizzati:** `RepeatModeler` (identificazione *de novo*) e `RepeatMasker` (mascheramento).

---

## 2. Configurazione ed Esecuzione di MAKER

**MAKER** è il framework principale utilizzato per l'annotazione. La procedura inizia con la generazione dei file di configurazione necessari tramite il comando:

```bash
maker -CTL

```

In questo progetto, l'annotazione è stata guidata utilizzando il proteoma di *A. stephensis* come evidenza di omologia. Di seguito sono riportati i parametri chiave configurati nel file `maker_opts.ctl`:

```bash
#-----Input Genomico
genome= <GENOME_FILE> 
organism_type=eukaryotic 

#-----Evidenze di Omologia
protein= <PROTEOME_FILE> # Proteoma di A. stephensis
protein2genome=1 

#-----Repeat Masking
rmlib= <REPEAT_LIBRARY> 
softmask=1 

#-----Ab-initio Gene Prediction
snaphmm= <SNAP_HMM_FILE>
augustus_species= <SPECIES_MODEL>

```

Al termine dell'elaborazione, per riorganizzare la directory di output ed eliminare i file temporanei generati durante il calcolo distribuito, si utilizza:

```bash
maker -base <OUTPUT_PREFIX>

```

---

## 3. Consolidamento dei Risultati

Una volta completati i calcoli, i dati frammentati devono essere unificati per generare il file `.gff` finale (contenente le coordinate e le annotazioni) e i relativi file FASTA (sequenze nucleotidiche e amminoacidiche).

```bash
# Unificazione dei risultati
fasta_merge -d <DATASTORE_INDEX_FILE>
gff3_merge -d <DATASTORE_INDEX_FILE>

```

> Il file **.gff** finale è strutturato in colonne standardizzate che riportano tutte le caratteristiche delle sequenze identificate come codificanti.

---

## 4. Algoritmi di Predizione Genica

La pipeline integra due software di predizione che operano con logiche complementari:

* **[SNAP](https://github.com/KorfLab/SNAP.git) (Semi-HMM-based Nucleic Acid Parser):** Un programma di ricerca genica basato su modelli HMM, addestrabile su set di geni noti per ottimizzare la specificità dell'annotazione.
* **[AUGUSTUS](https://github.com/Gaius-Augustus/Augustus.git):** Un software di predizione che opera in modalità **Ab initio**.
> *Nota:* I modelli Ab initio non si basano sulla somiglianza con sequenze di database, ma sfruttano le proprietà statistiche intrinseche del DNA (composizione in basi, siti di splicing, ecc.) per predire la struttura genica.



---

## 5. Sintesi del Workflow Iterativo

Il processo di annotazione segue un approccio ciclico per massimizzare la qualità del risultato:

1. **Masking:** Identificazione e mascheramento delle ripetizioni genomiche.
2. **Evidence Alignment:** Allineamento del proteoma di riferimento per guidare la ricerca iniziale.
3. **Model Training:** Estrazione di un set iniziale di modelli genici ad alta confidenza per l'addestramento di SNAP e AUGUSTUS.
4. **Iterazione:** Esecuzione di round successivi per raffinare i parametri dei predittori.
5. **Final Consensus:** Generazione dell'annotazione finale e valutazione della qualità.
