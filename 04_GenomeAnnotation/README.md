# Genome Annotation

L'annotazione rappresenta la fase cruciale che segue l'assemblaggio di un genoma. In questo stadio, la sequenza genomica grezza viene analizzata per identificare e descrivere le entità funzionali presenti, come sequenze codificanti (CDS), pseudogeni, trasposoni, ncRNA e tRNA. L'obiettivo operativo è trasformare una sequenza anonima di nucleotidi in una mappa biologica, associando coordinate precise a funzioni specifiche.

Il workflow viene solitamente eseguito in diversi round iterativi: ogni passaggio permette agli algoritmi di affinare i propri modelli statistici, migliorando progressivamente la sensibilità e la precisione nell'identificazione dei geni.

---

### 1. Annotazione degli Elementi Ripetuti

In questo progetto, il focus primario è rivolto alle regioni codificanti del genoma. È quindi indispensabile identificare e "mascherare" preventivamente i **repetitive elements** (elementi ripetuti), poiché le loro caratteristiche strutturali potrebbero trarre in inganno gli algoritmi di predizione genica.
* **Tool utilizzati:** Per la gestione di queste sequenze sono stati impiegati software come `RepeatModeler` (per l'identificazione *de novo*) e `RepeatMasker` (per il mascheramento).

---

### 2. Generazione dell'Annotazione con MAKER

Il cuore della pipeline di annotazione è **MAKER**. Per configurare correttamente l'analisi, il software richiede la generazione di file di controllo (`.ctl`) tramite il comando:

```bash
maker -CTL
Dopo aver configurato i file (integrando, in questo caso, il proteoma della specie A. stephensis come evidenza di omologia), si procede all'esecuzione del programma. Di seguito sono riportati i parametri principali impostati nel file di configurazione:

Plaintext
#-----Genome (Required)
genome= <GENOME_FILE> 
organism_type=eukaryotic 

#-----Protein Homology Evidence
protein= <PROTEOME_FILE> # Proteoma di A. stephensis

#-----Repeat Masking
rmlib= <REPEAT_LIBRARY> 
softmask=1 

#-----Gene Prediction
snaphmm= <SNAP_HMM_FILE>
augustus_species= <AUGUSTUS_MODEL>
protein2genome=1 
Al termine dell'elaborazione, per pulire la directory di lavoro dai file intermedi e organizzare l'output, si utilizza:

Bash
maker -base <OUTPUT_PREFIX>
3. Consolidamento dei Risultati
Per rendere i risultati fruibili, le informazioni frammentate generate da MAKER devono essere unificate in un unico file .gff (che contiene le coordinate e le annotazioni) e nei relativi file FASTA (sequenze proteiche e nucleotidiche).

Bash
# Unificazione dei dati estratti
fasta_merge -d <DATASTORE_INDEX_FILE>
gff3_merge -d <DATASTORE_INDEX_FILE>
4. Algoritmi di Predizione Genica
La pipeline integra due software specializzati per l'identificazione delle sequenze codificanti:

SNAP (Semi-HMM-based Nucleic Acid Parser): Un programma di ricerca genica adattabile che viene addestrato su set di geni noti per ottimizzare la qualità della predizione.

AUGUSTUS: Un software di predizione che può operare in modalità Ab initio.

Nota: I modelli Ab initio non si basano sulla somiglianza con sequenze già depositate in database, ma sfruttano esclusivamente le proprietà statistiche intrinseche del DNA (come la composizione in basi e i siti di splicing) per predire la struttura genica.

5. Sintesi del Workflow Iterativo
Il processo di annotazione può essere riassunto nei seguenti passaggi chiave:

Masking: Identificazione e mascheramento delle ripetizioni genomiche.

Evidence Alignment: Utilizzo di evidenze esterne (es. proteoma) per guidare la ricerca.

Model Training: Estrazione di un set iniziale di modelli genici ad alta confidenza.

Iterazione: Addestramento di SNAP e AUGUSTUS e ripetizione dei round per affinare la precisione.

Consensus & Evaluation: Creazione di un set di geni di consenso e valutazione finale dell'annotazione.
