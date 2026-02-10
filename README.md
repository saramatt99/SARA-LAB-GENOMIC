## Comparative Genomics

Questa repository GitHub raccoglie il workflow sviluppato durante il corso di **Genomica Comparata** del corso di laurea magistrale in *Biodiversity and Evolution* (Università di Bologna).
Il progetto ha avuto come obiettivo l’applicazione di approcci di genomica e filogenomica per analizzare **le relazioni evolutive all’interno del genere *Anopheles*** e per valutare se il **complesso *Anopheles gambiae*** mostri un pattern evolutivo distinto rispetto ad altre linee del genere.

In particolare, l’analisi mira a integrare **ricostruzione filogenetica** ed **evoluzione delle famiglie geniche**, al fine di comprendere se la posizione del complesso *gambiae* rifletta una storia evolutiva peculiare, riconoscibile sia a livello di divergenza filogenetica sia di dinamica genica (espansioni e contrazioni).

Il progetto è strutturato in due componenti principali.
La prima riguarda l’assemblaggio e l’annotazione del genoma di **Anopheles stephensi**, utilizzato come riferimento.
La seconda consiste in uno studio comparativo che coinvolge **cinque specie del genere *Anopheles*** (*Anogam, Anoalb, Anoste, Anoara, Anofun*), selezionate in modo da rappresentare differenti complessi e aree biogeografiche.

Ogni cartella della repository documenta un passaggio specifico dell’analisi, seguendo un flusso logico che va dalla preparazione dei dati grezzi fino all’interpretazione evolutiva dei risultati.

---

### 00_data

Questa sezione descrive la costruzione del dataset genomico comparativo.
I genomi di *Anopheles stephensi* e delle altre cinque specie di *Anopheles* sono stati scaricati da **NCBI** e sottoposti a un processo di curazione volto a garantire la comparabilità tra i dataset.

Le operazioni principali includono:

* selezione delle isoforme proteiche più lunghe tramite **AGAT**;
* rimozione di pseudogeni e sequenze contenenti stop codon interni;
* uniformazione degli header delle sequenze, necessaria per le analisi di ortologia e filogenomica successive.

---

### 01_ComputerEnvs

Questa cartella documenta l’infrastruttura computazionale utilizzata.
Sono riportati gli ambienti **Conda** creati per isolare le dipendenze software delle diverse fasi del workflow, assicurando riproducibilità e stabilità delle analisi.

---

### 02_KmerBased_GenomeSurvey

In questa fase è stata effettuata una valutazione preliminare della qualità dei dati di sequenziamento di *Anopheles stephensi*.
Le reads sono state analizzate con **FastQC** e filtrate con **Trimmomatic** per rimuovere adattatori e regioni a bassa qualità.
Successivamente, un’analisi basata sui **K-mer** (KAT e GenomeScope) ha permesso di stimare parametri genomici di base, come dimensione del genoma ed eterozigosità, prima dell’assemblaggio.

---

### 03_GenomeAssembly

Questa sezione descrive il processo di assemblaggio del genoma di *Anopheles stephensi*.
L’assemblaggio è stato valutato iterativamente utilizzando metriche strutturali (N50) e biologiche (**BUSCO**, **KAT**) per verificarne completezza e qualità.
Sono stati inoltre applicati step di polishing, decontaminazione tassonomica e scaffolding finale per ottenere un genoma di riferimento affidabile.

---

### 04_GenomeAnnotation

Qui viene illustrato il workflow di annotazione strutturale del genoma assemblato.
L’annotazione è stata eseguita tramite **MAKER**, iniziando con una prima fase basata su evidenze di omologia proteica, seguita dall’addestramento dei predittori *ab initio*.
Una seconda esecuzione di MAKER ha prodotto l’annotazione finale, successivamente validata tramite statistiche generate con **AGAT**.

---

### 05_OG.Inference_Phylogenomic

Questa sezione descrive il cuore filogenomico dell’analisi.
I geni delle sei specie di *Anopheles* sono stati raggruppati in ortogruppi utilizzando **OrthoFinder**, distinguendo ortologhi e paraloghi.
Le famiglie multi-copia sono state ulteriormente suddivise in sottogruppi strettamente ortologhi mediante **DISCO**.

L’albero delle specie è stato ricostruito a partire da ortologhi single-copy, utilizzando allineamenti multipli rifiniti e un approccio supermatrix con **IQ-TREE**.
L’albero è stato radicato utilizzando *Anopheles culicifacies* come outgroup, permettendo una corretta polarizzazione dei caratteri evolutivi.

---

### 06_DivergenceTime_Estimation

In questa fase è stata stimata la tempistica delle divergenze evolutive tra le specie analizzate.
È stato costruito un albero ultrametrico tramite **IQ-TREE** e l’algoritmo **LSD2 (Least Square Dating)**, utilizzando calibrazioni temporali derivate da **TimeTree.org**.
Questa analisi consente di contestualizzare la posizione del complesso *Anopheles gambiae* nel tempo evolutivo del genere.

---

### 07_GeneFamilies_Evolution

Questa sezione è dedicata all’analisi dell’evoluzione delle famiglie geniche tramite **CAFE**.
Sono stati applicati modelli stocastici di nascita-morte per stimare variazioni nel numero di copie geniche lungo i rami dell’albero filogenetico.

Attraverso un processo di selezione dei modelli basato su **AIC** e **BIC**, è stato identificato il modello ottimale per individuare famiglie geniche significativamente espanse o contratte.
Questa analisi permette di valutare se il complesso *Anopheles gambiae* mostri un pattern di evoluzione delle famiglie geniche distinto rispetto alle altre linee del genere.

---

### 09_GeneAnnotation_functional_enrichment

Questa cartella descrive l’analisi di annotazione funzionale e di arricchimento.
Per ciascun ortogruppo è stata selezionata la proteina rappresentativa più lunga, utilizzata per ottenere annotazioni funzionali (GO terms).

Le famiglie geniche risultate significative dall’analisi CAFE sono state sottoposte a **Gene Ontology enrichment** mediante **topGO** in R, al fine di identificare processi biologici, funzioni molecolari e componenti cellulari sovra-rappresentati.
Questa integrazione consente di collegare i risultati filogenetici e di dinamica genica a potenziali differenze funzionali associate alle diverse linee evolutive di *Anopheles*, con particolare attenzione al complesso *gambiae*.
