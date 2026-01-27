# Genome Annotation Pipeline

L'annotazione rappresenta la fase cruciale successiva all'assemblaggio: è il processo in cui la sequenza genomica grezza viene analizzata per identificare e descrivere le entità biologiche presenti, come geni codificanti, pseudogeni, trasposoni, ncRNA e tRNA. In termini pratici, l'obiettivo è assegnare una funzione biologica e una coordinata precisa a quella che, inizialmente, è solo una stringa anonima di nucleotidi.

Il workflow si basa su cicli iterativi (round) in cui algoritmi di predizione vengono progressivamente addestrati per affinare la capacità di identificazione dei geni.

---

## 1. Identificazione degli Elementi Ripetuti
Poiché il focus del progetto è sulle regioni codificanti, è fondamentale individuare e "mascherare" gli elementi ripetitivi, che presentano caratteristiche strutturali diverse dai geni.
* **Tool utilizzati:** `RepeatModeler` per la creazione di librerie di ripetizioni e `RepeatMasker` per il mascheramento.

## 2. Configurazione ed Esecuzione di MAKER
**MAKER** è l'orchestratore principale utilizzato per l'annotazione. Il processo inizia con la generazione dei file di controllo necessari per definire i parametri di analisi:

```bash
maker -CTL
Configurazione del file maker_opts.ctl
Per questo progetto, il software è stato configurato integrando il proteoma di A. stephensis come evidenza di omologia. Di seguito i parametri principali impostati:

Plaintext
#-----Genome
genome= <GENOME_FILE> 
organism_type=eukaryotic

#-----Protein Homology Evidence
protein= <PROTEOME_FILE> 

#-----Repeat Masking
rmlib= <REPEAT_LIBRARY>
softmask=1 

#-----Gene Prediction
snaphmm= <SNAP_HMM>
augustus_species= <SPECIES_MODEL>
protein2genome=1
Per ottimizzare la directory di lavoro ed eliminare i file intermedi dopo il run, si utilizza:

Bash
maker -base <OUTPUT_PREFIX>
3. Consolidamento dei Risultati
Dopo l'elaborazione, è necessario estrarre e unificare le informazioni frammentate nei vari file di output in un unico database consultabile (formato .gff) e nelle relative sequenze fasta (proteine e trascritti).

Bash
# Merge dei file di output
fasta_merge -d <DATASTORE_INDEX_FILE>
gff3_merge -d <DATASTORE_INDEX_FILE>
Il file .gff finale contiene la mappatura dettagliata di ogni sequenza riconosciuta come codificante, suddivisa in colonne standardizzate.

4. Algoritmi di Predizione Genica
La pipeline sfrutta due software principali che operano secondo logiche complementari:

SNAP (Semi-HMM-based Nucleic Acid Parser): Un programma flessibile per genomi eucariotici e procariotici. Viene addestrato su modelli genetici noti per elevare la qualità dell'annotazione.

AUGUSTUS: Un predittore che può operare in modalità ab initio. A differenza dei metodi basati sulla somiglianza, i modelli ab initio si affidano esclusivamente a proprietà statistiche intrinseche della sequenza di DNA.

Sintesi del Processo Iterativo
L'annotazione non è un passaggio unico, ma un ciclo di miglioramento continuo:

Mascheramento degli elementi ripetitivi.

Allineamento di evidenze esterne (es. proteomi).

Estrazione di modelli genici iniziali ad alta affidabilità.

Training di SNAP e AUGUSTUS sui modelli estratti.

Esecuzione di round successivi per raffinare la predittività.

Validazione e generazione del consenso finale.
