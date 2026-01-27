# Genome Annotation

## Introduzione

L’annotazione genomica rappresenta una delle fasi fondamentali successive all’assemblaggio di un genoma.  
In questo passaggio vengono individuate, classificate e descritte le diverse componenti funzionali presenti nella sequenza genomica ottenuta.

Tra gli elementi annotati rientrano:

- sequenze codificanti per proteine  
- pseudogeni  
- elementi trasponibili  
- RNA non codificanti (ncRNA)  
- tRNA e altri RNA strutturali  

Il processo consiste nell’attribuire una funzione biologica a regioni che inizialmente rappresentano soltanto sequenze anonime di nucleotidi, associando ciascuna funzione a una specifica posizione genomica.

Per svolgere questa attività vengono utilizzati software specializzati nella predizione genica, come **SNAP** e **AUGUSTUS**, che sono addestrati per riconoscere le regioni codificanti sulla base di modelli genetici noti.

Generalmente l’annotazione non viene eseguita in un’unica fase, ma attraverso diversi cicli successivi (round).  
Ogni round consente di migliorare progressivamente le prestazioni degli algoritmi, aumentando l’accuratezza nella identificazione dei geni.

---

## Annotazione delle sequenze ripetute

Nel contesto di questo progetto l’attenzione è rivolta principalmente alle regioni codificanti del genoma.  
Per questo motivo risulta necessario individuare preliminarmente tutti gli elementi classificati come **repetitive elements**, che presentano caratteristiche strutturali differenti rispetto alle regioni geniche.

Una volta identificati, questi elementi vengono mascherati in modo da non interferire con la predizione dei geni.

I principali strumenti utilizzati per questa operazione sono:

- **RepeatModeler**, per la costruzione di una libreria specifica di elementi ripetuti  
- **RepeatMasker**, per il mascheramento delle sequenze ripetitive nel genoma  

---

## Generazione dell’annotazione con MAKER

Il programma utilizzato per l’annotazione genomica è **MAKER**, una pipeline che integra evidenze esterne e software di predizione genica.

Per iniziare è necessario generare i file di configurazione tramite il comando:

```bash
maker -CTL
Questo comando crea una serie di file di controllo che permettono di impostare le opzioni di esecuzione.

Una volta configurati i parametri, MAKER viene eseguito fornendo:

il genoma assemblato

informazioni provenienti dal proteoma della specie A. stephensis

la libreria degli elementi ripetuti

i modelli per la predizione genica

Di seguito sono riportate le principali sezioni dei file di configurazione:

#----- Genome
genome=<GENOME>
organism_type=eukaryotic

#----- Re-annotation options
maker_gff=
est_pass=0
altest_pass=0
protein_pass=0
rm_pass=0
model_pass=0
pred_pass=0
other_pass=0

#----- EST evidence
est=
altest=
est_gff=
altest_gff=

#----- Protein homology
protein=<PROTEOME>
protein_gff=

#----- Repeat masking
model_org=<EMPTY>
rmlib=<RepeatModeler_library>
repeat_protein=
rm_gff=
prok_rm=0
softmask=1

#----- Gene prediction
snaphmm=
gmhmm=
augustus_species=
fgenesh_par_file=
pred_gff=
model_gff=
est2genome=0
protein2genome=<0 OR 1>
trna=0
snoscan_rrna=
unmask=0

#----- Additional features
other_gff=

#----- External applications
alt_peptide=C
cpus=<CPUS>

#----- MAKER behavior
max_dna_len=100000
min_contig=1

pred_flank=200
pred_stats=<0 OR 1>
AED_threshold=1
min_protein=50
alt_splice=0
always_complete=0
map_forward=0
keep_preds=0

split_hit=10000
single_exon=0
single_length=250
correct_est_fusion=0

tries=2
clean_try=0
clean_up=0
TMP=
Riorganizzazione della directory di lavoro
Una volta completata l’esecuzione di MAKER, è possibile riorganizzare i file di output ed eliminare quelli intermedi mediante il comando:

maker -base <OUTPUT_PREFIX>
Questo passaggio consente di ottenere una struttura di directory più ordinata e facilmente consultabile.

Ricostruzione del file GFF finale
Per unire tutte le informazioni prodotte da MAKER in un unico file utilizzabile nelle analisi successive, vengono eseguiti i seguenti comandi:

fasta_merge -d <DATASTORE_INDEX_FILE>
gff3_merge -d <DATASTORE_INDEX_FILE>
Questi comandi permettono di:

generare le sequenze nucleotidiche e amminoacidiche dei geni predetti

costruire il file .gff finale contenente tutte le annotazioni

Il file GFF è organizzato in colonne che riportano informazioni come:

coordinate genomiche

tipo di feature (gene, CDS, exon, ecc.)

attributi funzionali

Predizione dei geni
Il processo di predizione genica viene effettuato tramite due software principali altamente specializzati:

SNAP (Semi-HMM-based Nucleic Acid Parser)
SNAP è un programma generico per la predizione dei geni, applicabile sia a genomi eucariotici che procariotici.
Si basa su modelli Hidden Markov (HMM) e può essere addestrato utilizzando set di geni noti, migliorando così la qualità delle predizioni.

AUGUSTUS
AUGUSTUS è un software di predizione genica di tipo ab initio.
Anche questo programma utilizza modelli statistici e può essere allenato per aumentare l’accuratezza dei risultati.

I metodi ab initio non si basano su omologie con geni conosciuti, ma esclusivamente sulle caratteristiche intrinseche della sequenza di DNA e sui modelli probabilistici dei geni.

Riassunto del workflow e iterazioni successive
L’annotazione genomica è un processo complesso e dispendioso in termini di tempo, ma rappresenta uno dei passaggi più importanti nell’analisi di un genoma.

Il flusso di lavoro principale può essere riassunto nei seguenti punti:

Identificazione degli elementi ripetuti presenti nel genoma

Mascheramento delle sequenze ripetitive

Utilizzo di evidenze esterne (come proteomi) per migliorare la predizione

Estrazione di un primo set di modelli genici affidabili

Addestramento dei software di predizione (SNAP e AUGUSTUS)

Esecuzione di ulteriori round di annotazione ripetendo i passaggi di training

Creazione di un insieme di risultati di consenso

Valutazione finale della qualità dell’annotazione
