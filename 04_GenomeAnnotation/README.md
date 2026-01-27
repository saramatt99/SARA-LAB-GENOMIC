Genome Annotation

L’annotazione genomica è la fase successiva all’assemblaggio di un genoma e rappresenta il momento in cui vengono individuate e descritte tutte le componenti funzionali presenti nella sequenza ottenuta. Durante questo processo vengono riconosciuti diversi elementi, tra cui regioni codificanti per proteine, pseudogeni, elementi trasponibili, RNA non codificanti (ncRNA), tRNA e altre strutture genomiche rilevanti.

L’obiettivo principale dell’annotazione consiste nell’associare una funzione biologica a porzioni di DNA che inizialmente rappresentano solo sequenze anonime di nucleotidi, assegnando a ciascun elemento una posizione specifica all’interno del genoma.

Per svolgere questa attività vengono utilizzati software specializzati nella predizione genica, come SNAP e AUGUSTUS, che si basano su modelli statistici e possono essere addestrati per migliorare l’identificazione delle regioni codificanti. L’annotazione viene generalmente eseguita attraverso più round successivi, poiché ogni iterazione consente di affinare progressivamente l’accuratezza dei modelli predittivi.

Identificazione e mascheramento degli elementi ripetuti

Poiché il progetto è focalizzato principalmente sulle regioni codificanti, è necessario individuare preliminarmente tutte le sequenze classificate come elementi ripetuti. Questi presentano strutture differenti rispetto ai geni e possono interferire con la predizione genica.

Per questa fase vengono utilizzati:

RepeatModeler, per la costruzione di una libreria specifica di sequenze ripetitive

RepeatMasker, per il mascheramento delle regioni ripetitive nel genoma

Annotazione genomica con MAKER

Il programma utilizzato per l’annotazione è MAKER, una pipeline che integra evidenze esterne e strumenti di predizione genica.

Per iniziare, vengono generati i file di configurazione tramite il comando:

maker -CTL


Questo comando produce una serie di file di controllo che permettono di impostare le opzioni di esecuzione. All’interno di questi file vengono specificati:

il genoma assemblato

il tipo di organismo (eucariotico)

il proteoma di riferimento della specie A. stephensis

la libreria di elementi ripetuti ottenuta con RepeatModeler

i parametri per i software di predizione genica

Tra le principali sezioni configurabili rientrano:

definizione del genoma di input

opzioni per eventuali ri-annotazioni

inserimento di evidenze basate su EST o mRNA

utilizzo di omologie proteiche

mascheramento delle sequenze ripetitive

parametri per SNAP e AUGUSTUS

impostazioni di performance (CPU, memoria)

Una volta configurati i file di controllo, MAKER viene eseguito per produrre i risultati di annotazione.

Riorganizzazione dei file di output

Al termine dell’esecuzione di MAKER, è possibile riorganizzare la directory di lavoro ed eliminare i file intermedi mediante il comando:

maker -base <OUTPUT_PREFIX>


Questo passaggio consente di ottenere una struttura di output più ordinata e facilmente consultabile.

Generazione del file GFF finale

Per unire tutte le informazioni prodotte da MAKER in un unico file utilizzabile nelle analisi successive, vengono eseguiti i seguenti comandi:

fasta_merge -d <DATASTORE_INDEX_FILE>
gff3_merge -d <DATASTORE_INDEX_FILE>


Questi comandi permettono di:

generare le sequenze nucleotidiche e amminoacidiche dei geni predetti

costruire il file .gff finale contenente tutte le annotazioni

Il file GFF è organizzato in colonne che riportano informazioni quali:

coordinate genomiche

tipo di elemento annotato (gene, CDS, exon, ecc.)

attributi funzionali associati

Predizione genica

La predizione delle regioni codificanti viene effettuata principalmente tramite due software:

SNAP (Semi-HMM-based Nucleic Acid Parser)

SNAP è un programma basato su modelli Hidden Markov (HMM), utilizzabile per genomi eucariotici e procariotici. Può essere addestrato utilizzando set di geni noti, migliorando la qualità delle predizioni.

AUGUSTUS

AUGUSTUS è un software di predizione genica di tipo ab initio che utilizza modelli statistici per individuare le regioni codificanti direttamente dalla sequenza genomica. Anche questo programma può essere allenato per aumentare l’accuratezza dei risultati.

I metodi ab initio non si basano sulla similarità con geni conosciuti, ma esclusivamente sulle caratteristiche intrinseche della sequenza di DNA.

Workflow generale e iterazioni successive

L’annotazione genomica è un processo lungo e complesso, ma rappresenta uno dei passaggi fondamentali nell’analisi di un genoma.

Il flusso di lavoro principale può essere sintetizzato nei seguenti punti:

Identificazione degli elementi ripetuti nel genoma

Mascheramento delle sequenze ripetitive

Integrazione di evidenze esterne (es. proteomi)

Generazione di un primo set di modelli genici affidabili

Addestramento dei software di predizione genica

Esecuzione di ulteriori round di annotazione per migliorare la predittività

Creazione di un insieme di risultati di consenso

Valutazione finale della qualità dell’annotazione
