# Genome Annotation

## Panoramica

L’annotazione genomica è la fase successiva all’assemblaggio e ha l’obiettivo di individuare e descrivere le diverse componenti funzionali presenti nella sequenza ottenuta.  
Tra queste rientrano:

- regioni codificanti per proteine  
- pseudogeni  
- elementi trasponibili  
- RNA non codificanti (ncRNA)  
- tRNA e altri elementi strutturali  

In pratica, a una sequenza di nucleotidi inizialmente priva di informazioni biologiche viene assegnata una funzione specifica e una posizione precisa all’interno del genoma.

Per questo processo vengono impiegati software di predizione genica come **SNAP** e **AUGUSTUS**, che utilizzano modelli addestrati per riconoscere le regioni codificanti.

L’annotazione viene solitamente eseguita in più cicli successivi (round), poiché ogni iterazione consente di affinare i modelli predittivi e migliorare progressivamente l’accuratezza dei risultati.

---

## Identificazione degli elementi ripetuti

Poiché il progetto è focalizzato principalmente sulle regioni codificanti, è necessario individuare e mascherare preventivamente le sequenze ripetitive del genoma.  

Questi elementi presentano strutture differenti rispetto ai geni e possono interferire con la predizione.

I principali strumenti utilizzati sono:

- **RepeatModeler**  
- **RepeatMasker**

---

## Esecuzione dell’annotazione con MAKER

Il software scelto per l’annotazione è **MAKER**, che richiede la creazione preliminare di file di configurazione:

```bash
maker -CTL
