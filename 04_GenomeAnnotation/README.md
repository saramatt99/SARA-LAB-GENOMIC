# Genome Annotation

L’annotazione genomica è la fase che segue l’assemblaggio di un genoma e rappresenta il momento in cui vengono individuate e descritte tutte le componenti funzionali presenti nella sequenza ottenuta. Durante questo processo vengono riconosciuti diversi elementi, tra cui **regioni codificanti per proteine**, **pseudogeni**, **elementi trasponibili**, **RNA non codificanti (ncRNA)**, **tRNA** e altre strutture genomiche rilevanti.

L’obiettivo principale dell’annotazione consiste nell’associare una **funzione biologica** a porzioni di DNA che inizialmente rappresentano soltanto sequenze anonime di nucleotidi, assegnando a ciascun elemento una posizione specifica all’interno del genoma.

Per svolgere questa attività vengono utilizzati software specializzati nella **predizione genica**, come **SNAP** e **AUGUSTUS**, basati su modelli statistici che possono essere addestrati per migliorare l’identificazione delle regioni codificanti.  

L’annotazione viene generalmente eseguita attraverso **più round successivi**, poiché ogni iterazione consente di affinare progressivamente la qualità delle predizioni.

---

## Identificazione e mascheramento degli elementi ripetuti

Poiché il progetto è focalizzato principalmente sulle regioni codificanti, è necessario individuare preliminarmente tutte le sequenze classificate come **repetitive elements**, che presentano strutture differenti rispetto ai geni e possono interferire con la predizione genica.

Per questa fase vengono utilizzati:

- **RepeatModeler** → costruzione di una libreria specifica di elementi ripetuti  
- **RepeatMasker** → mascheramento delle regioni ripetitive nel genoma  

---

## Annotazione genomica con MAKER

Il programma utilizzato per l’annotazione è **MAKER**, una pipeline che integra evidenze esterne e strumenti di predizione genica.

Per iniziare vengono generati i file di configurazione tramite:

```bash
maker -CTL
