# Genome Annotation

L'annotazione è il processo di caratterizzazione funzionale di un genoma assemblato. L'obiettivo è mappare elementi biologici (CDS, pseudogeni, trasposoni, ncRNA, tRNA) su sequenze nucleotidiche anonime, assegnando loro una posizione e una funzione specifica.
Il workflow utilizza software come **SNAP** e **AUGUSTUS** e viene eseguito in più **round iterativi** per affinare progressivamente i modelli statistici e migliorare l'accuratezza della predizione.

### Annotazione di elementi ripetuti
Per evitare interferenze con la predizione dei geni codificanti, è necessario identificare e mascherare preventivamente i **repetitive elements**.
> Tool utilizzati: `RepeatModeler` (identificazione) e `RepeatMasker` (mascheramento).

-----

### Pipeline MAKER
**MAKER** è il framework principale utilizzato per l'annotazione. La configurazione avviene tramite la generazione e la modifica dei file di controllo:

```bash
maker -CTL
```

In questo progetto, l'annotazione è stata guidata dal proteoma di riferimento di A. stephensis. I parametri principali nel file maker_opts.ctl includono:
#-----Input Genomico
    genome= <GENOME> 
    organism_type=eukaryotic 

    #-----Evidenze di Omologia (Proteoma A. stephensis)
    protein= <PROTEOME> 
    protein2genome=1 

    #-----Repeat Masking
    rmlib= <RepeatModeler_library> 
    softmask=1 

    #-----Ab-initio Gene Prediction
    snaphmm= <SNAP_HMM>
    augustus_species= <SPECIES>
   
  Per ottimizzare la directory di output e rimuovere i file temporanei:
  ```bash
maker -base <OUTPUT_PREFIX>
  ```
Consolidamento dei Risultati
Al termine dei round di computazione, i dati vengono unificati per generare il file .gff finale e i database FASTA delle sequenze predette.
```bash
# Merge dei risultati distribuiti
fasta_merge -d <DATASTORE_INDEX_FILE>
gff3_merge -d <DATASTORE_INDEX_FILE>
```
Il file .gff risultante contiene le coordinate genomiche e le annotazioni funzionali strutturate in colonne standardizzate.

Software di Predizione Genica
Il processo si affida a due motori di predizione complementari:

SNAP: Predittore basato su modelli HMM, addestrato su set di geni noti per aumentare la specificità.

AUGUSTUS: Software di predizione Ab initio che utilizza modelli statistici intrinseci del DNA (siti di splicing, contenuto GC) indipendentemente dalle evidenze esterne di omologia.
Workflow Sintetico
Masking: Mascheramento delle ripetizioni genomiche.

Evidence Alignment: Allineamento del proteoma di riferimento.

Training: Estrazione di modelli ad alta confidenza per addestrare SNAP e AUGUSTUS.

Iterazione: Esecuzione di round successivi per ottimizzare la sensibilità predittiva.

Final Consensus: Generazione dell'annotazione finale e valutazione della qualità.
