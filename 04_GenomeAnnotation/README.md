## Genome Annotation

La fase di **genome annotation** ha come obiettivo l’identificazione e la localizzazione delle caratteristiche funzionali all’interno dell’assemblaggio genomico. Il risultato finale di questo processo è generalmente un file in formato **GFF** (General Feature Format) o **GTF**, che riporta le coordinate genomiche di geni, esoni, introni e altre feature annotate rispetto al genoma di riferimento.
Accanto all’annotazione strutturale, l’annotazione funzionale permette di associare a ciascuna sequenza una possibile funzione biologica, fornendo un primo livello di interpretazione del contenuto genomico.

---

## Mascheramento di elementi ripetitivi e trasposoni

Prima di procedere con la predizione genica, è fondamentale identificare ed escludere **regioni ripetitive e trasposoni**. Questi elementi possono infatti generare **falsi positivi**, poiché i software di gene prediction potrebbero erroneamente interpretare sequenze trasponibili come geni codificanti dell’ospite.

Per questo motivo, è stato adottato un workflow standard basato su **RepeatModeler** e **RepeatMasker**:

* **RepeatModeler** analizza il genoma in modalità *de novo* per costruire una libreria specifica di sequenze ripetitive.
* **RepeatMasker** utilizza questa libreria per individuare e mascherare le regioni ripetitive, impedendo che interferiscano con la predizione dei geni.

Poiché RepeatModeler è particolarmente oneroso dal punto di vista computazionale, la libreria di consenso è stata fornita direttamente dal docente ed è disponibile al seguente percorso:

```
/home/PERSONALE/mirko.martini3/01_2024/00_Data/02_annotation/Anoste_RepeatModeler_library.fa
```

Il mascheramento delle ripetizioni è stato eseguito internamente da MAKER durante l’annotazione e i risultati sono stati successivamente riassunti.

---

## Prima annotazione – MAKER (Round 1)

Per l’annotazione genomica è stato utilizzato **MAKER**, una pipeline modulare che integra evidenze sperimentali e predizioni *ab initio*.
A differenza di molti software che utilizzano lunghe stringhe di opzioni da linea di comando, MAKER si basa su **file di configurazione** dedicati, che rendono il workflow più chiaro e riproducibile.

Il comando seguente genera i tre file di controllo necessari:

```bash
maker -CTL
```

I file prodotti sono:

* **maker_exe.ctl**: contiene i percorsi agli eseguibili esterni (BLAST, RepeatMasker, ecc.). Generalmente viene riconosciuto automaticamente.
* **maker_bopts.ctl**: definisce le soglie di filtraggio per gli allineamenti. In questo lavoro sono stati mantenuti i parametri di default.
* **maker_opts.ctl**: rappresenta il file principale di configurazione ed è l’unico che richiede modifiche manuali.

### Parametri principali modificati in `maker_opts.ctl`

I parametri più rilevanti impostati per questa analisi sono:

* `genome=`: percorso al file FASTA dell’assemblaggio genomico.
* `protein=`: file FASTA contenente evidenze proteiche esterne.
* `model_org=`: lasciato vuoto per evitare bias nel mascheramento delle ripetizioni.
* `rmlib=`: percorso alla libreria personalizzata di RepeatModeler.
* `protein2genome=1`: abilita la predizione genica basata sull’allineamento delle proteine.
* `cpus=`: numero di core utilizzati.
* `pred_stats=1`: abilita il calcolo delle metriche di qualità (AED e QI).
* `min_protein=50`: lunghezza minima delle proteine predette.
* `alt_splice=1`: abilita la predizione di isoforme alternative.
* `split_hit=`: valore rimosso per evitare limiti predefiniti sugli allineamenti.

Una volta completata la configurazione, MAKER è stato avviato con:

```bash
maker -base <OUTPUT_PREFIX>
```

---

## Unione dei risultati e valutazione

Al termine della prima annotazione, i file GFF e FASTA prodotti sono stati uniti utilizzando:

```bash
fasta_merge -d <datastore_index_file>
gff3_merge -d <datastore_index_file>
```

La qualità dell’annotazione è stata quindi valutata tramite **AGAT**, che fornisce statistiche descrittive utili per verificare la coerenza biologica del gene set ottenuto.

```bash
agat_sp_statistics.pl --gff file.gff -o output_file
agat_sq_repeats_analyzer.pl -i input_file -o output_file
```

In questa fase è inoltre possibile valutare il proteoma risultante tramite BUSCO, per verificare la completezza dei geni codificanti.

---

## Seconda annotazione – addestramento SNAP e Augustus

### SNAP

**SNAP** è un gene predictor probabilistico basato su Hidden Markov Models (HMM), particolarmente efficace quando addestrato su dati specifici dell’organismo in analisi.
L’addestramento è stato effettuato a partire dai modelli genici ad alta confidenza ottenuti nel primo round di MAKER, utilizzando gli strumenti della suite SNAP (in particolare **Fathom**).

I principali passaggi sono stati:

```bash
maker2zff -c 0 -e 0 -l 80 -x 0.1 -d <datastore_index>
fathom *.ann *.dna -gene-stats
fathom *.ann *.dna -validate
fathom *.ann *.dna -categorize 1000
fathom uni.ann uni.dna -export 1000 -plus
forge export.ann export.dna
hmm-assembler.pl <NAME> <FORGE_DIR> > snap.hmm
```

---

### Augustus

**Augustus** è un gene predictor per genomi eucariotici basato su Generalized Hidden Markov Models.
Poiché l’addestramento manuale può essere computazionalmente impegnativo, è stato utilizzato **BUSCO** per addestrare automaticamente Augustus sfruttando geni altamente conservati.

```bash
busco -i Anoste_chr.fasta -c 30 -l $BUSCO/culicidae_odb12 \
--augustus --long -m genome --out Anoste_cu \
--augustus_parameters='--progress=true'
```

---

## Secondo round di MAKER

Per il secondo round di annotazione, il file `maker_opts.ctl` è stato duplicato e aggiornato con i seguenti cambiamenti:

* `protein_gff` e `rm_gff`: aggiunti i file GFF prodotti nel round precedente.
* `protein=`: svuotato.
* `snaphmm=`: percorso al modello SNAP addestrato.
* `augustus_species=`: impostato sul nome della specie addestrata.
* `protein2genome=0` e `est2genome=0`: disabilitati.
* `pred_stats=1`: mantenuto attivo.

Dopo l’esecuzione di MAKER, i risultati sono stati nuovamente uniti e analizzati.

---

## Valutazione finale del gene set

La qualità complessiva dell’annotazione è stata valutata confrontando le statistiche finali (lunghezza media dei geni, numero medio di esoni e introni, distribuzione delle feature) con valori noti in letteratura per specie affini. Questo confronto rappresenta uno dei criteri più affidabili per verificare la plausibilità biologica dell’annotazione.


