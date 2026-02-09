## Genome Annotation and Functional Enrichment

### Annotazione funzionale

L’annotazione funzionale rappresenta un passaggio fondamentale per attribuire un significato biologico ai dati genomici ottenuti. Questo processo consente di collegare le sequenze proteiche a funzioni molecolari, processi biologici e componenti cellulari noti.
Esistono diversi strumenti per l’annotazione basata su similarità di sequenza, tra cui **BLAST**, **DIAMOND** e **HMMER**. BLAST utilizza un approccio di allineamento classico, mentre DIAMOND ne rappresenta un’alternativa estremamente veloce, ottimizzata per dataset di grandi dimensioni, con una minima perdita di sensibilità. HMMER, invece, si basa su modelli probabilistici (Hidden Markov Models) ed è particolarmente efficace nel riconoscere omologie remote interrogando database di profili come **Pfam** o quelli inclusi nel consorzio **InterPro**.

In questo lavoro è stato utilizzato **DIAMOND** come motore principale per l’annotazione basata su similarità, grazie alla sua efficienza computazionale e all’elevata sensibilità nel confronto con database proteici di grandi dimensioni.

---

### File di input per l’annotazione

Per garantire un’annotazione funzionale coerente e informativa, è stato generato un file di input contenente una singola sequenza proteica rappresentativa per ciascun ortogruppo.
In particolare, è stata selezionata **la proteina più lunga all’interno di ogni ortogruppo trim-mato**, evitando l’inclusione dei gap. Questa strategia consente di massimizzare l’informazione funzionale mantenendo una rappresentazione non ridondante delle famiglie geniche.

Il file è stato ottenuto utilizzando lo script:

```bash
bash longest_protein_OG.sh
```

Il risultato finale è il file:

```text
longest_sara.faa
```

che contiene le sequenze proteiche utilizzate per tutte le analisi successive.

---

### Database utilizzati

Le sequenze proteiche sono state confrontate con diversi database di riferimento:

* **Nr (Non-redundant)**: collezione non ridondante di proteine provenienti da GenPept, Swiss-Prot, PIR, PDB e RefSeq
* **Swiss-Prot**: database curato manualmente di proteine annotate (parte di UniProt)
* **Pfam**: raccolta di famiglie proteiche rappresentate da allineamenti multipli e profili HMM

---

### Annotazione con DIAMOND

DIAMOND è stato utilizzato per identificare omologhi proteici nei database di riferimento.
A causa di limitazioni computazionali del server, la costruzione del database e la ricerca DIAMOND sono state eseguite dal docente a partire dal file `longest_sara.faa`. I risultati dell’annotazione sono stati restituiti sotto forma dei seguenti file:

```text
sara_diamond.tsv
sara_diamond_names.tsv
```

Questi file contengono, rispettivamente, i risultati grezzi del confronto e una versione semplificata con le migliori annotazioni per ciascun ortogruppo.

---

### Annotazione Gene Ontology (GO)

Per associare le sequenze proteiche a termini funzionali standardizzati, è stata utilizzata l’annotazione Gene Ontology (GO).
L’annotazione GO è stata ottenuta tramite **InterProScan**, che integra diversi database (Pfam, PRINTS, SUPERFAMILY, ecc.) e consente di assegnare termini GO e pathway Reactome.

Anche in questo caso, l’esecuzione di InterProScan è stata effettuata dal docente a partire dal file `longest_sara.faa`. Il file di output principale è:

```text
longest_sara.tsv
```

che contiene, tra le altre informazioni, i termini GO associati a ciascuna proteina.

---

### Preparazione del background GO

Per l’analisi di arricchimento funzionale è stato costruito un **gene universe**, costituito da tutti i geni annotati con almeno un termine GO.
Il file è stato ottenuto rimuovendo metadati non necessari e riformattando le annotazioni GO nel formato richiesto da **topGO**:

```bash
awk -F'\t' '{
  gsub(/@.*/,"",$1);
  gsub(/\([^)]*\)/,"",$2);
  gsub(/\|/,",",$2);
  split($2,a,",");
  for(i in a) if(a[i]!="") seen[$1,a[i]]=1
}
END {
  for(k in seen){
    split(k,b,SUBSEP)
    groups[b[1]] = (groups[b[1]] ? groups[b[1]] "," b[2] : b[2])
  }
  for(g in groups) print g "\t" groups[g]
}' sara_id_go.tsv | grep -v "-" > go_back.tsv
```

Successivamente è stata creata una versione collassata del background:

```text
go_back_collapsed.tsv
```

---

### Gene Ontology Functional Enrichment

L’analisi di arricchimento funzionale è stata utilizzata per identificare termini GO sovra-rappresentati in un insieme di geni di interesse rispetto al background genomico.

I **geni di interesse** sono stati definiti come quelli appartenenti a **ortogruppi significativamente espansi o contratti**, identificati tramite l’analisi di evoluzione delle famiglie geniche con **CAFE5**.
L’elenco finale è contenuto nel file:

```text
interesting.txt
```

L’analisi di arricchimento è stata condotta in **R** utilizzando il pacchetto **topGO**, applicando il test di Fisher con algoritmo *elim*, che tiene conto della struttura gerarchica dei termini GO.

Sono state analizzate separatamente le tre ontologie:

* **BP** (Biological Process)
* **MF** (Molecular Function)
* **CC** (Cellular Component)

I risultati significativi (p-value < 0.05) sono riportati nei file:

```text
topGOe_name_interest_BP.txt
topGOe_name_interest_MF.txt
topGOe_name_interest_CC.txt
```

---

### Visualizzazione e interpretazione

I termini GO arricchiti ottenuti sono stati successivamente utilizzati per la visualizzazione tramite **REVIGO**, al fine di ridurre la ridondanza e facilitare l’interpretazione biologica dei risultati.

