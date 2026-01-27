# Analisi dell’evoluzione delle famiglie geniche

Per studiare l’espansione e la contrazione delle famiglie geniche utilizziamo **CAFE**, basato su un modello di **birth-and-death**. Per eseguire l’analisi servono due elementi principali:

1. Un **albero temporale** (time tree) in formato `.nwk`
2. Una **tabella con il numero di geni per famiglia e specie**, che possiamo prendere da Orthofinder (`Orthogroups.GeneCount.tsv`) ma va prima modificata.

---

## Preparazione dei dati

L’albero originale è in formato `.nex`. Lo abbiamo caricato su **iTOL**, esportato come `.nwk` e salvato in `timetree.nwk`.

La tabella di Orthofinder va modificata per rimuovere l’ultima colonna e aggiungere una prima colonna di valori casuali `NONE`:

```bash
sed $'s/^/NONE\t/g' Orthogroups.GeneCount.tsv | rev | cut -f 2- | rev > ../../../../07_GeneFamilies_Evolution/GeneCount_CAFE.tsv
```

---

## Avvio di CAFE

Per stimare il **modello di errore generale**, lanciamo:

```bash
cafe5 -i GeneCount_CAFE.tsv -t timetree.nwk -o Error_model -e
```

Di default, CAFE assume lo stesso turnover per tutte le specie (una sola λ). Possiamo però esplorare modelli più complessi con λ diverse per specie o clade, oppure usare il parametro **gamma** per variare il turnover tra famiglie geniche.

---

## Replicati tecnici

Per testare la variabilità dello strumento facciamo **10 replicati** per ciascun livello di complessità `k = 1..5`. Creiamo le cartelle e lanciamo l’analisi 1-lambda:

```bash
for k in {1..5}; do for n in {1..10}; do mkdir -p 00_1L/${k}K/${n}N; cafe5 -i GeneCount_CAFE.tsv -t timetree.nwk -o 00_1L/${k}K/${n}N -eError_model/Base_error_model.txt -k ${k}; done; done
```

---

## Analisi con 2-lambda

Per fare analisi con λ diverse tra clade, prepariamo un albero `.nwk` con le etichette dei lambda:

```
((Anogam:2,Anoste:2):1,(Culqui:2,(Sabcya:1,(Aedalb:1,Aedaeg:1):1):1):1);
```

E lanciamo i replicati:

```bash
for k in {1..5}; do for n in {1..10}; do mkdir -p 00_2L/${k}K/${n}N; cafe5 -i GeneCount_CAFE.tsv -t timetree.nwk -o 00_2L/${k}K/${n}N -y timetree2L.nwk -eError_model/Base_error_model.txt -k ${k}; done; done
```

---

## Scelta del modello

* L’aumento dei parametri (`k` o multi-lambda) aumenta la **likelihood**, ma bisogna verificare se il modello più complesso è giustificato.
* È consigliato fare **replicati multipli** e confrontare lnL tra run.

---

## Risultati di interesse

* **Base_results.txt** → lnL finale, lambda (turnover per milione di anni), epsilon (percentuale famiglie stabili), lambda massima, numero di iterazioni
* **Base_asr.tre** → stato ancestrale dei nodi; numeri in `<>` = nodi, numeri senza `<>` = membri della famiglia per specie
* **Base_change.tab** → differenze tra membri di un nodo e il nodo precedente per ogni famiglia
* **Base_clade_results.txt** → numero di famiglie significativamente aumentate o diminuite per clade
* **Base_count.tab** → conteggio dei membri per famiglia e specie
* **Base_family_results.txt** → colonne: famiglia, p-value per contrazione/espansione, significatività (y/n)

> Con gamma al posto di Base (quando viene usato il parametro gamma) si producono gli stessi file con prefisso `Gamma_`, con in più il parametro **alpha** che regola la distribuzione dei tassi di evoluzione tra famiglie.

---


