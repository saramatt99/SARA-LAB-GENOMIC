# Inferenza Filogenomica dagli Ortogruppi

Questa sezione descrive la pipeline utilizzata per l'inferenza filogenomica e l'analisi dell'evoluzione delle famiglie geniche.

## 1. Preparazione dei Dati e Ortologia

I genomi di partenza sono stati scaricati da NCBI e processati (pulizia header) come descritto nel `README` della directory `00_data`. Sebbene BUSCO sia uno standard per il controllo qualità, in questa specifica pipeline il passaggio è stato omesso.

L'identificazione dei gruppi di ortologia è stata effettuata con **OrthoFinder**. I risultati (spostati in questa directory) forniscono gruppi che includono sia ortologhi che paraloghi. Il file `Statistics_Overall.csv` generato da OrthoFinder contiene il riepilogo delle statistiche per specie (numero di geni, ortogruppi, ecc.).

## 2. Filtraggio dei Paraloghi (DISCO)

Per districare le relazioni evolutive complesse e risolvere i gruppi contenenti paraloghi, è stato utilizzato **DISCO**. Questo tool scompone ("spacchetta") gli ortogruppi complessi in gruppi di ortologhi a singola copia (anche se non necessariamente presenti in tutte le specie), trasformando un albero genico complesso in sotto-alberi risolti.

### Esecuzione di DISCO

A differenza di OrthoFinder, DISCO richiede un file separato per ogni albero. Partendo dal file `Resolved_Gene_Trees.txt` (dove gli alberi sono nella seconda colonna), abbiamo iterato su ogni riga per isolare gli alberi e lanciato DISCO.

```bash
# Ciclo per estrarre gli alberi e lanciare DISCO
while IFS=' ' read -r OG tree; do 
    python3 ../99_scripts/disco.py \
    -i <(echo "$tree") \
    -o ../01_Disco/${OG/:/}.nwk \
    -d "|" \
    -m 4 \
    --remove_in_paralogs \
    --keep-labels \
    --verbose >> ../01_Disco/disco.log
done < <(sed -E 's/[A-Z][a-z]{5}_//g; s/\)n[0-9]*+/\)/g' Resolved_Gene_Trees.txt)

```

### Pulizia e Organizzazione Output

Alcuni file generati potrebbero risultare vuoti. Per pulizia, vengono identificati e rimossi:

```bash
# Identificazione file vuoti (backup lista)
find . -size 0 -print > empty_disco.txt

# Rimozione file vuoti
find . -size 0 -delete

```

Successivamente, l'output di DISCO è stato processato per separare le sequenze fasta corrispondenti ai nuovi gruppi risolti:

```bash
# Split delle sequenze basato sui risultati DISCO
bash ../99_scripts/split_disco_output.sh <PATH_TO_ORTHOGROUP_SEQUENCES>

```

---

## 3. Costruzione dello Species Tree

Per studiare l'espansione e la contrazione delle famiglie geniche (tramite CAFE) è necessario un albero di specie calibrato (Time Tree). Per costruirlo, abbiamo selezionato un subset di ortologhi a singola copia (Single Copy Orthologues - SCO) completi.

> **Nota:** In assenza di SCO completi da OrthoFinder, si sarebbero potuti utilizzare i risultati di DISCO o procedere con un numero ridotto di specie.

### Selezione e Allineamento

Sono stati estratti casualmente 200 ortogruppi a singola copia per la costruzione dell'albero:

```bash
# Estrazione randomica di 200 ortogruppi
ls 00_OrthoFinder/Results/Single_Copy_Orthologue_Sequences | shuf -n 200 > species_tree_OG.txt

# Allineamento con MAFFT
for OG in $(cat species_tree_OG.txt); do 
    mafft --auto --anysymbol "$OG" > ../03_aligned/${OG/.fa/_aligned.faa}
done

```

### Trimming e Concatenazione

Le sequenze allineate sono state "trimmate" con **BMGE** per rimuovere regioni scarsamente conservate, e gli header sono stati puliti per compatibilità con i software successivi.

```bash
# Trimming con BMGE
for OG in *; do 
    bmge -i "$OG" -t AA -m BLOSUM62 -e 0.5 -g 0.4 -of ../04_trimmed/${OG/_aligned.faa/_trimmed.faa}
done

# Pulizia header (rimozione pipe e caratteri successivi)
sed -i.old -E 's/\|.+$//' *

# Concatenazione delle sequenze (AMAS)
python3 ../99_scripts/AMAS.py concat -y nexus -i *.faa -f fasta -d aa -t conc_species_tree

```

### Inferenza dell'Albero (IQ-TREE)

L'albero filogenetico è stato calcolato utilizzando **IQ-TREE** con 100 bootstrap:

```bash
iqtree -m TESTNEW -b 100 -s conc_species_tree --prefix species_tree -nt 9

```

---

## 4. Preparazione per l'Analisi delle Famiglie Geniche (CAFE)

Per l'analisi di espansione e contrazione, non ci limitiamo ai single copy, ma utilizziamo l'intero dataset derivato da DISCO (filtrato dai paraloghi).

### Allineamento e Trimming del Dataset Completo

La procedura è analoga a quella dello species tree, ma applicata a tutti i file fasta generati nella fase `02_disco`.

```bash
# Allineamento di tutte le famiglie geniche
for file in *faa; do 
    mafft --auto --anysymbol "$file" > ../03_aligned/prova/${file/.faa/_aligned.faa}
done

# Trimming
for file in *; do 
    bmge -i "$file" -t AA -m BLOSUM62 -e 0.5 -g 0.4 -of ../../04_trimmed/prova/${file/_aligned.faa/_trimmed.faa}
done

```
