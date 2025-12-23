# Analisi genomica comparativa nel genere *Anopheles*

## Dataset

Per rispondere alla domanda biologica di interesse, è stato costruito un dataset comprendente **sei specie appartenenti al genere *Anopheles*** (famiglia Culicidae). Le specie sono state selezionate in modo da rappresentare **linee evolutive differenti**, includendo membri del **complesso *Anopheles gambiae***, specie filogeneticamente più distanti e taxa provenienti da differenti aree biogeografiche (Africa, Asia e Nuovo Mondo).

## Domanda biologica

L’analisi mira a investigare le **relazioni filogenetiche tra specie del genere *Anopheles*** e a valutare se il **complesso *Anopheles gambiae*** presenti un pattern evolutivo distinto rispetto alle altre linee del genere, in termini di divergenza filogenetica ed evoluzione delle famiglie geniche.

---

##Struttura del dataset

La tabella seguente riassume le informazioni associate a ciascuna specie inclusa nel dataset.

| Accession number | Nome specie              | ID     |
| ---------------- | ------------------------ | ------ |
| GCF_943734735.2  | *Anopheles gambiae*      | Anogam |
| GCF_943734845.2  | *Anopheles arabiensis*   | Anoara |
| GCF_013758885.1  | *Anopheles albimanus*    | Anoalb |
| GCF_013141755.1  | *Anopheles stephensi*    | Anoste |
| GCF_943734845.2  | *Anopheles funestus*     | Anofun |
| GCF_943734845.2  | *Anopheles culicifacies* | Anocul |

L’**Accession number** rappresenta l’identificativo univoco del genoma nella banca dati NCBI ed è stato utilizzato come riferimento per il download delle sequenze genomiche e dei file di annotazione associati.

Il dataset è stato salvato in formato **TSV**, così da garantire una struttura tabulare adatta alle successive fasi di analisi.

---

## Download dei dati

Per ciascuna specie sono stati scaricati il genoma completo e il relativo file di annotazione (`.gff`), utilizzando gli Accession number riportati nella tabella. L’intero processo è stato automatizzato tramite uno script dedicato.

```bash
bash ../99_scripts/download_dataset.sh dataset.tsv
```

---

## Pre-processing dei dati

Una volta acquisiti i genomi e le annotazioni, i dati sono stati sottoposti a una fase di pre-processing, finalizzata alla **standardizzazione dei file** e alla **rimozione di informazioni ridondanti**, in preparazione alle analisi comparative successive.

### Selezione dell’isoforma più lunga

Poiché a ciascun gene possono essere associate più isoforme, è stata mantenuta esclusivamente **l’isoforma più lunga**, al fine di evitare ridondanze nei dataset proteici.

```bash
for gff in *.gff; do
  agat_sp_keep_longest_isoform.pl --gff "$gff" -o ${gff/.gff/_longest.gff}
done
```

### Creazione dei proteomi

A partire dai file `.gff` filtrati e dai genomi associati, sono stati generati i **proteomi** per ciascuna specie, estraendo le sequenze codificanti e traducendole in amminoacidi.

```bash
for gff in *_longest.gff; do
  agat_sp_extract_sequences.pl \
  -g "$gff" \
  -f ../00_genome/${gff/_longest.gff/.fna} \
  -t cds -p --cfs \
  --output ../02_Proteome/${gff/_longest.gff/.faa}
done
```

### Rimozione dei pseudogeni

I genomi scaricati possono contenere sequenze annotate come **pseudogeni**, che non rappresentano CDS funzionali. Tali sequenze sono state identificate ed eliminate tramite uno script dedicato.

```bash
bash ../../99_scripts/pseudogene_find_eliminate.sh
```

### Standardizzazione delle intestazioni

Per rendere i file di proteoma più leggibili e compatibili con le analisi successive, è stata effettuata una **modifica delle intestazioni (header)**, uniformando l’identificativo dei geni per ciascuna specie.

```bash
for prot in *.faa; do
  ID=$(basename -s .faa "$prot")
  sed -i.old -E "s/>(rna-XM_[0-9]+\.[0-9]) (gene=gene-.[^ ]+) name=(.[^ ]+) .+$/>${ID}\|\3/" "$prot"
done
```

---

## Analisi filogenetica

L’albero filogenetico è stato ricostruito utilizzando un approccio **Maximum Likelihood (ML)**.
Tra le specie incluse nel dataset, ***Anopheles culicifacies*** è stata utilizzata come **outgroup**, in quanto filogeneticamente più distante rispetto alle altre specie analizzate, consentendo una corretta radicazione dell’albero.

---

## Analisi dell’evoluzione delle famiglie geniche (CAFE)

L’evoluzione delle famiglie geniche è stata analizzata utilizzando **CAFE**, implementando un modello con **due λ**, distinguendo il **complesso *Anopheles gambiae*** (*A. gambiae* e *A. arabiensis*) dalle restanti specie del genere *Anopheles*.
Questa scelta riflette l’ipotesi di differenti tassi di guadagno e perdita genica associati a linee evolutive con storie adattative distinte.

---

Informazioni sull’ambiente di lavoro

La cartella Environment Information fornisce una panoramica chiara e organizzata degli ambienti di lavoro (conda environment) utilizzati nel progetto.

Gestione degli environment

Il passaggio da un environment a un altro può essere effettuato tramite:

conda activate <nome_environment>
Environment utilizzati

base

R 4.5.1

tree

cafe 5.1.0

disco 1.4.1

ete3 3.1.3

gotree 0.4.5

HYPHY 2.5.71

iqtree 3.0.1

mafft 7.526

paml 4.10.7

orthofinder 2.5.5

raxml-ng 1.2.2

treeswift 1.1.45

sequence

agat 1.4.1

blobtools 1.1.1

bmge 1.12

BUSCO 6.0.0

edirect 24.0

mafft 7.526

ncbi-datasets 18.3.1

SRA-tools 3.2.1

assembly

python 3.13.5

assembly-stats 1.0.1

augustus 3.1

blast 2.16.0

diamond 2.1.10

fastqc 0.12.1

hypo 1.0.3

maker 3.01.04

minimap2 2.28

mosdepth 0.3.10

multiqc 1.31

r-base 4.3.3

samtools 1.21

spades 4.2.0

trimmomatic 0.40

kat

kat 2.4.2

Note sugli environment

La suddivisione dei software in environment distinti garantisce la gestione efficiente delle dipendenze e riduce potenziali conflitti tra versioni incompatibili.

Software e requisiti

AGAT

CAFE

Software per inferenza filogenetica ML

Bash

Note finali

Le lunghezze dei rami dell’albero ML rappresentano il numero di sostituzioni per sito e non devono essere interpretate come tempi di divergenza.
