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

## Download dei genomi

I genomi utilizzati nell’analisi sono stati scaricati da **NCBI**, selezionando esclusivamente assembly che rispettassero alcuni criteri minimi di qualità, in modo da garantire coerenza e affidabilità nelle analisi comparative.
I filtri applicati sono stati:

* genoma di riferimento
* genoma annotato
* livello minimo di assemblaggio: *scaffold*
* data di rilascio successiva al 2010

Il download è stato automatizzato tramite uno script Bash che utilizza il tool `ncbi-datasets` per scaricare sia il genoma che il file di annotazione GFF.

### Script di download

```bash
#!/bin/bash

# Script per scaricare genomi e file GFF da NCBI datasets
# Crea cartelle pronte all'uso per le analisi successive

AN2name=$1

mkdir 00_genome
mkdir 01_gff

while IFS=$'\t' read -r AN sname ID; do
    echo $AN
    datasets download genome accession "$AN" --filename "$ID".zip --include genome,gff3
    unzip "$ID".zip -d "$ID"
    mv "$ID"/ncbi_dataset/data/"$AN"/*.fna 00_genome/"$ID".fna
    mv "$ID"/ncbi_dataset/data/"$AN"/*.gff 01_gff/"$ID".gff
    rm -r "$ID"/
done < "$AN2name"
```

Per avviare il download:

```bash
nano dataset.tsv
conda activate sequence
bash download_dataset.sh dataset.tsv
```

---

## Selezione dell’isoforma più lunga

In seguito all’annotazione genomica, uno stesso gene può essere rappresentato da più isoforme. Tuttavia, mantenere tutte le isoforme può introdurre ridondanza e bias nelle analisi comparative.
Per questo motivo è stata selezionata **una sola isoforma per gene**, mantenendo esclusivamente la più lunga.

A tale scopo è stato utilizzato **AGAT (Another GFF Analysis Toolkit)**.
Il file GFF è stato filtrato con `agat_sp_keep_longest_isoform.pl` e, successivamente, le sequenze codificanti sono state estratte e tradotte in proteine.

```bash
cd 01_gff
for gff in *.gff; do
    agat_sp_keep_longest_isoform.pl -gff "$gff" -o ${gff/.gff/_longest.gff}
done
```

```bash
cd 00_data
mkdir 02_raw_proteoms
for gff in 01_gff/*_longest.gff; do
    agat_sp_extract_sequences.pl \
    -g "$gff" \
    -f ../00_genome/${gff/_longest.gff/.fna} \
    -t cds -p --cfs \
    --output 02_raw_proteoms/${gff/_longest.gff/faa}
done
```

---

## Rimozione dei pseudogeni

Per assicurare che il dataset finale includesse solo proteine funzionali, sono state eliminate tutte le sequenze contenenti **codoni di stop interni**, indicative di pseudogeni.

È stato utilizzato uno script Bash che:

1. converte i file FASTA in formato a singola riga,
2. individua le sequenze contenenti `*`,
3. rimuove automaticamente le sequenze corrispondenti.

### Script per identificare ed eliminare pseudogeni

```bash
#!/bin/bash

mkdir raw_proteomes
mv *.faa raw_proteomes/

cd raw_proteomes
for proteome in *.faa; do
    awk '{if(NR==1) {print $0} else {if($0 ~ /^>/) {print "\n"$0} else {printf $0}}}' "$proteome" > ../${proteome/.faa}".faa"
done
cd ..

mkdir 00_pseudogene_name

for proteome in *.faa; do
    species=$(basename -s .faa "$proteome")
    grep -B1 '*' "$proteome" | grep ">" >> 00_pseudogene_name/"$species"_pseudogenes_name.txt
done

for pseudo_file in 00_pseudogene_name/*_pseudogenes_name.txt; do
    species=$(basename -s _pseudogenes_name.txt "$pseudo_file")
    while IFS=$'\t' read -r header; do
        sed -E -i "/${header}/{N;d;}" "$species".faa
    done < "$pseudo_file"
done

mv 00_pseudogene_name ../00_genome
```

Esecuzione dello script:

```bash
cd 02_raw_proteoms
bash ../../99_scripts/pseudogenes_find_eliminate.sh
```

---

## Modifica degli header FASTA

Per facilitare le analisi successive (OrthoFinder, DISCO, CAFE), gli header delle proteine sono stati semplificati mantenendo solo l’identificativo di specie e di gene.
Questa operazione è fondamentale per evitare ambiguità nei passaggi di inferenza ortologica.

```bash
for prot in *.faa; do
    ID=$(basename -s .faa "$prot")
    sed -i.old -E "s/>(rna-XM_[0-9]+\.[0-9]) (gene=gene-.[^ ]+) name=(.[^ ]+) .+$/>${ID}\|\3/" "$prot"
done
```
