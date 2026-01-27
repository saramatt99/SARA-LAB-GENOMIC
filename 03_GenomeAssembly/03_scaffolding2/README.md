# Scaffolding del Genoma

## 1. Creazione degli scaffold

In questa fase dell’assemblaggio genomico, si costruiscono gli **scaffold**.  
Uno scaffold rappresenta una versione virtuale dei contig, ordinati secondo criteri specifici per riflettere la struttura cromosomica della specie studiata.

---

## 2. Miglioramento e correzione dell’assemblaggio

Prima di creare gli scaffold, è possibile migliorare la qualità dell’assemblaggio confrontandolo con un **genoma di riferimento**.  
Questo passaggio permette di correggere eventuali errori nei contig.

> Per assemblaggi de novo, è consigliabile utilizzare un genoma di riferimento di una specie filogeneticamente vicina a quella in studio.

```bash
ragtag.py correct -t 20 <REFERENCE_GENOME> <DRAFT_GENOME>
