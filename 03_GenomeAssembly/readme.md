## Genome Assembly

All'interno di questa cartella vengono gestite tutte le fasi necessarie per l'assemblaggio del genoma.

Il workflow principale è suddiviso in quattro passaggi fondamentali:

1. **Assemblaggio a livello contig**: creazione dei contig e del consensus, utilizzando le reads disponibili.

2. **Pulizia**: rimozione di artefatti ed errori accumulatisi durante l'assemblaggio, sfruttando sia short che long reads.

3. **Decontaminazione**: identificazione e rimozione di porzioni del genoma assemblato provenienti da organismi diversi da *A. stephensi*.

4. **Assemblaggio a livello scaffolding**: unione dei contig ottenuti in scaffold, per ottenere un genoma più completo e continuo.
