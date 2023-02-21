# Mouse-GEM data

This directory contains datasets that were used for the generation and curation of Mouse-GEM. The files and their contents are briefly summarized below.


- `human2MouseOrthologs.tsv`: The human-mouse ortholog pairs downloaded from the [Alliance of Genomes Resources](https://www.alliancegenome.org) using following code:
```bash
curl -X \
GET "https://www.alliancegenome.org/api/homologs/9606/10090?filter.stringency=stringent&limit=50000&page=1" \
-H  "accept: application/json" | \
jq -r \
'["fromGeneId", "fromSymbol", "toGeneId", "toSymbol", "best", 
"bestReverse", "methodCount", "totalMethodCount"], (.results[] | 
[.gene["id"], .gene["symbol"], .homologGene["id"], .homologGene["symbol"], 
.best, .bestReverse, .methodCount, .totalMethodCount]) 
| @tsv' > human2MouseOrthologs.tsv
```
- `mouseSpecificMets.tsv` and `mouseSpecificRxns.tsv`: The curated metabolic network that is not part of human metabolism but specific to mouse.
- `MGI_gene_ID_mapping.tsv`: The information of gene ID mapping between various identifiers that was extracted from the `MGI_Gene_Model_Coord.rpt` [file](http://www.informatics.jax.org/downloads/reports/index.html) in MGI.
- `MSigDB`: The mouse gene sets that were retrieved from the [Molecular Signatures Database](http://www.gsea-msigdb.org/gsea/index.jsp) and followed by replacing human genes with the ortholog counterparts of mouse.


