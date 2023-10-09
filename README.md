This is an R package developed to integrate analysis of multiple types of molecular profiling data (transcriptomics, proteomics, metabolomics, and epigenomics). It loads and deduplicates data, converts identifiers, performs enrichment/overrepresentation analysis on individual datasets, and combines these analyses with rank aggregation. There are multiple options for each step.

This software depends on the following R packages:

From CRAN:

	igraph
	httr
	RankAggreg
	RobustRankAggreg

From BioConductor:

	org.Mm.eg.db
	AnnotationDbi
	EnrichmentBrowser
	rWikiPathways
	PGSEA

It comes bundled with GSEA_4.0.3 and JDK11 for running GSEA.

To install:

git clone https://github.com/jamieobrett/multiomics.git

Then open R

	devtools::document(pkg="path/to/multiomics")
	devtools::install(pkg="path/to/multiomics")
	library(multiomics)
