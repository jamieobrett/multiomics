#' Function: makeKEGGpathwayTables
#'
#' Gets pathway, module, gene, and compound data from the KEGG API. Pathways named "mmu00010", "map00010", or "ko00010" are synonymous for these purposes; the number after the prefix refers to the same pathway for organism-specific ("mmu") and universal ("map" or "ko"). The same for modules such as "mmu_M00001" and "M00001".
#' @param updateDesired TRUE or FALSE. If TRUE, will regenerate the tables from a new web GET; each new csv file will be written to the cwd (it is up to the developer to move this to ./metabolomics/inst/extdata and overwrite prior tables if desired). If FALSE, will use existing tables in that folder.
#' @return Named list of dataframes: "GENE_KEGG-KEGGPATHWAY", "GENE_KEGG-KEGGMODULE", "METABOLITE_KEGG-KEGGPATHWAY", "METABOLITE_KEGG-KEGGMODULE". Each dataframe contains two columns titled "ENTITY" and "PATHWAY". KEGG identifiers are of the form org:idnumber for genes ("mmu:103988") or idnumber for metabolites ("C00022"). Pathways are named with the universal KEGG Map prefix ("map00010"); modules are named with the universal KEGG Module prefix ("M00001"). Each dataframe has an attribute "tablename" matching the name used in the named list.
#' @export
makeKEGGpathwayTables <- function(updateDesired=FALSE)
{
result <- list();
URLS <- list(
	"GENE_KEGG-KEGGPATHWAY"="http://rest.kegg.jp/link/pathway/mmu",
	"GENE_KEGG-KEGGMODULE"="http://rest.kegg.jp/link/module/mmu",
	"METABOLITE_KEGG-KEGGPATHWAY"="http://rest.kegg.jp/link/pathway/cpd",
	"METABOLITE_KEGG-KEGGMODULE"="http://rest.kegg.jp/link/module/cpd"
);
if(updateDesired) {
	for(i in 1:length(URLS)) {
		getRowed <- getReadHttp(txtURL=URLS[[i]]);
		# Parse output into a dataframe. Input content has each entry as follows:
		# GENE-KEGGPATHWAY = "mmu:103988", "\t", "path:mmu00010"
		# GENE-KEGGMODULE = "mmu:103988", "\t", "md:mmu_M00001"
		# METABOLITE-KEGGPATHWAY = "cpd:C00022", "\t", "path:map00010"
		# METABOLITE-KEGGMODULE = "cpd:C00022", "\t", "md:M00001"
		resultTable <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("ENTITY", "PATHWAY"));
		for(j in 1:length(getRowed)) {
			thisRow <- strsplit(getRowed[j],"\t")[[1]];
			resultTableEntity <- gsub("cpd:","",thisRow[1]); # Remove "cpd:"
			resultTablePathway <- gsub("path:mmu|path:map","map",gsub("md:mmu_|md:","",thisRow[2]) ); # Remove "md:mmu_" or "md:" from modules and give pathways the "map" prefix
			resultTable[j,] <- c(resultTableEntity, resultTablePathway);
		}
		# Add to the final list of dataframes
		result[[names(URLS[i])]] <- resultTable;
		# Update the written tables
		write.csv(resultTable,file=paste0(names(URLS[i]),".csv"),row.names=FALSE);
	}
} else {
	result <- readTables(tableNames=names(URLS),tableType="csv",extdataFolder="PATHWAYs");
}
for(i in 1:length(result)) {
	attributes(result[[i]])$tablename <- names(result)[i];
}
return(result);
}

#' Function: makeKEGGpathwaynameTables
#'
#' Gets pathway and module names from the KEGG API.
#' @param updateDesired TRUE or FALSE. If TRUE, will regenerate the table from a new web GET; each new csv file will be written to the cwd (it is up to the developer to move this to ./metabolomics/inst/extdata and overwrite prior tables if desired). If FALSE, will use existing tables in that folder.
#' @return Named list dataframes "KEGGPATHWAYPATHWAYNAMES" and "KEGGMODULEPATHWAYNAMES" containing two columns titled "PATHWAY" and "NAME". Pathways are named with the universal KEGG Map prefix ("map00010"); modules are named with the universal KEGG Module prefix ("M00001").
#' @export
makeKEGGpathwaynameTables <- function(updateDesired=FALSE)
{
result <- list();
URLS <- list(
	"KEGGPATHWAYPATHWAYNAMES"="http://rest.kegg.jp/list/pathway",
	"KEGGMODULEPATHWAYNAMES"="http://rest.kegg.jp/list/module"
);
if(updateDesired) {
	for(i in 1:length(URLS)) {
		getRowed <- getReadHttp(txtURL=URLS[[i]]);
		resultTable <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("PATHWAY", "NAME"));
		# Parse output into a dataframe. Input content has each entry as follows:
		# "path:mmu00010", "\t", "Glycolysis / Gluconeogenesis - Mus musculus (mouse)"
		for(j in 1:length(getRowed)) {
			thisRow <- strsplit(getRowed[j],"\t")[[1]];
			resultTablePathway <- gsub("path:mmu|path:map","map",gsub("md:mmu_|md:","",thisRow[1])); # Remove "md:mmu_" or "md:" from modules and give pathways the "map" prefix
			resultTableName <- thisRow[2];
			resultTable[dim(resultTable)[1]+1,] <- c(resultTablePathway, resultTableName);
		}
	# Update the written tables
	write.table(resultTable,file=paste0(names(URLS[i]),".txt"),row.names=FALSE,quote=FALSE,sep="\t");
	# Add to the final list of dataframes
	result[[names(URLS[i])]] <- resultTable;
	}
} else {
	result <- readTables(tableNames=names(URLS),tableType="txt",extdataFolder="PATHWAYNAMEs");
}
return(result);
}

#' Function: makeREACTOMEpathwayTables
#'
#' Gets pathway, gene, and compound data from the Reactome website (https://reactome.org/download-data). Pathways named "mmu00010", "map00010", or "ko00010" are synonymous for these purposes; the number after the prefix refers to the same pathway for organism-specific("mmu") and universal ("map" or "ko"). The same for modules such as "mmu_M00001" and "M00001".
#' @param updateDesired TRUE or FALSE. If TRUE, will regenerate the tables from a new web GET; each new csv file will be written to the cwd (it is up to the developer to move this to ./metabolomics/inst/extdata and overwrite prior tables if desired). If FALSE, will use existing tables in that folder.
#' @return Named list of dataframes: "GENE_ENTREZID-REACTOMEPATHWAY", "METABOLITE_CHEBI-REACTOMEPATHWAY". Each dataframe contains two columns titled "ENTITY" and "PATHWAY". Genes are identified by Entrez ID. Metabolites are identified by ChEBI ID. Reactome pathway identifiers are the Reactome Stable Identifier ("R-MMU-166663"). Each dataframe has an attribute "tablename" matching the name used in the named list.
#' @export
makeREACTOMEpathwayTables <- function(updateDesired=FALSE) {
result <- list();
URLS <- list(
	"GENE_ENTREZID-REACTOMEPATHWAY"="https://reactome.org/download/current/NCBI2Reactome_All_Levels.txt",
	"METABOLITE_CHEBI-REACTOMEPATHWAY"="https://reactome.org/download/current/ChEBI2Reactome_All_Levels.txt"
);
if(updateDesired) {
	for(i in 1:length(URLS)) {
		getRowed <- getReadHttp(txtURL=URLS[[i]]);
		# Parse output into a dataframe. Input content has each entry as follows: Tab-delimited: gene/compound ID, Reactome pathway ID, URL, pathway/reaction name, evidence code, species
		resultTable <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("ENTITY", "PATHWAY"));
		for(j in 1:length(getRowed)) {
			thisRow <- strsplit(getRowed[j],"\t")[[1]];
			if(thisRow[6] != "Mus musculus") { next; } #Filter for mouse pathways
			resultTableEntity <- thisRow[1]; #Gene (Entrez) or compound (ChEBI)
			resultTablePathway <- thisRow[2]; #Reactome stable ID pathway
			resultTable[dim(resultTable)[1]+1,] <- c(resultTableEntity, resultTablePathway);
		}
		# Add to the final list of dataframes
		result[[names(URLS[i])]] <- resultTable;
		# Update the written tables
		write.csv(resultTable,file=paste(names(URLS[i]),".csv",sep=""),row.names=FALSE);
	}
} else {
	result <- readTables(tableNames=names(URLS),tableType="csv",extdataFolder="PATHWAYs");
}
for(i in 1:length(result)) {
	attributes(result[[i]])$tablename <- names(result)[i];
}
return(result);
}

#' Function: makeREACTOMEpathwaynameTable
#'
#' Gets pathway names from the Reactome website.
#' @param updateDesired TRUE or FALSE. If TRUE, will regenerate the table from a new web GET; each new csv file will be written to the cwd (it is up to the developer to move this to ./metabolomics/inst/extdata and overwrite prior tables if desired). If FALSE, will use existing tables in that folder.
#' @return Named list of single dataframe "REACTOMEPATHWAYNAMES" containing two columns titled "PATHWAY" and "NAME". Reactome pathway identifiers are the Reactome Stable Identifier ("R-MMU-166663").
#' @export
makeREACTOMEpathwaynameTable <- function(updateDesired=FALSE)
{
TABLENAME <- "REACTOMEPATHWAYNAMES";
if(updateDesired) {
	REACTOME_URL <- "https://reactome.org/download/current/ReactomePathways.txt";
	resultTable <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("PATHWAY", "NAME"));
	getRowed <- getReadHttp(txtURL=REACTOME_URL);
	# Parse output into a dataframe. Input content has each entry as follows:
	# "R-BTA-73843", "\t", "5-Phosphoribose 1-diphosphate biosynthesis", "\t", "Bos taurus"
	for(j in 1:length(getRowed)) {
		thisRow <- strsplit(getRowed[j],"\t")[[1]];
		if(thisRow[3] != "Mus musculus") { next; } #Filter for mouse pathways
		resultTablePathway <- thisRow[1];
		resultTableName <- thisRow[2];
		resultTable[dim(resultTable)[1]+1,] <- c(resultTablePathway, resultTableName);
	}
	# Update the written tables
	write.table(resultTable,file=paste0(TABLENAME,".txt"),row.names=FALSE,quote=FALSE,sep="\t");
	result <- list(resultTable); names(result) <- TABLENAME;
} else {
	result <- readTables(tableNames=c(TABLENAME),tableType="txt",extdataFolder="PATHWAYNAMEs");
}
return(result);
}

#' Function: makeBIOCYCpathwayTables
#'
#' Gets pathway, gene, and compound data from BioCyc based on pre-created Smart Tables.
#' @param updateDesired TRUE or FALSE. If TRUE, will regenerate the tables from a new web GET; each new csv file will be written to the cwd (it is up to the developer to move this to ./metabolomics/inst/extdata and overwrite prior tables if desired). If FALSE, will use existing tables in that folder.
#' @return Named list of dataframes: "GENE_MGI-BIOCYCPATHWAY", "METABOLITE_BIOCYC-BIOCYCPATHWAY". Each dataframe contains two columns titled "ENTITY" and "PATHWAY". Pathway names are the BioCyc Object ID ("PWY-5657" or "GLYCOGENINS"). Genes names are the MGI ID ("MGI:97362"). Metabolite names are the BioCyc Object ID ("L-CYSTATHIONINE"). Each dataframe has an attribute "tablename" matching the name used in the named list.
#' @export
makeBIOCYCpathwayTables <- function(updateDesired=FALSE)
{
GENETABLENAME <- "GENE_MGI-BIOCYCPATHWAY";
METABTABLENAME <- "METABOLITE_BIOCYC-BIOCYCPATHWAY";
resultTableGenes <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("ENTITY", "PATHWAY"));
resultTableMetabolites <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("ENTITY", "PATHWAY"));
if(updateDesired) {
	BioCycURL <- "https://websvc.biocyc.org/st-service-get?id=biocyc17-40784-3786381030&format=tsv";
	BioCycData <- getReadHttp(txtURL=BioCycURL);
	# Each element of this vector is a pathway: pathwayid, "\t", MGI:gene1 // MGI:gene2 // ... // MGI:geneX, "\t", CPD1 // CPD2 // ... // CPDX, "\t", pathway common name
	for(k in 2:length(BioCycData)) { # Skip the header
		mySplit <- strsplit(BioCycData[k],"\t")[[1]]; # pathwayid, '//' delimited list of genes, '//' delimited list of compounds, common name of pathway
		myPathway <- mySplit[1];
		myGenes <- strsplit(mySplit[2]," // ")[[1]];
		myMetabolites <- strsplit(mySplit[3]," // ")[[1]];
		resultTableGenes <- rbind(resultTableGenes, data.frame("ENTITY"=myGenes,"PATHWAY"=myPathway));
		resultTableMetabolites <- rbind(resultTableMetabolites, data.frame("ENTITY"=myMetabolites, "PATHWAY"=myPathway));
	}
	# Update the written tables
	write.csv(resultTableGenes,file=paste0(GENETABLENAME,".csv"),row.names=FALSE);
	write.csv(resultTableMetabolites,file=paste0(METABTABLENAME,".csv"),row.names=FALSE);
	result <- list(resultTableGenes, resultTableMetabolites); names(result) <- c(GENETABLENAME, METABTABLENAME);
} else {
	result <- readTables(tableNames=c(GENETABLENAME,METABTABLENAME),tableType="csv",extdataFolder="PATHWAYs");
}
for(i in 1:length(result)) {
	attributes(result[[i]])$tablename <- names(result)[i];
}
return(result);
}

#' Function: makeBIOCYCpathwaynameTable
#'
#' Gets pathway names from the BioCyc website.
#' @param updateDesired TRUE or FALSE. If TRUE, will regenerate the table from a new web GET; each new csv file will be written to the cwd (it is up to the developer to move this to ./metabolomics/inst/extdata and overwrite prior tables if desired). If FALSE, will use existing tables in that folder.
#' @return Named list of single dataframe "BIOCYCPATHWAYNAMES" containing two columns titled "PATHWAY" and "NAME". Pathway names are the BioCyc Object ID ("PWY-5657" or "GLYCOGENINS"). 
#' @export
makeBIOCYCpathwaynameTable <- function(updateDesired=FALSE)
{
TABLENAME <- "BIOCYCPATHWAYNAMES";
if(updateDesired) {
	BioCycURL <- "https://websvc.biocyc.org/st-service-get?id=biocyc17-40784-3786381030&format=tsv";
	resultTable <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("PATHWAY", "NAME"));
	BioCycData <- getReadHttp(txtURL=BioCycURL);
	# Each element of this vector is a pathway: pathwayid, "\t", MGI:gene1 // MGI:gene2 // ... // MGI:geneX, "\t", CPD1 // CPD2 // ... // CPDX, "\t", quoted "pathway common name"
	for(k in 2:length(BioCycData)) { # Skip the header
		mySplit <- strsplit(BioCycData[k],"\t")[[1]]; # pathwayid, '//' delimited list of genes, '//' delimited list of compounds, common name of pathway
		resultTablePathway <- mySplit[1];
		resultTableName <- gsub("^\"|\"$","",mySplit[4]); # Remove quotes
		resultTable[dim(resultTable)[1]+1,] <- c(resultTablePathway, resultTableName);
	}
	# Update the written tables
	write.table(resultTable,file=paste0(TABLENAME,".txt"),row.names=FALSE,quote=FALSE,sep="\t");
	result <- list(resultTable); names(result) <- TABLENAME;
} else {
	result <- readTables(tableNames=c(TABLENAME),tableType="txt",extdataFolder="PATHWAYNAMEs");
}
return(result);
}

#' Function: makeWIKIPATHWAYSpathwayTables
#'
#' Gets pathway, gene, and compound data from WikiPathways, using rWikiPathways.
#' @param updateDesired TRUE or FALSE. If TRUE, will regenerate the tables from a new web queries; each new csv file will be written to the cwd (it is up to the developer to move this to ./metabolomics/inst/extdata and overwrite prior tables if desired). If FALSE, will use existing tables in that folder.
#' @param CONVERSION_CUTOFF Default 0.6. At least CONVERSION_CUTOFF metabolites in each pathway must be converted to use the pathway.
#' @return Named list of dataframes: "GENE_ENTREZID-WIKIPATHWAYS", "METABOLITE_KEGG-WIKIPATHWAYS". Each dataframe contains two columns titled "ENTITY" and "PATHWAY". Pathway names are the WikiPathways ID ("WP1"). Genes names are the Entrez ID ("23887"), and metabolite names are the KEGG ID ("C00037"). Any metabolites with IDs that could not be converted to these forms are omitted. Each dataframe has an attribute "tablename" matching the name used in the named list.
#' @export
makeWIKIPATHWAYSpathwayTables <- function(updateDesired=FALSE,CONVERSION_CUTOFF=0.6)
{
GENETABLENAME <- "GENE_ENTREZID-WIKIPATHWAYS";
METABTABLENAME <- "METABOLITE_KEGG-WIKIPATHWAYS";
resultTableGenes <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("ENTITY", "PATHWAY"));
resultTableMetabolites <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("ENTITY", "PATHWAY"));
if(updateDesired) {
	# METABOLITES
	wikipathways <- rWikiPathways::listPathwayIds(organism="Mus musculus");
	badPathways <- vector();
	for(i in 1:length(wikipathways)) {
		gpml <- rWikiPathways::getPathway(wikipathways[i]); # gpml format of pathway
		gpmlLines <- strsplit(gpml, "\n")[[1]];
		# exclude pathways where not enough elements convert for use
		totalMetabolites <- 0;
		convertedMetabolites <- 0;
		# Parse gpml
		j <- 1;
		while(j <= length(gpmlLines)) {
			if(grepl("<DataNode",gpmlLines[j])) { # gene or metabolite
				if(grepl("Type=\"Metabolite",gpmlLines[j])) { # metabolite
					while(!grepl("</DataNode>",gpmlLines[j])) { # find the metabolite info
						j <- j+1;
						if(grepl("<Xref Database=\".+\" ID=\".+\" />",gpmlLines[j])) { # there's an ID entry
							totalMetabolites <- totalMetabolites+1;
							preConversionSource <- trimws(gsub("<Xref Database=\"(.+)\" ID=\"(.+)\" />","\\1",gpmlLines[j])); # get the ID database
							preConversionID <- trimws(gsub("<Xref Database=\"(.+)\" ID=\"(.+)\" />","\\2",gpmlLines[j])); # get the ID
							convertible <- TRUE;
							if(preConversionSource=="KEGG Compound") {
								preConversionSource <- "KEGG";
							} else if(preConversionSource=="PubChem-compound") {
								preConversionSource <- "PUBCHEM_CID";
							} else if(preConversionSource=="HMDB") {
								preConversionSource <- "HMDB";
							} else if(preConversionSource=="ChEBI") {
								preConversionSource <- "CHEBI";
								preConversionID <- strsplit(preConversionID,":")[[1]][2];
							}
							else if(preConversionSource=="CAS") {
								preConversionSource <- "CAS";
							}
							else if(preConversionSource=="LIPID MAPS") {
								preConversionSource <- "LIPIDMAPS";
							} else {
								print(paste0(preConversionSource," type of metabolite ID and ",preConversionID," ID cannot be converted for pathway ",wikipathways[i]));
								convertible <- FALSE;
							}
							if(convertible) {
								conversionRes <- NULL;
								conversionRes <- convertMetaboliteIDs(tableWithIDCol=data.frame("ID"=preConversionID,"dummy"=0),inputID=preConversionSource,outputID="KEGG",quieter=TRUE); # convert the ID
								if(!is.null(conversionRes)) {
									postConversionID <- as.character(conversionRes$ID[1]);
									resultTableMetabolites[dim(resultTableMetabolites)[1]+1,] <- c(postConversionID,wikipathways[i]); # add to the table
									convertedMetabolites <- convertedMetabolites+1;
								} else {
									print(paste0(preConversionSource," type of metabolite ID and ",preConversionID," ID cannot be converted for pathway ",wikipathways[i]));
								}
							}
							break;
						}
					}
				}
			}
			j <- j+1;
		}
		if(totalMetabolites>0 && convertedMetabolites/totalMetabolites<CONVERSION_CUTOFF) {
			badPathways <- c(badPathways,wikipathways[i]);
		}
	}
	# Remove pathways that didn't mean CONVERSION_CUTOFF
	print(paste0("Poorly converted pathways: ",paste0(badPathways,collapse=" ")));
	resultTableMetabolites <- resultTableMetabolites[!(resultTableMetabolites$PATHWAY %in% badPathways),]
	# Remove any metabolite NA identifiers
	idfound <- !is.na(resultTableMetabolites$ENTITY);
	resultTableMetabolites <- resultTableMetabolites[idfound,];
	# Remove duplicates
	resultTableMetabolites <- resultTableMetabolites[!duplicated(resultTableMetabolites),];
	# GENES
	wikigmtfile <- rWikiPathways::downloadPathwayArchive(organism="Mus musculus",format="gmt");
	wikigmt <- scan(file=wikigmtfile,what="",sep="\n",quote=NULL); #First column: pathway name%WikiPathways_date%WPnnn%Mus musculus. Second column: url. Third+ columns: Entrez IDs.
	for(i in 1:length(wikigmt)) {
		geneset <- strsplit(wikigmt[i],"\t")[[1]];
		genesetID <- gsub(".+%(WP[[:digit:]]+)%Mus musculus$","\\1",geneset[1]); #WPnnn
		genesetgenes <- geneset[seq(3,length(geneset))];
		for(j in 1:length(genesetgenes)) {
			resultTableGenes[dim(resultTableGenes)[1]+1,] <- c(genesetgenes[j],genesetID);
		}
	}
	# Update the written tables
	write.csv(resultTableGenes,file=paste0(GENETABLENAME,".csv"),row.names=FALSE);
	write.csv(resultTableMetabolites,file=paste0(METABTABLENAME,".csv"),row.names=FALSE);
	result <- list(resultTableGenes, resultTableMetabolites); names(result) <- c(GENETABLENAME, METABTABLENAME);
} else {
	result <- readTables(tableNames=c(GENETABLENAME,METABTABLENAME),tableType="csv",extdataFolder="PATHWAYs");
}
for(i in 1:length(result)) {
	attributes(result[[i]])$tablename <- names(result)[i];
}
return(result);
}

#' Function: makeWIKIPATHWAYSpathwaynameTable
#'
#' Gets pathway names using rWikiPathways.
#' @param updateDesired TRUE or FALSE. If TRUE, will regenerate the tables from a new web queries; each new csv file will be written to the cwd (it is up to the developer to move this to ./metabolomics/inst/extdata and overwrite prior tables if desired). If FALSE, will use existing tables in that folder.
#' @return Named list of single dataframe "WIKIPATHWAYSPATHWAYNAMES" containing two columns titled "PATHWAY" and "NAME". Reactome pathway identifiers are the Reactome Stable Identifier ("R-MMU-166663").
#' @export
makeWIKIPATHWAYSpathwaynameTable <- function(updateDesired=FALSE)
{
TABLENAME <- "WIKIPATHWAYSPATHWAYNAMES";
wikipathways <- rWikiPathways::listPathwayIds(organism="Mus musculus");
if(updateDesired) {
	resultTable <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("PATHWAY", "NAME"));
	pathways <- rWikiPathways::listPathways(organism="Mus musculus");
	for(i in 1:length(pathways)) {
		resultTablePathway <- pathways[[i]][["id"]];
		resultTableName <- pathways[[i]][["name"]];
		resultTable[dim(resultTable)[1]+1,] <- c(resultTablePathway, resultTableName);
	}
	# Update the written tables
	write.table(resultTable,file=paste0(TABLENAME,".txt"),row.names=FALSE,quote=FALSE,sep="\t");
	result <- list(resultTable); names(result) <- TABLENAME;
} else {
	result <- readTables(tableNames=c(TABLENAME),tableType="txt",extdataFolder="PATHWAYNAMEs");
}
return(result);
}

#' Function: tableToGMT
#'
#' Generates gmt tables for GSEA use from pathway tables. Each new file will be written to the cwd (it is up to the developer to move this to ./metabolomics/inst/extdata and overwrite prior tables if desired).
#' @param pathwayTable Dataframe to be converted into a gmt. The first column should have entity (gene, metabolite) identifiers. The second column should have pathway identifiers. The attribute "tablename" will be used at the prefix for the file.
#' @export
tableToGMT <- function(pathwayTable=NULL)
{
pathways <- NULL;
gmtprefix <- attributes(pathwayTable)$tablename;
pathways[levels(factor(pathwayTable$PATHWAY))] <- list("NA"); #List with each element named as the pathway
for(i in 1:dim(pathwayTable)[1]) {
	pathways[[pathwayTable$PATHWAY[i]]] <- c(pathways[[pathwayTable$PATHWAY[i]]], pathwayTable$ENTITY[i]);
}
outputString <- NULL;
for(i in 1:length(pathways)) {
	outputString <- paste(c(outputString, paste(c(names(pathways)[i], pathways[[i]]),sep="",collapse="\t")), sep="", collapse="\n");
}
outfilehandle <- file(paste0(gmtPrefix,".gmt"));
writeLines(outputString, con=outfilehandle);
close(outfilehandle);
}