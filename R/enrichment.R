#' Function: indivEnrichment
#'
#' Calculates enrichment of pathways in a dataset with one type of data using various methods
#' @param analysistype GSEA, parametricGSEA, overrepresentation
#' @param pathwaysDb KEGGPATHWAY, KEGGMODULE, BIOCYC, WIKIPATHWAYS, REACTOME
#' @param entityType GENE or COMPOUND
#' @param metadata Metadata dataframe, of the format returned by readMetadata
#' @param inputdata Expression values dataframe, of the format returned by readData
#' @param group1 Experimental group name ("Old")
#' @param group2 Control group name ("Young")
#' @param outputDir Path to directory for output files to be written
#' @param minSamples Parameter for filtering data(minimum number of samples that have at least minLevel measured level)
#' @param minLevel Parameter for filtering data (minimum level for filtering)
#' @param minVariables Parameter for filter pathways (min number of variables in a pathway that pass minSamples and minLevel)
#' @param maxVariables Parameter for filter pathways (max number of variables in a pathway that pass minSamples and minLevel)
#' @param GSEA_weight Parameter for weighting the running enrichment score using p=1 ("weighted") or p=2 ("weighted_p2")
#' @param GSEA_sortType Parameter for ranking genes based on signed S2N metric ("real") or absolute value of the S2N metric ("abs")
#' @param GSEA_nperm Parameter for the number of permutations done in GSEA. Default = 1000 (recommended). Use 100 for testing.
#' @param GSEA_seed Parameter for seed set in GSEA (default = "timestamp").
#' @param overrepType Parameter for overrepresentation analysis. "Group1" for overrepresentation of Group1-unique genes. "Group2" for overrepresentation of Group2-unique genes. "Either" for overrepresentation among group-unique genes (genes unique to Group1 or Group2).
#' @param overrepScoreType Parameter for overrepresentation analysis. "OR" for odds ratio (which may be Inf or 0; R can sort Inf but Inf cannot be used for aggregation weights), "lowCI" for the lower bound of the 95% confidence interval on the OR, or "OtoE" for the observed-to-expected ratio (which also may be Inf).
#' @param doLog Take the log2 of abundance values prior to analysis? This only affects GSEA, since PGSEA always takes the fold-change of logged values.
#' @return Dataframe with rows = pathways, columns = ("PATHID","PATHNAME","PATHSIZEPOSTFILTER","PATHSIZEPREFILTER","PVAL","QVAL","SCORE"). Pathways are sorted in descending order of SCORE, with ties broken by QVAL. The dataframe has attributes "datacategory" equivalent to the input data and "pathwaysDb" equivalent to pathwaysDb.
#' @export
indivEnrichment <- function(analysistype="GSEA",pathwaysDb="KEGGPATHWAY",entityType="GENE",metadata=NULL,inputdata=NULL,group1=NULL,group2=NULL,outputDir=getwd(),minSamples=2,minLevel=5,minVariables=2,maxVariables=500,GSEA_weight="weighted",GSEA_sortType="real",overrepType="Either",GSEA_nperm=1000,GSEA_seed="timestamp",doLog=FALSE,overrepScoreType="OR")
{
print(paste("Starting an individual enrichment analysis with parameters:","analysistype",analysistype,"pathwaysDb",pathwaysDb,"entityType",entityType,"datacategory",attributes(inputdata)$datacategory,"group1",group1,"group2",group2,"minSamples",minSamples,"minLevel",minLevel,"minVariables",minVariables,"maxVariables",maxVariables,"GSEA_weight",GSEA_weight,"GSEA_sortType",GSEA_sortType,"overrepType",overrepType,"overrepScoreType",overrepScoreType,"GSEA_nperm",GSEA_nperm,"GSEA_seed",GSEA_seed,"doLog",doLog));
if (!file.exists(outputDir)) {dir.create(outputDir)};
resultTable <- NULL;
# Convert IDs of the input data to match those in the pathways
inputIdType <- attributes(inputdata)$idtype;
if(entityType=="GENE") {
	inputdata <- convertGeneIDs(tableWithIDCol=inputdata,inputID=inputIdType,outputID=retrievePathwaysNameAndType(pathwaysDb=pathwaysDb,entityType=entityType)["type"],quieter=FALSE);
} else if(entityType=="COMPOUND") {
	inputdata <- convertMetaboliteIDs(tableWithIDCol=inputdata,inputID=inputIdType,outputID=retrievePathwaysNameAndType(pathwaysDb=pathwaysDb,entityType=entityType)["type"],quieter=FALSE);
}
print(paste0("Converted identifiers from ",inputIdType," to ",attributes(inputdata)$idtype));
# Filter data
preFilterRows <- dim(inputdata)[1];
inputdata <- inputdata[rowSums(inputdata >= minLevel) >= minSamples,];
print(paste0("Filtered ",preFilterRows," ",entityType,"s based on minLevel ",minLevel," in minSamples ",minSamples,", resulting in ",dim(inputdata)[1], " final rows"));
# GSEA, which itself filters pathways
if(analysistype=="GSEA") {
	resultTable <- enrichGSEA(pathwaysDb=pathwaysDb,entityType=entityType,metadata=metadata,inputdata=inputdata,group1=group1,group2=group2,outputDir=outputDir,minVariables=minVariables,maxVariables=maxVariables,DELETE=FALSE,GSEA_weight=GSEA_weight,GSEA_sortType=GSEA_sortType,GSEA_nperm=GSEA_nperm,GSEA_seed=GSEA_seed,doLog=doLog,minSamples=minSamples,minLevel=minLevel);
} else if(analysistype=="parametricGSEA") {
	resultTable <- enrichPGSEA(pathwaysDb=pathwaysDb,entityType=entityType,metadata=metadata,inputdata=inputdata,group1=group1,group2=group2,minVariables=minVariables,maxVariables=maxVariables);
} else {
	# Get and filter pathways
	pathwaysName <- retrievePathwaysNameAndType(pathwaysDb=pathwaysDb,entityType=entityType)["name"];
	if(pathwaysDb=="KEGGPATHWAY" || pathwaysDb=="KEGGMODULE") {
		pathwaysDf <- makeKEGGpathwayTables(updateDesired=FALSE)[[pathwaysName]];
	} else if(pathwaysDb=="BIOCYC") {
		pathwaysDf <- makeBIOCYCpathwayTables(updateDesired=FALSE)[[pathwaysName]];
	} else if(pathwaysDb=="REACTOME") {
		pathwaysDf <- makeREACTOMEpathwayTables(updateDesired=FALSE)[[pathwaysName]];
	} else if(pathwaysDb=="WIKIPATHWAYS") {
		pathwaysDf <- makeWIKIPATHWAYSpathwayTables(updateDesired=FALSE)[[pathwaysName]];
	} else {
		stop(paste0("Could not retrieve that pathway dataframe for ",pathwaysName));
	}
	pathwaysFilteredDf <- filterPathways(minVariables=minVariables,maxVariables=maxVariables,inputPathways=pathwaysDf,variables=inputdata$ID);
	## Other enrichment types here
	if(analysistype=="overrepresentation") {
		resultTable <- overrepresentation(pathwaysFilteredDf=pathwaysFilteredDf,pathwaysDb=pathwaysDb,entityType=entityType,metadata=metadata,inputdata=inputdata,group1=group1,group2=group2,overrepType=overrepType,overrepScoreType=overrepScoreType);
	}
}
resultTable <- resultTable[order(resultTable$SCORE,(-1)*resultTable$QVAL,decreasing=TRUE),];
rownames(resultTable) <- 1:nrow(resultTable);
attributes(resultTable)$datacategory <- attributes(inputdata)$datacategory;
attributes(resultTable)$pathwaysDb <- pathwaysDb;
return(resultTable);
}

#' Function: retrievePathwaysNameAndType
#'
#' Helper function for enrichment tools. Obtains the name of the datatable appropriate for the pathway database type and the entity type.
#' @param pathwaysDb KEGGPATHWAY, KEGGMODULE, BIOCYC, WIKIPATHWAYS, REACTOME
#' @param entityType GENE or COMPOUND
#' @return Named vector of strings: "name" to use when obtaining the dataframe from a named list or from extdata (example: "GENE_KEGG-KEGGMODULE"), "type" is the type of entity ID (examples: "KEGG", "MGI").
#' @export
retrievePathwaysNameAndType <- function(pathwaysDb=NULL,entityType=NULL) {
if(pathwaysDb=="KEGGPATHWAY" && entityType=="GENE") {
	result <- c("GENE_KEGG-KEGGPATHWAY","KEGG");
} else if(pathwaysDb=="KEGGPATHWAY" && entityType=="COMPOUND") {
	result <- c("METABOLITE_KEGG-KEGGPATHWAY","KEGG");
} else if(pathwaysDb=="KEGGMODULE" && entityType=="GENE") {
	result <- c("GENE_KEGG-KEGGMODULE","KEGG");
} else if(pathwaysDb=="KEGGMODULE" && entityType=="COMPOUND") {
	result <- c("METABOLITE_KEGG-KEGGMODULE","KEGG");
} else if(pathwaysDb=="BIOCYC" && entityType=="GENE") {
	result <- c("GENE_MGI-BIOCYCPATHWAY","MGI");
} else if(pathwaysDb=="BIOCYC" && entityType=="COMPOUND") {
	result <- c("METABOLITE_BIOCYC-BIOCYCPATHWAY","BIOCYC");
} else if(pathwaysDb=="WIKIPATHWAYS" && entityType=="GENE") {
	result <- c("GENE_ENTREZID-WIKIPATHWAYS","ENTREZID");
} else if(pathwaysDb=="WIKIPATHWAYS" && entityType=="COMPOUND") {
	result <- c("METABOLITE_KEGG-WIKIPATHWAYS","KEGG");
} else if(pathwaysDb=="REACTOME" && entityType=="GENE") {
	result <- c("GENE_ENTREZID-REACTOMEPATHWAY","ENTREZID");
} else if(pathwaysDb=="REACTOME" && entityType=="COMPOUND") {
	result <- c("METABOLITE_CHEBI-REACTOMEPATHWAY","CHEBI");
} else {
	stop(paste0("Could not retrieve the pathways name for pathwaysDb ",pathwaysDb," and entityType ",entityType));
}
names(result) <- c("name","type");
return(result);
}

#' Function: filterPathways
#'
#' Filters input pathways for representation in the data. The gene/metabolite identifier type of the pathway dataframe and variables vector must be the same (example: both ENTREZID).
#' @param minVariables Min number of variables in a pathway that must be present
#' @param maxVariables Max number of variables in a pathway that can be present
#' @param inputPathways Dataframe of the format in the lists returned by make[DB]pathwayTables (with rows ENTITY and PATHWAY)
#' @param variables Vector of variables (that are considered expressed in a dataset) to check for representation in pathways. For example, the ID column of a loaded dataframe after filtering for expressing
#' @return Dataframe identical to inputPathways except that the pathways that did not meet criteria have been removed, and rows where the ENTITY is not in the input variables vector have been removed.
#' @export
filterPathways <- function(minVariables=2,maxVariables=500,inputPathways=NULL,variables=NULL)
{
filteredPathways1 <- inputPathways[inputPathways$ENTITY %in% variables,]; #Remove all pathways with 0 variables from the variables vector
if(dim(filteredPathways1)[1]==0) { stop("All pathways were filtered out. Check min-max parameters and ID conversions."); }
representation <- colSums(table(filteredPathways1)); #Vector with names = pathways, elements = counts of variables represented in each pathway
pathwaysMeetingCriteria <- names(representation[representation >= minVariables & representation <= maxVariables]);
result <- filteredPathways1[filteredPathways1$PATHWAY %in% pathwaysMeetingCriteria,];
print(paste0("Filtered from ",length(unique(inputPathways$PATHWAY))," pre-filter input pathways to ",length(unique(result$PATHWAY))," post-filter pathways based on ",minVariables," minVariables and ",maxVariables," maxVariables"));
return(result);
}