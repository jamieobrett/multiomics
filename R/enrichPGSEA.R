#' Function: enrichPGSEA
#' 
#' Runs and parses parametric GSEA enrichment using the PGSEA package.
#' @param pathwaysDb KEGGPATHWAY, KEGGMODULE, BIOCYC, WIKIPATHWAYS, REACTOME
#' @param entityType GENE or COMPOUND
#' @param metadata Metadata dataframe, of the format returned by readMetadata
#' @param inputdata Expression values dataframe, of the format created by readData
#' @param minVariables Parameter for filter pathways (min number of variables in a pathway that pass minSamples and minLevel)
#' @param maxVariables Parameter for filter pathways (max number of variables in a pathway that pass minSamples and minLevel)
#' @param group1 Experimental group name ("Old")
#' @param group2 Control group name ("Young")
#' @return Matrix with rows = pathways, columns = ("PATHID","PATHNAME","PATHSIZEPOSTFILTER","PATHSIZEPREFILTER","PVAL","QVAL","SCORE"). The SCORE is the Z-score associated with each p-value, which is converted into an FDR q-value by the Benjamini-Hochberg method.
#' @export
enrichPGSEA <- function(pathwaysDb=NULL,entityType=NULL,metadata=NULL,inputdata=NULL,group1=NULL,group2=NULL,minVariables=2,maxVariables=10000)
{
# Get GMT data
GMTFILE <- system.file("extdata","GMTs",paste0(retrievePathwaysNameAndType(pathwaysDb=pathwaysDb,entityType=entityType)["name"],".gmt"),package="multiomics");
GMT <- PGSEA::readGmt(GMTFILE);
# Prepare sample data. For input into PGSEA (which is designed for individual microarray analysis), take the median of each group's log2-transformed expression values and calculate the difference between groups (this is equivalent to taking the median of each pairwise comparison). For DNAm data, however, instead take the median of each group's absolute DNAm level and calculate the differences between groups (thus this is absolute differences, not FC of logged values).
sampleNames <- colnames(inputdata)[seq(2,length(colnames(inputdata)))];
metadata <- metadata[metadata$Name %in% sampleNames,]; #Subset the metadata to only the relevant samples
group1names <- metadata[metadata$Group==group1,"Name"];
group2names <- metadata[metadata$Group==group2,"Name"];
if(attributes(inputdata)$datacategory=="DNAm") {
	expressionData <- inputdata[,c(group1names,group2names)];
} else {
	expressionData <- log2(inputdata[,c(group1names,group2names)]);
}
rownames(expressionData) <- inputdata$ID;
fcData <- setNames(data.frame(matrix(ncol = 2, nrow = nrow(expressionData))), c(group1, group2));
rownames(fcData) <- rownames(expressionData);
for (i in 1:nrow(expressionData)) { 
	fcData[i,] <- tapply(as.numeric(expressionData[i,]),metadata$Group[match(colnames(expressionData),metadata$Name)],median);
}
fcData[,paste0(group2,"v",group1)] <- fcData[,group2] - fcData[,group1];
fcData[,c(group1,group2)] <- NULL;
# Run PGSEA
resultPGSEA <- PGSEA::PGSEA(fcData,cl=GMT,range=c(minVariables,maxVariables),ref=NULL,p.value=TRUE);
# Process results
resultTable <- setNames(data.frame(matrix(ncol = 7, nrow = 0)), c("PATHID","PATHNAME","PATHSIZEPOSTFILTER","PATHSIZEPREFILTER","PVAL","QVAL","SCORE"));
resultTable[seq(1,nrow(resultPGSEA$result)),] <- NA;
resultTable$PATHID <- trimws(rownames(resultPGSEA$result));
resultTable$PATHNAME <- convertPathwayInfo(pathwayIDs=resultTable$PATHID,type=pathwaysDb,outputInfo="NAME");
# Get and filter pathways to obtain postfilter size
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
sizes <- table(pathwaysFilteredDf$PATHWAY);
resultTable$PATHSIZEPOSTFILTER <- sizes[resultTable$PATHID];
if(entityType=="GENE") {
	resultTable$PATHSIZEPREFILTER <- convertPathwayInfo(pathwayIDs=resultTable$PATHID,type=pathwaysDb,outputInfo="NGENES");
} else if(entityType=="COMPOUND") {
	resultTable$PATHSIZEPREFILTER <- convertPathwayInfo(pathwayIDs=resultTable$PATHID,type=pathwaysDb,outputInfo="NCOMPOUNDS");
}
resultTable$PVAL <- resultPGSEA$p.result[,1];
resultTable$SCORE <- resultPGSEA$result[,1];
resultTable <- resultTable[!is.na(resultTable$PVAL),]; #PGSEA makes filtered-out pathways into NA rather than omitting them
resultTable$QVAL <- p.adjust(resultTable$PVAL,method="BH");
return(resultTable);
}