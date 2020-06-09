#' Function: overrepresentation
#' 
#' Performs overrepresentation analysis using Fisher's exact test.
#' @param pathwaysFilteredDf Prefiltered-for-size dataframe of the type returned by filterPathways.
#' @param pathwaysDb KEGGPATHWAY, KEGGMODULE, BIOCYC, WIKIPATHWAYS, REACTOME.
#' @param entityType GENE or COMPOUND
#' @param metadata Metadata dataframe, of the format returned by readMetadata
#' @param inputdata Binary values dataframe, of the format created by readData. Only one column per group (no replicates).
#' @param outputDir Path to directory for output files to be written
#' @param group1 Experimental group name ("Old")
#' @param group2 Control group name ("Young")
#' @param overrepType "Group1" for overrepresentation of Group1-unique genes. "Group2" for overrepresentation of Group2-unique genes. "Either" for overrepresentation among group-unique genes (genes unique to Group1 or Group2).
#' @param overrepScoreType "OR" for odds ratio (which may be Inf or 0; R can sort Inf but Inf cannot be used for aggregation weights), "lowCI" for the lower bound of the 95% confidence interval on the OR, or "OtoE" for the observed-to-expected ratio (which also may be Inf).
#' @return Matrix with rows = pathways, columns = ("PATHID","PATHNAME","PATHSIZEPOSTFILTER","PATHSIZEPREFILTER","PVAL","QVAL","SCORE"), with "SCORE" being the odds ratio (which can be Inf or 0).
#' @export
overrepresentation <- function(pathwaysFilteredDf=NULL,pathwaysDb=NULL,entityType=NULL,metadata=NULL,inputdata=NULL,group1=NULL,group2=NULL,overrepType=overrepType,overrepScoreType="OR")
{
# Prepare sample data
sampleNames <- colnames(inputdata)[seq(2,length(colnames(inputdata)))];
metadata <- metadata[metadata$Name %in% sampleNames,]; #Subset the metadata to only the relevant samples
group1name <- metadata[metadata$Group==group1,"Name"];
group2name <- metadata[metadata$Group==group2,"Name"];
uniqueTo1 <- inputdata$ID[inputdata[,group1name]==1 & inputdata[,group2name]==0];
uniqueTo2 <- inputdata$ID[inputdata[,group1name]==0 & inputdata[,group2name]==1];
uniqueTo1Or2 <- c(uniqueTo1,uniqueTo2);
all12 <- inputdata$ID;
# Determine overlaps between elements of interest (Group1-unique, Group2-unique, or Either) with pathway elements
myPathways <- unique(pathwaysFilteredDf$PATHWAY);
resultTable <- setNames(data.frame(matrix(ncol = 7, nrow = 0)), c("PATHID","PATHNAME","PATHSIZEPOSTFILTER","PATHSIZEPREFILTER","PVAL","QVAL","SCORE"));
resultTable[seq(1,length(myPathways)),] <- NA;
for(i in 1:length(myPathways)) {
	pathwayElements <- pathwaysFilteredDf$ENTITY[pathwaysFilteredDf$PATHWAY==myPathways[i]];
	if(overrepType=="Group1") {
		contingencyA <- length(intersect(uniqueTo1,pathwayElements));
		contingencyB <- length(uniqueTo1)-contingencyA;
		contingencyC <- length(pathwayElements)-contingencyA;
		contingencyD <- length(all12)-contingencyB-contingencyC-contingencyA;
	} else if(overrepType=="Group2") {
		contingencyA <- length(intersect(uniqueTo2,pathwayElements));
		contingencyB <- length(uniqueTo2)-contingencyA;
		contingencyC <- length(pathwayElements)-contingencyA;
		contingencyD <- length(all12)-contingencyB-contingencyC-contingencyA;
	} else { #"Either"
		contingencyA <- length(intersect(uniqueTo1Or2,pathwayElements));
		contingencyB <- length(uniqueTo1Or2)-contingencyA;
		contingencyC <- length(pathwayElements)-contingencyA;
		contingencyD <- length(all12)-contingencyB-contingencyC-contingencyA;
	}
	myTest <- fisher.test(matrix(c(contingencyA,contingencyB,contingencyC,contingencyD),2,2),alternative='greater');
	resultTable[i,"PVAL"] <- myTest$p.value;
	if (overrepScoreType=="OR") {
		resultTable[i,"SCORE"] <- myTest$estimate;
	} else if (overrepScoreType=="lowCI") {
		resultTable[i,"SCORE"] <- myTest$conf.int[1];
	} else {
		expectedOverlap <- (contingencyA + contingencyB)*length(pathwayElements)/length(all12);
		if(expectedOverlap==0) {
			resultTable[i,"SCORE"] <- Inf;
		} else {
			resultTable[i,"SCORE"] <- round(contingencyA/expectedOverlap,digits=3);
		}
	}
	resultTable[i,"PATHID"] <- myPathways[i];
	resultTable[i,"PATHSIZEPOSTFILTER"] <- length(pathwayElements);
}
resultTable$PATHNAME <- convertPathwayInfo(pathwayIDs=resultTable$PATHID,type=pathwaysDb,outputInfo="NAME");
if(entityType=="GENE") {
	resultTable$PATHSIZEPREFILTER <- convertPathwayInfo(pathwayIDs=resultTable$PATHID,type=pathwaysDb,outputInfo="NGENES");
} else { #COMPOUND
	resultTable$PATHSIZEPREFILTER <- convertPathwayInfo(pathwayIDs=resultTable$PATHID,type=pathwaysDb,outputInfo="NCOMPOUNDS");
}
resultTable$QVAL <- p.adjust(resultTable$PVAL,method="BH");
return(resultTable);
}