#' Function: enrichGSEA
#' 
#' Runs and parses standard GSEA enrichment using the external GSEA java-based tool.
#' @param pathwaysDb KEGGPATHWAY, KEGGMODULE, BIOCYC, WIKIPATHWAYS, REACTOME
#' @param entityType GENE or COMPOUND
#' @param metadata Metadata dataframe, of the format returned by readMetadata
#' @param inputdata Expression values dataframe, of the format created by readData
#' @param outputDir Path to directory for output files to be written
#' @param minVariables Parameter for filter pathways (min number of variables in a pathway that pass minSamples and minLevel)
#' @param maxVariables Parameter for filter pathways (max number of variables in a pathway that pass minSamples and minLevel)
#' @param minSamples Parameter for naming the GCT file based on already-filtered data
#' @param minLevel Parameter for naming the GCT file based on already-filtered data
#' @param GSEA_weight Parameter for weighting the running enrichment score using p=1 ("weighted") or p=2 ("weighted_p2")
#' @param GSEA_sortType Parameter for ranking genes based on signed S2N metric ("real") or absolute value of the S2N metric ("abs")
#' @param GSEA_nperm Parameter for the number of permutations done in GSEA. Default = 1000 (recommended). Use 100 for testing.
#' @param GSEA_seed Parameter for seed set in GSEA (default = "timestamp").
#' @param group1 Experimental group name ("Old")
#' @param group2 Control group name ("Young")
#' @param doLog Take the log2 of abundance values prior to analysis?
#' @param DELETE FALSE (default) or TRUE for deleting all generated external files.
#' @param quieter If TRUE, will suppress printing output from GSEA.
#' @return Matrix with rows = pathways, columns = ("PATHID","PATHNAME","PATHSIZEPOSTFILTER","PATHSIZEPREFILTER","PVAL","QVAL","SCORE")
#' @export
enrichGSEA <- function(pathwaysDb=NULL,entityType=NULL,metadata=NULL,inputdata=NULL,outputDir=getwd(),group1=NULL,group2=NULL,minVariables=1,maxVariables=10000,GSEA_weight="weighted",GSEA_sortType="real",DELETE=FALSE,quieter=TRUE,GSEA_nperm=1000,GSEA_seed="timestamp",minSamples=NULL,minLevel=NULL,doLog=FALSE)
{
# Prepare files and tool
GSEABAT <- system.file("extdata","GSEAJAVA","GSEA_4.0.3","gsea-cli.bat",package="multiomics");
CLSFILE <- makeCLS(metadata=metadata,inputdata=inputdata,filepath=file.path(outputDir,paste0(attributes(inputdata)$datacategory,".cls")));
if(doLog) {
	inputdata[,colnames(inputdata)!="ID"] <- log2(inputdata[,colnames(inputdata)!="ID"]);
}
GCTFILE <- makeGCT(metadata=metadata,inputdata=inputdata,filepath=file.path(outputDir,paste0(attributes(inputdata)$datacategory,"_",attributes(inputdata)$idtype,"_","minSamp",minSamples,"minLevel",minLevel,(if(doLog) "log2" else "nolog"),".gct")));
GMTFILE <- system.file("extdata","GMTs",paste0(retrievePathwaysNameAndType(pathwaysDb=pathwaysDb,entityType=entityType)["name"],".gmt"),package="multiomics");
# Set up parameters
MODE_PARAM <- "GSEA";
RES_PARAM <- paste0("-res ",GCTFILE);
CLS_PARAM <- paste0("-cls ",CLSFILE,"#",group1,"_versus_",group2);
GMX_PARAM <- paste0("-gmx ",GMTFILE);
COLLAPSE_PARAM <- paste0("-collapse ","No_Collapse");
NORM_PARAM <- paste0("-norm ","meandiv");
NPERM_PARAM <- paste0("-nperm ",GSEA_nperm);
PERMUTE_PARAM <- paste0("-permute ","gene_set");
SCORE_PARAM <- paste0("-scoring_scheme ",GSEA_weight);
RPT_PARAM <- paste0("-rpt_label ",group1,"v",group2,"_",pathwaysDb,"_",attributes(inputdata)$datacategory); # Label for files
METRIC_PARAM <- paste0("-metric ","Signal2Noise"); #Other options: "tTest", "Ratio_of_Classes", "Diff_of_Classes"
SORT_PARAM <- paste0("-sort ",GSEA_sortType);
ORDER_PARAM <- paste0("-order ","descending");
OMIT_NOMATCHES_PARAM <- paste0("-include_only_symbols ","true");
MAKESETS_PARAM <- paste0("-make_sets ","false");
MEDIAN_PARAM <- paste0("-median ","false");
SEED_PARAM <- paste0("-rnd_seed ",GSEA_seed);
MAX_PARAM <- paste0("-set_max ",maxVariables);
MIN_PARAM <- paste0("-set_min ",minVariables);
OUT_PARAM <- paste0("-out ",outputDir);
MYARGS <- paste(c(MODE_PARAM,RES_PARAM,CLS_PARAM,GMX_PARAM,COLLAPSE_PARAM,NORM_PARAM,NPERM_PARAM,PERMUTE_PARAM,SCORE_PARAM,RPT_PARAM,METRIC_PARAM,SORT_PARAM,ORDER_PARAM,OMIT_NOMATCHES_PARAM,MAKESETS_PARAM,MEDIAN_PARAM,SEED_PARAM,MAX_PARAM,MIN_PARAM,OUT_PARAM),sep="",collapse=" ");
GSEACMD <- gsub("/","\\\\",GSEABAT);
MYARGS <- gsub("/","\\\\",MYARGS);
# Run GSEA
commandOut <- NULL;
if(format(Sys.time(),"%H%M")>2350) {Sys.sleep(900)} #Avoid confusion about how GSEA names files around midnight
commandOut <- NULL;
if(quieter) {
	try (commandOut <- system2(GSEACMD,MYARGS,stderr=FALSE,stdout=TRUE));
} else {
	try (commandOut <- system2(GSEACMD,MYARGS,stderr="",stdout=""));
}
# Process output
resultTable <- setNames(data.frame(matrix(ncol = 7, nrow = 0)), c("PATHID","PATHNAME","PATHSIZEPOSTFILTER","PATHSIZEPREFILTER","PVAL","QVAL","SCORE"));
if( length(attr(commandOut,"status"))>0 && attr(commandOut,"status")>0 ) {
	print(GSEACMD); # Something didn't work; debug
	stop("GSEA command did not run properly");
} else {
	outFolderContents <- list.files(outputDir,full.names=T);
	latestGSEA <- outFolderContents[match(max(file.info(outFolderContents)$ctime),file.info(outFolderContents)$ctime)];
	fileEnding <- strsplit(latestGSEA,".Gsea.")[[1]][-1];
	resultsFile1 <- file.path(latestGSEA,paste0("gsea_report_for_",group1,"_",fileEnding,".xls"));
	resultsFile2 <- file.path(latestGSEA,paste0("gsea_report_for_",group2,"_",fileEnding,".xls"));
	gseaResults <- read.delim(resultsFile1,stringsAsFactors=FALSE,header=TRUE);
	gseaResults <- rbind(gseaResults,read.delim(resultsFile2,stringsAsFactors=FALSE,header=TRUE));
	# Column names of interest in gseaResults: NAME, SIZE, NES, NOM.p.val, FDR.q.val
	resultTable[seq(1,dim(gseaResults)[1]),] <- NA;
	resultTable$PATHID <- gseaResults$NAME;
	if(pathwaysDb=="KEGGPATHWAY") { # GSEA will make pathway names in all caps, which interferes with future conversions only for KEGGPATHWAY IDs
		resultTable$PATHID <- gsub("MAP","map",resultTable$PATHID);
	}
	resultTable$PATHNAME <- convertPathwayInfo(pathwayIDs=resultTable$PATHID,type=pathwaysDb,outputInfo="NAME");
	resultTable$PATHSIZEPOSTFILTER <- gseaResults$SIZE;
	if(entityType=="GENE") {
		resultTable$PATHSIZEPREFILTER <- convertPathwayInfo(pathwayIDs=resultTable$PATHID,type=pathwaysDb,outputInfo="NGENES");
	} else if(entityType=="COMPOUND") {
		resultTable$PATHSIZEPREFILTER <- convertPathwayInfo(pathwayIDs=resultTable$PATHID,type=pathwaysDb,outputInfo="NCOMPOUNDS");
	}
	resultTable$PVAL <- gseaResults$NOM.p.val;
	resultTable$QVAL <- gseaResults$FDR.q.val;
	resultTable$SCORE <- gseaResults$NES;
	resultTable <- resultTable[!is.na(resultTable$PVAL),]; #https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/FAQ#What_does_it_mean_for_a_gene_set_to_have_NES_and_nominal_p-values_of_NaN_.28also_shown_as_blanks.29.3F For GSEA_sortType="abs" especially for GSEA_nperm with low values, there will not be any permutations that end up with "negatively enriched" (depleted of any changed genes) results for a denomintor. In general we don't care about these "negatively enriched" (unchanged) gene sets, so remove.
	}
	if (DELETE) { unlink(latestGSEA,recursive=T) }
return(resultTable);
}

#' Function: makeCLS
#'
#' Makes a .cls file for use in GSEA. If the file already exists, does not make a new one and prints a message indicating this.
#' @param metadata Metadata dataframe, of the format returned by readMetadata.
#' @param inputdata Data dataframe, of the format returned by readData.
#' @param filepath Full path to and including the file name to be generated.
#' @return The filepath is returned.
#' @export
makeCLS <- function(metadata=NULL,inputdata=NULL,filepath=NULL) {
if (file.exists(filepath)) {
	print(paste0("Using pre-existing file ",filepath));
	return(filepath);
}
sampleNames <- colnames(inputdata)[seq(2,length(colnames(inputdata)))];
metadata <- metadata[metadata$Name %in% sampleNames,]; #Subset the metadata to only the relevant samples
# Header
outputString <- as.character(dim(inputdata)[2]-1);
outputString <- paste(outputString,as.character(length(levels(metadata$Group))),"1\n",sep="\t");
# Classes list, note this must following the order of samples if uniquely taken
outputString <- paste0(outputString,"#\t");
sampleClassesVector <- metadata[match(sampleNames,metadata$Name),"Group"];
outputString <- paste0(outputString,paste(unique(sampleClassesVector),sep="",collapse="\t"),"\n");
# Classes for each sample
outputString <- paste0(outputString,paste(sampleClassesVector,sep="",collapse="\t"));
# Write
outfilehandle <- filepath;
writeLines(outputString, con=outfilehandle);
print(paste0("Wrote ",filepath));
return(filepath);
}

#' Function: makeGCT
#'
#' Makes a .gct file for use in GSEA. If the file already exists, does not make a new one and prints a message indicating this.
#' @param metadata Metadata dataframe, of the format returned by readMetadata.
#' @param inputdata Data dataframe, of the format returned by readData.
#' @param filepath Full path to and including the file name to be generated.
#' @return The filepath is returned.
#' @export
makeGCT <- function(metadata=NULL,inputdata=NULL,filepath=NULL) {
if (file.exists(filepath)) {
	print(paste0("Using pre-existing file ",filepath));
	return(filepath);
}
sampleNames <- colnames(inputdata)[seq(2,length(colnames(inputdata)))];
metadata <- metadata[metadata$Name %in% sampleNames,]; #Subset the metadata to only the relevant samples
# Header
outputString <- "#1.2\n";
outputString <- paste0(outputString,as.character(dim(inputdata)[1]),"\t",as.character(length(sampleNames)),"\n");
outputString <- paste0(outputString,"NAME\tDescription\t",paste(sampleNames,sep="",collapse="\t"));
# Values
for(i in 1:dim(inputdata)[1]) {
	outputString <- paste(c(outputString,paste(as.character(inputdata[i,]),sep="",collapse="\t")),sep="",collapse="\n");
}
# Write
outfilehandle <- filepath;
writeLines(outputString, con=outfilehandle);
print(paste0("Wrote ",filepath));
return(filepath);
}
