#' Function: readMetadata
#'
#' Reads in the metadata file.
#' @param metadatafile Full path and filename of the metadata, which must be a tab-delimited text file. \cr
#'  row 1: headers of the columns. \cr
#'  column Name: UNIQUE sample names (example: YQ1), which must match the sample names in the data files exactly. \cr
#'  column Group: group name of each sample (example: Young). \cr
#'  column Color: color of each sample for plotting (example: blue). \cr
#'  column DataCategory: category of data (examples: RNA, Protein, Metabolite, DNAm, NFR). \cr
#'  column Remove: 1 if outlier or test sample to be removed from all analysis.
#' @return Dataframe with columns: Name, Group, Color, DataCategory, Remove.
#' @export
readMetadata <- function(metadatafile="metadata.txt")
{
result <- read.delim(metadatafile,stringsAsFactors=FALSE);
if(anyDuplicated(result$Name)) { stop("Sample names are not unique in the metadata file."); } # Check for unique sample names
result$Group <- factor(result$Group);
print("Done reading in metadata");
return(result);
}

#' Function: readFilelist
#'
#' Reads in a list of input datafiles.
#' @param filelistfile Full path and filename of the filelist, which must be a tab-delimited text file. \cr
#'  row 1: headers of the columns. \cr
#'  column DataCategory: category of data (examples: RNA, Protein, Metabolite, DNAm, NFR). \cr
#'  column Path: Full path and filename of the input file. \cr
#'  column Format: Either "absolute" or "fraction" or "binary" [absolute: RNA-Seq read counts, mass spec intensities; fraction: DNA methylation ratios; binary: NFR presence, variable beyond cutoff threshold]. \cr
#'  column Type: "GENE" or "COMPOUND". \cr
#'  column IDType: One of "BIOCYC", "KEGG", "HMDB", "PUBCHEM_SID", "PUBCHEM_CID", "CHEBI", "METLIN", "NAME", "COMPOUND", "ENTREZID", "ENTREZPROTEIN", "UNIPROT", "SYMBOL".
#' @return Dataframe with columns: File, Format, Type, IDType.
#' @export
readFilelist <- function(filelistfile="filelist.txt")
{
result <- read.delim(filelistfile,stringsAsFactors=FALSE);
print("Done reading in list of files");
return(result);
}

#' Function: readData
#'
#' Reads in input data, merging duplicate variables via geometric mean.
#' @param dataformat Either "absolute" or "fraction" or "binary" [absolute: RNA-Seq read counts, mass spec intensities; fraction: DNA methylation ratios; binary: NFR presence, variable beyond cutoff threshold].
#' @param datatype "GENE" or "COMPOUND".
#' @param datafile Full path and filename of the raw data, which must be a tab-delimited text file. \cr
#'   row 1: names of the columns
#'   column 1 (must be named: ID): variable names (genes, metabolites, etc.). \cr
#'   column 2 onwards, named by sample names, which must match those in the metadata file: measured values for each sample.
#' @return Dataframe with: \cr
#'   column ID: identifiers. \cr
#'   other columns: named as each sample, contain values. \cr
#'   each row: the values for one variable (gene, metabolite, etc.) across all samples. \cr
#'   attributes: include dataformat, datatype, datacategory, and idtype.
#' @export
readData <- function(dataformat="absolute",datatype="GENE",datacategory="RNA",idtype="KEGG",datafile="input.txt")
{
result <- read.delim(datafile,stringsAsFactors=FALSE);
attributes(result)$datatype <- datatype;
attributes(result)$dataformat <- dataformat;
attributes(result)$datacategory <- datacategory;
attributes(result)$idtype <- idtype;
result <- deduplicateByID(result);
print(paste("Done reading in raw data from ",datafile,sep=""));
return(result);
}

#' Function: removeSamples
#'
#' Removes outlier or test samples from subsequent analysis.
#' @param metadata Dataframe from readMetadata with columns Name, Group, Color, Datatype, Remove.
#' @param inputdatalist List of dataframes from readData (within each dataframe, each column represents one sample and is named to match the metadata dataframe).
#' @return List: The first element is the metadata dataframe with indicated samples (if any) removed. The second element is a list identical to inputdatalist except for the outlier samples (if any) removed.
#' @export
removeSamples <- function(metadata=data.frame(),inputdatalist=list())
{
toRemove <- metadata$Name[metadata$Remove==1]; # Get the names of samples to be removed
if(length(toRemove)==0) { return(list(metadata,inputdatalist)); } # NTD if no samples to remove
for(i in 1:length(toRemove)) { # Remove samples from inputdatalist
	for(j in 1:length(inputdatalist)) {
		if(toRemove[i] %in% colnames(inputdatalist[[j]])) {	inputdatalist[[j]][,toRemove[i]] <-NULL; break; }
	}
}
metadata <- metadata[metadata$Remove!=1,]; # Remove samples from metadata
rownames(metadata) <- 1:dim(metadata)[1];
return(list(metadata,inputdatalist));
}
