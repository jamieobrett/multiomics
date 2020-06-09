#' Function: writeOutput
#'
#' Writes the input dataframe to a csv file with the name "fileprefix_timestamp.csv".
#' @param fileprefix Name of final file (with path, otherwise writes to current directory) to precede the timestamp and extension.
#' @param writerownames Whether to write the row names of the data frame to the output.
#' @param input Input dataframe.
#' @param timestamp Timestamp to follow the file name before the extension.
#' @export
writeOutput <- function(fileprefix="output",writerownames=FALSE,input=NULL,timestamp=constants$timestamp)
{
write.csv(input,file=paste(fileprefix,timestamp,".csv",sep=''),row.names=writerownames)
}

#' Function: makeConstants
#'
#' Sets up a named list of constants: timestamp, keggGeneConvertFile.
#' @export
makeConstants <- function()
{
result <- list(
	"timestamp"=paste(format(Sys.time(),"%d%b%y_%I%p%M"),sep="")
	);
return(result);
}

#' Function: readTables
#'
#' Helper function for functions making tables if update is not desired.
#' @param tableNames Vector containing table names (example: "GENE_KEGG-KEGGPATHWAY") from which files will be read from the extdata folder and for which R tables will be named.
#' @param tableType "csv" or "txt", which applies to all table names in tableNames.
#' @param extdataFolder Name of subfolder within extdata containing the file ("GENECONVERSIONs").
#' @return Named list of dataframes.
#' @export
readTables <- function(tableNames=NULL,tableType=NULL,extdataFolder=NULL) {
result <- list();
for(i in 1:length(tableNames)) {
	if(tableType=="csv") {
		result[[tableNames[i]]] <- read.csv(system.file("extdata",extdataFolder,paste0(tableNames[i],".csv"),package="multiomics"),stringsAsFactors=FALSE);
	} else if(tableType=="txt") {
		result[[tableNames[i]]] <- read.delim(system.file("extdata",extdataFolder,paste0(tableNames[i],".txt"),package="multiomics"),stringsAsFactors=FALSE,header=TRUE);
	} else {
		stop(paste0("Could not find the table for ",tableNames[i]," of type ",tableType));
	}
}
return(result);
}

#' Function: getReadHttp
#'
#' Helper function for getting txt-format HTTP content.
#' @param txtURL Web URL from which to get content.
#' @return Vector whose elements are the lines (newline-delimited) of text from the URL.
#' @export
getReadHttp <- function(txtURL=NULL) {
getResponse <- httr::GET(txtURL);
if(httr::http_status(getResponse)$category != "Success") { stop(paste0("Could not httr::GET data for ",txtURL)); }
getContent <- suppressMessages(httr::content(getResponse,as="text"));
getLines <- strsplit(getContent,"\n")[[1]];
return(getLines);
}