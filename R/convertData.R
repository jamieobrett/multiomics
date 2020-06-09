#' Function: IDtoRowName
#'
#' Names the rows as "ID" as removes this column.
#' @param input Dataframe with an "ID" column names. The other columns contain sample values. If "ID" is not a vector of unique elements, deduplicates first.
#' @return Dataframe identical to the input except: rows have been deduplicated, "ID" column is now row names.
#' @export
IDtoRowName <- function(input=NULL)
{
result <- deduplicateByID(input=input);
rownames(result) <- result[,1]; result[,1] <- NULL;
return(result);
}

#' Function: deduplicateByID
#'
#' Replaces duplicate rows (based on "ID" column value) with the single row in which the sum of the two highest expression levels is the highest.
#' @param input Dataframe with an "ID" column of (not necessarily unique) item names. The other columns contain sample values.
#' @param quieter If TRUE, does not print information about duplicates.
#' @return Dataframe identical to the input except: duplicate rows have been eliminated.
#' @export
deduplicateByID <- function(input=NULL,quieter=FALSE)
{
result <- data.frame(matrix(nrow=0,ncol=ncol(input)));
uniqueIDs <- unique(input$ID);
for(id in uniqueIDs) {
	allIDRows <- input[input$ID==id,];
	toptwo <- apply(allIDRows[,2:ncol(allIDRows)],1,function(x) {sort(x,TRUE)[1]+sort(x,TRUE)[2]} );
	rowToKeep <- names(sort(toptwo,TRUE))[1];
	result[nrow(result)+1,] <- input[rowToKeep,];
}
colnames(result) <- colnames(input); # Rename columns
attributes(result)$datatype <- attributes(input)$datatype;
attributes(result)$dataformat <- attributes(input)$dataformat;
attributes(result)$datacategory <- attributes(input)$datacategory;
attributes(result)$idtype <- attributes(input)$idtype;
if(!quieter) {
	print(paste(nrow(input) - nrow(result)," rows removed during deduplication from table of datatype ",attributes(result)$datatype,sep=""));
}
return(result);
}

#' Function: convertGeneIDs
#'
#' Converts identifiers for genes, transcripts, and proteins/enzymes. After conversion, removes items that do not convert and merges duplicates.
#' @param tableWithIDCol Dataframe to have identifiers converted. Column "ID" has the identifiers. The other columns contain each sample's values.
#' @param inputID ENTREZID, SYMBOL, ACCNUM (GenBank), ENSEMBL, ENSEMBLPROT, ENSEMBLTRANS (Ensembl Transcript), ENZYME (Enzyme Commission), ALIAS (all gene symbols, including unofficial), UNIPROT, KEGG, MGI, ENZYME, REFSEQ.
#' @param outputID Same list as inputID.
#' @param quieter If TRUE, does not print information about conversion or duplicate status.
#' @return Dataframe identical to the input except: rows whose input ID did not map to an output ID are removed, multiple rows that mapped to the same output identifier have been merged.
#' @export
convertGeneIDs <- function(tableWithIDCol=NULL,inputID="SYMBOL",outputID="ENTREZID",quieter=FALSE)
{
if(inputID=="MGI") { # org.Mm.eg.db has MGI genes in the format like MGI:MGI:99607 rather than just MGI:99607
	tableWithIDCol$ID <- gsub("(MGI:)+","MGI:MGI:",tableWithIDCol$ID);
}
if(outputID=="KEGG") { # org.Mm.eg.db doesn't have KEGG identifiers, so first convert to Entrez ID and then use KEGG's conversion table.
	if(inputID!="ENTREZID") { tableWithIDCol <- convertGeneIDs(tableWithIDCol=tableWithIDCol,inputID=inputID,outputID="ENTREZID"); }
	result <- convertKEGG(tableWithIDCol=tableWithIDCol,inputID="ENTREZID",outputID="KEGG");
} else if(inputID=="KEGG") { # org.Mm.eg.db doesn't have KEGG identifiers, so first convert to Entrez ID intermediate and then use org.Mm.eg.db.
	result <- convertKEGG(tableWithIDCol=tableWithIDCol,inputID="KEGG",outputID="ENTREZID");
	if(outputID!="ENTREZID") { result <- convertGeneIDs(tableWithIDCol=result,inputID="ENTREZID",outputID=outputID); }
} else { # No KEGG IDs involved
	if(inputID==outputID) { # AnnotationDbi does not allow this
		idConversion <- tableWithIDCol$ID;
	} else {
		if(any(match(as.character(tableWithIDCol$ID), AnnotationDbi::keys(org.Mm.eg.db::org.Mm.eg.db,keytype=inputID)),na.rm=TRUE)) { # Avoid errors if no IDs can be converted
			idConversion <- suppressMessages(AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db,keys=as.character(tableWithIDCol$ID),keytype=inputID,column=outputID));
		} else { # No IDs can be converted
			idConversion <- NULL;
		}
	}
	result <- mergeConversionVectorTable(idConversion=idConversion,tableWithIDCol=tableWithIDCol,quieter=quieter);
}
if(outputID=="MGI") { # org.Mm.eg.db has MGI genes in the format like MGI:MGI:99607 rather than just MGI:99607
	result$ID <- gsub("MGI:MGI:","MGI:",result$ID);
}
attributes(result)$idtype <- outputID;
if(dim(result)[1]==0) { # If no IDs convert at all, skip the deduplication step
	return(result);
}
result <- deduplicateByID(input=result,quieter=quieter);
return(result);
}

#' Function: simpleGeneIDConversion
#'
#' Takes as input a vector of gene identifiers and returns a dataframe containing the conversion, no other processing done.
#' @param inputvector Vector of identifiers, all of the same type.
#' @param inputID ENTREZID, SYMBOL, ACCNUM (GenBank), ENSEMBL, ENSEMBLPROT, ENSEMBLTRANS (Ensembl Transcript), ENZYME (Enzyme Commission), ALIAS (all gene symbols, including unofficial), UNIPROT, KEGG, MGI, ENZYME, REFSEQ.
#' @param outputID Same list as inputID.
#' @return Dataframe with columns INPUT (containing the input identifiers) and OUTPUT (containing the converted identifiers). If an input ID did not map, "UNMAPPED" is the output value. Duplicates are retained.
#' @export
simpleGeneIDConversion <- function(inputvector=NULL,inputID=NULL,outputID=NULL)
{
inputvector[is.na(inputvector)] <- "UNMAPPED"; # change any NA to "UNMAPPED" to avoid errors by AnnotationDbi returning NULL
if(inputID=="MGI") { # org.Mm.eg.db has MGI genes in the format like MGI:MGI:99607 rather than just MGI:99607
	inputvector <- gsub("(MGI:)+","MGI:MGI:",inputvector);
}
if(outputID=="KEGG") { # org.Mm.eg.db doesn't have KEGG identifiers, so first convert to Entrez ID and then use KEGG's conversion table.
	if(inputID!="ENTREZID") { intermediatevector <- simpleGeneIDConversion(inputvector=inputvector,inputID=inputID,outputID="ENTREZID")[,"OUTPUT"]; }
	else { intermediatevector <- inputvector; }
	keggTables <- makeKEGGgeneTables(updateDesired=FALSE);
	m <- match(intermediatevector,table=keggTables$KEGGandENTREZID[,"OTHERID"]);
	outputvector <- keggTables$KEGGandENTREZID[m,outputID];
} else if(inputID=="KEGG") { # org.Mm.eg.db doesn't have KEGG identifiers, so first convert to Entrez ID intermediate and then use org.Mm.eg.db.
	keggTables <- makeKEGGgeneTables(updateDesired=FALSE);
	m <- match(inputvector,table=keggTables$KEGGandENTREZID[,inputID]);
	intermediatevector <- keggTables$KEGGandENTREZID[m,"OTHERID"];
	if(outputID!="ENTREZID") { outputvector <- simpleGeneIDConversion(inputvector=intermediatevector,inputID="ENTREZID",outputID=outputID)[,"OUTPUT"]; }
	else { outputvector <- intermediatevector; }
} else { # No KEGG IDs involved
	if(inputID==outputID) { # AnnotationDbi does not allow this
		outputvector <- inputvector;
	} else {
		if(any(match(as.character(inputvector), AnnotationDbi::keys(org.Mm.eg.db::org.Mm.eg.db,keytype=inputID)),na.rm=TRUE)) { # Avoid errors if no IDs can be converted
			outputvector <- suppressMessages(AnnotationDbi::mapIds(org.Mm.eg.db::org.Mm.eg.db,keys=as.character(inputvector),keytype=inputID,column=outputID));
			outputvector <- as.character(outputvector);
		} else { # No IDs can be converted
			outputvector <- "UNMAPPED";
		}
	}
}
if(outputID=="MGI") { # org.Mm.eg.db has MGI genes in the format like MGI:MGI:99607 rather than just MGI:99607
	outputvector <- gsub("MGI:MGI:","MGI:",outputvector);
}
outputvector <- as.character(outputvector);
outputvector[is.na(outputvector)] <- "UNMAPPED"; # change any NA to "UNMAPPED" to avoid errors by AnnotationDbi returning NULL
result <- data.frame("INPUT"=inputvector,"OUTPUT"=outputvector);
return(result);
}

#' Function: convertKEGG
#'
#' Helper function for convertGeneIDs to interconvert between KEGG identifier and ENTREZID.
#' @param tableWithIDCol Dataframe to have identifiers converted. Row names are the identifiers. The columns contain each sample's values.
#' @param inputID KEGG or ENTREZID.
#' @param outputID ENTREZID or KEGG.
#' @return Dataframe identical to the input except: rows whose input ID did not map to an output ID are removed. Duplicates are retained.
#' @export
convertKEGG <- function(tableWithIDCol=NULL,inputID="KEGG",outputID="ENTREZID")
{
keggTables <- makeKEGGgeneTables(updateDesired=FALSE);
if(inputID=="ENTREZID") { inputID <- "OTHERID"; } # To match the columns headers of the KEGG tables.
if(outputID=="ENTREZID") { outputID <- "OTHERID"; } # To match the columns headers of the KEGG tables.
if(inputID=="KEGG" && outputID=="KEGG") {
	idConversion <- tableWithIDCol$ID;
} else {
	m <- match(tableWithIDCol$ID,table=keggTables$KEGGandENTREZID[,inputID]);
	idConversion <- keggTables$KEGGandENTREZID[m,outputID];
}
names(idConversion) <- tableWithIDCol$ID; # Formatted to match the mapIds character vector from org.Mm.eg.db.
result <- mergeConversionVectorTable(idConversion=idConversion,tableWithIDCol=tableWithIDCol);
return(result);
}

#' Function: mergeConversionVectorTable
#'
#' Helper function for convertGeneIDs and convertMetaboliteIDs, to process results of the AnnotationDbi "select" function.
#' @param idConversion Named character vector with names = input IDs, elements = converted IDs, order preserved, no duplicates (first match selected). If no match was found, the element is NA.
#' @param tableWithIDCol Dataframe to have identifiers converted. Row names are the identifiers. The columns contain each sample's values.
#' @param quieter If TRUE, does not print information about ability to perform ID conversion.
#' @return Dataframe identical to the input except: column "ID" is now the converted identifiers, any unmapped rows are removed. There may be multiple rows with the same ID (no deduplication in this function).
#' @export
mergeConversionVectorTable <- function(idConversion=c(),tableWithIDCol=NULL,quieter=FALSE)
{
idfound <- !is.na(idConversion);
result <- tableWithIDCol[idfound,]; # Remove unmapped rows.
result$ID <- idConversion[idfound];
if(!quieter) {
	print(paste(sum(is.na(idConversion))," rows had no ID conversion, from table of datatype ",attributes(tableWithIDCol)$datatype,sep=""));
}
return(result);
}

#' Function: convertMetaboliteIDs
#'
#' Converts identifiers for metabolites. After conversion, removes items that do not convert and merges duplicates.
#' @param tableWithIDCol Dataframe to have identifiers converted. Column "ID" has the identifiers. The other columns contain each sample's values.
#' @param inputID KEGG, BIOCYC, LIPIDMAPS, HMDB, PUBCHEM_SID, PUBCHEM_CID, CHEBI, METLIN, NAME.
#' @param outputID Same list of options as for inputID.
#' @param quieter If TRUE, does not print information about conversion or duplicate status.
#' @return Dataframe identical to the input except: rows whose input ID did not map to an output ID are removed, multiple rows that mapped to the same output identifier have been merged.
#' @export
convertMetaboliteIDs <- function(tableWithIDCol=NULL,inputID="KEGG",outputID="BIOCYC",quieter=FALSE)
{
metabConversionTable <- makeMetabolitesConversionTable(updateDesired=FALSE)[[1]]; # Columns: KEGG, OTHERID, OTHERTYPE
if(inputID!="KEGG") {
	# Convert to a KEGG ID as an intermediate
	inputIDsubset <- metabConversionTable[metabConversionTable$OTHERTYPE==inputID,];
	m <- match(tableWithIDCol$ID,table=inputIDsubset$OTHERID);
	keggIDs <- inputIDsubset[m,"KEGG"];
} else {
	keggIDs <- tableWithIDCol$ID;
}
# Convert to the Output ID
if(outputID!="KEGG") {
	outputIDsubset <- metabConversionTable[metabConversionTable$OTHERTYPE==outputID,];
	m <- match(keggIDs,table=outputIDsubset$KEGG);
	idConversion <- outputIDsubset[m,"OTHERID"];
} else {
	idConversion <- keggIDs;
}
# Prepare output
names(idConversion) <- tableWithIDCol$ID; # Formatted to match the mapIds character vector from org.Mm.eg.db.
result <- mergeConversionVectorTable(idConversion=idConversion,tableWithIDCol=tableWithIDCol,quieter=quieter);
attributes(result)$idtype <- outputID;
if(dim(result)[1]==0) { # If no IDs convert at all, skip the deduplication step
	return(result);
}
result <- deduplicateByID(input=result,quieter=quieter);
return(result);
}

#' Function: simpleMetabIDConversion
#'
#' Takes as input a vector of metabolite identifiers and returns a dataframe containing the conversion, no other processing done.
#' @param inputvector Vector of identifiers, all of the same type.
#' @param inputID KEGG, BIOCYC, LIPIDMAPS, HMDB, PUBCHEM_SID, PUBCHEM_CID, CHEBI, METLIN, NAME.
#' @param outputID Same list as inputID.
#' @return Dataframe with columns INPUT (containing the input identifiers) and OUTPUT (containing the converted identifiers). If an input ID did not map, "UNMAPPED" is the output value. Duplicates are retained.
#' @export
simpleMetabIDConversion <- function(inputvector=NULL,inputID=NULL,outputID=NULL)
{
inputvector[is.na(inputvector)] <- "UNMAPPED"; # change any NA to "UNMAPPED" to avoid errors by AnnotationDbi returning NULL
if(inputID==outputID) { #save time
	outputvector <- inputvector;
} else {
	metabConversionTable <- makeMetabolitesConversionTable(updateDesired=FALSE)[[1]]; # Columns: KEGG, OTHERID, OTHERTYPE
	if(inputID!="KEGG") {
		# Convert to a KEGG ID as an intermediate
		inputIDsubset <- metabConversionTable[metabConversionTable$OTHERTYPE==inputID,];
		m <- match(inputvector,table=inputIDsubset$OTHERID);
		keggIDs <- inputIDsubset[m,"KEGG"];
	} else {
		keggIDs <- inputvector;
	}
	# Convert to the Output ID
	if(outputID!="KEGG") {
		outputIDsubset <- metabConversionTable[metabConversionTable$OTHERTYPE==outputID,];
		m <- match(keggIDs,table=outputIDsubset$KEGG);
		idConversion <- outputIDsubset[m,"OTHERID"];
	} else {
		idConversion <- keggIDs;
	}
	# Prepare output
	outputvector <- idConversion;
	outputvector[is.na(outputvector)] <- "UNMAPPED"; # change any NA to "UNMAPPED"
}
result <- data.frame("INPUT"=inputvector,"OUTPUT"=outputvector);
return(result);
}

#' Function: convertPathwayInfo
#'
#' Converts Pathway identifiers into other pathway information.
#' @param pathwayIDs Vector with pathway identifiers.
#' @param type KEGGPATHWAY, KEGGMODULE, REACTOME, WIKIPATHWAYS, BIOCYC
#' @param outputInfo NAME, NGENES, NCOMPOUNDS
#' @return Vector with the information, in the same order as the input vector.
#' @export
convertPathwayInfo <- function(pathwayIDs=NULL,type="KEGGPATHWAY",outputInfo="NAME") {
infoConversionTable <- makePathwaysInfoTable(updateDesired=FALSE)[[1]];
typesubset <- infoConversionTable[infoConversionTable$TYPE==type,];
result <- typesubset[match(tolower(pathwayIDs),tolower(typesubset$PATHWAY)),outputInfo];
return(result);
}