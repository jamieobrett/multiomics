#' Function: makeKEGGgeneTables
#'
#' Gets data from the KEGG API that maps outside IDs to KEGG IDs. KEGG identifiers are of the form org:idnumber. OTHERID identifiers have any "prefix:" removed.
#' @param updateDesired TRUE or FALSE. If TRUE, will regenerate the tables from a new web GET; each new csv file will be written to the cwd (it is up to the developer to move this to ./metabolomics/inst/extdata and overwrite prior tables if desired). If FALSE, will use existing tables in that folder.
#' @return Named list of dataframes: "KEGG-ENTREZ", "KEGG-ENTREZPROTEIN", and "KEGG-UNIPROT". Each dataframe contains two columns titled "KEGG" and "OTHERID". Mappings may not be 1:1 (multiple KEGG IDs for one ENTREZ ID are possible); this maximizes conversion ability.
#' @export
makeKEGGgeneTables <- function(updateDesired=FALSE)
{
result <- list();
URLS <- list(
	"KEGGandENTREZID"="http://rest.kegg.jp/conv/mmu/ncbi-geneid",
	"KEGGandENTREZPROTEIN"="http://rest.kegg.jp/conv/mmu/ncbi-proteinid",
	"KEGGandUNIPROT"="http://rest.kegg.jp/conv/mmu/uniprot"
);
if(updateDesired) {
	for(i in 1:length(URLS)) {
		getRowed <- getReadHttp(txtURL=URLS[[i]]);
		# Parse output into a dataframe. Input content has each entry as "ncbi-geneid:100009600", "\t", "mmu:100009600"
		resultTable <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("KEGG", "OTHERID"));
		for(j in 1:length(getRowed)) {
			thisRow <- strsplit(getRowed[j],"\t")[[1]];
			resultTableKEGG <- thisRow[2];
			resultTableOTHERID <- strsplit(thisRow[1],":")[[1]][2]; # Remove "prefix:"
			resultTable[j,] <- c(resultTableKEGG, resultTableOTHERID);
		}
		# Add to the final list of dataframes
		result[[names(URLS[i])]] <- resultTable;
		# Update the written tables
		write.csv(resultTable,file=paste(names(URLS[i]),".csv",sep=""),row.names=FALSE);
	}
} else {
	result <- readTables(tableNames=names(URLS),tableType="csv",extdataFolder="GENECONVERSIONs");
}
return(result);
}

#' Function: makeKEGGcompoundTable
#'
#' Gets data from the KEGG API that maps KEGG IDs to compound names. Converts this into a dataframe for use; note that multiple compound names are semicolon-delimited.
#' @param updateDesired TRUE or FALSE. If TRUE, will regenerate the table from a new web GET; each new tsv file will be written to the cwd (it is up to the developer to move this to ./metabolomics/inst/extdata and overwrite the prior table if desired). If FALSE, will use existing table in that folder.
#' @return Named list of single dataframe "KEGGandCOMPOUNDSandPUBCHEM" with the following columns: KEGG, OTHERID, and OTHERTYPE. OTHERTYPE may be: COMPOUND (semicolon list of multiple), PUBCHEM_SID.
#' @export
makeKEGGcompoundTable <- function(updateDesired=FALSE)
{
TABLENAME <- "KEGGandCOMPOUNDSandPUBCHEM";
if(updateDesired) {
	# Get the list of KEGG IDs mapping to compounds
	KEGGURL <- "http://rest.kegg.jp/list/compound";
	getRowed <- getReadHttp(txtURL=KEGGURL);
	# Parse the output into a dataframe with columns KEGG, OTHERID, and OTHERTYPE (OTHERTYPE is COMPOUND).
	compoundTable <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("KEGG", "OTHERID", "OTHERTYPE"));
	for(j in 1:length(getRowed)) {
		thisRow <- strsplit(getRowed[j],"\t")[[1]];
		compoundTable[j,] <- list(strsplit(thisRow[1],":")[[1]][2], thisRow[2], "COMPOUND"); # Remove "cmp:" from the compound.
	}
	# Get the list of KEGG IDs mapping to PubChem SIDs
	KEGGURL <- "http://rest.kegg.jp/conv/pubchem/compound";
	getRowed <- getReadHttp(txtURL=KEGGURL);
	# Parse the output into a dataframe with columns KEGG, OTHERID, and OTHERTYPE (OTHERTYPE is PUBCHEM_SID).
	pubchemTable <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("KEGG", "OTHERID", "OTHERTYPE"));
	for(j in 1:length(getRowed)) {
		thisRow <- strsplit(getRowed[j],"\t")[[1]];
		pubchemTable[j,] <- list(strsplit(thisRow[1],":")[[1]][2], strsplit(thisRow[2],":")[[1]][2], "PUBCHEM_SID"); # Remove "cpd:" from the KEGG id; remove "pubchem:" from the PubChem id.
	}
	resultTable <- rbind(compoundTable,pubchemTable);
	write.table(resultTable,file=paste0(TABLENAME,".txt"),row.names=FALSE,quote=FALSE,sep="\t");
	result <- list(resultTable); names(result) <- TABLENAME;
} else {
	result <- readTables(tableNames=c(TABLENAME),tableType="txt",extdataFolder="METABCONVERSIONs");
}
return(result);
}

#' Function: makeBIOCYCcompoundTable
#'
#' Gets data from the BIOCYC Web Service that maps other IDs to BioCyc IDs.
#' @param updateDesired TRUE or FALSE. If TRUE, will regenerate the table from a new web GET; each new tsv file will be written to the cwd (it is up to the developer to move this to ./metabolomics/inst/extdata and overwrite the prior table if desired). If FALSE, will use existing table in that folder.
#' @return Named list of single dataframe "KEGGandBIOCYC" with the following columns: KEGG, OTHERID, and OTHERTYPE. OTHERTYPE is: BIOCYC. Mappings may not be 1:1 (multiple KEGG IDs for one BIOCYC ID are possible); this maximizes conversion ability.
#' @export
makeBIOCYCcompoundTable <- function(updateDesired=FALSE)
{
TABLENAME <- "KEGGandBIOCYC";
if(updateDesired) {
	# KEGG
	keggTable <- makeKEGGcompoundTable(updateDesired=FALSE)[[1]];
	biocyctable <- getBIOCYC(queries=unique(keggTable$KEGG),inputID="KEGG");
	resultTable <- setNames(data.frame(matrix(ncol = 3, nrow = dim(biocyctable)[1])), c("KEGG", "OTHERID", "OTHERTYPE"));
	resultTable$KEGG <- biocyctable[,1];
	resultTable$OTHERID <- biocyctable[,2];
	resultTable$OTHERTYPE <- "BIOCYC";
	idfound <- !is.na(resultTable$OTHERID);
	resultTable <- resultTable[idfound,];
	write.table(resultTable,file=paste0(TABLENAME,".txt"),row.names=FALSE,quote=FALSE,sep="\t");
	result <- list(resultTable); names(result) <- TABLENAME;
} else {
	result <- readTables(tableNames=c(TABLENAME),tableType="txt",extdataFolder="METABCONVERSIONs");
}
return(result);
}

#' Function: getBIOCYC
#'
#' Helper for makeBIOCYCcompoundTable. Gets data from the BIOCYC Web Service (https://biocyc.org/web-services.shtml) that maps other IDs to BioCyc IDs.
#' @param queries Character vector of identifiers to match to BioCyc identifiers.
#' @param maxQueries Maximimum number of identifiers to put in the GET URL at a time.
#' @param pauseTime Seconds to pause between batch queries.
#' @param inputID "KEGG" or "PUBCHEM"
#' @return Dataframe with columns: inputID (KEGG, with later function possibly PUBCHEM), BIOCYC
#' @export
getBIOCYC <- function(queries=c(),maxQueries=10,pauseTime=1,inputID="KEGG")
{
BioCycResult <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c(inputID, "BIOCYC"));
for(i in seq(from=1, to=length(queries), by=maxQueries)) {
	maxJ <- min(i+maxQueries-1,length(queries));
	BioCycURL <- "https://websvc.biocyc.org/META/foreignid?ids=";
	BioCycURL <- paste(BioCycURL,inputID,":",queries[i],sep=""); #First ID without preceding comma.
	for(j in (i+1):maxJ) { BioCycURL <- paste(BioCycURL,",",inputID,":",queries[j],sep=""); }
	BioCycData <- getReadHttp(txtURL=BioCycURL);
	# Each element of this vector is a compound: inputid, \t, 0 or 1, \t, {outputid only if mapped}
	for(k in 1:length(BioCycData)) {
		cmpdSplit <- strsplit(BioCycData[k],"\t")[[1]]; # inputid, 0 or 1 {, outputid only if mapped}
		cmpdSplit[1] <- strsplit(cmpdSplit[1],":")[[1]][2]; # Remove "inputID:" from the inputid.
		if(cmpdSplit[2] == 0) {	cmpdSplit[3] <- NA;	}
		BioCycResult[dim(BioCycResult)[1]+1,] <- c(cmpdSplit[1],cmpdSplit[3]);
	}
	Sys.sleep(pauseTime);
	print(paste("Done with ",maxJ," BioCyc queries",sep=""));
}
return(BioCycResult);
}

#' Function: makeMBROLEcompoundTable
#'
#' Gets data from MBROLE 2.0 (http://rest.kegg.jp/list/compound) which is accessed manually using all KEGG IDs from http://rest.kegg.jp/list/compound.
#' @param updateDesired TRUE or FALSE. If TRUE, will regenerate the table from a new web GET; each new tsv file will be written to the cwd (it is up to the developer to move this to ./metabolomics/inst/extdata and overwrite the prior table if desired). If FALSE, will use existing table in that folder.
#' @return Named list of single dataframe "MBROLE_TABLE" with the following columns: KEGG, OTHERID, and OTHERTYPE. OTHERTYPE is: HMDB, LIPIDMAPS, CHEBI, CAS, PUBCHEM_CID. Mappings may not be 1:1 (multiple KEGG IDs for one BIOCYC ID are possible); this maximizes conversion ability.
#' @export
makeMBROLEcompoundTable <- function(updateDesired=FALSE) {
TABLENAME <- "MBROLE_TABLE";
if(updateDesired) {
	mbrolefile <- system.file("extdata","METABCONVERIONs","mbrole_conversion.tsv",package="multiomics");
	resultTable <- read.delim(mbrolefile,stringsAsFactors=FALSE,header=TRUE); # Dataframe with headers: Input ID, IDsource, Output ID, source of the output ID
	colnames(resultTable) <- c("KEGG","IGNORE","OTHERID","OTHERTYPE");
	resultTable$IGNORE <- NULL;
	resultTable$OTHERTYPE <- gsub("CAS registry number", "CAS", resultTable$OTHERTYPE);
	resultTable$OTHERTYPE <- gsub("ChEBI", "CHEBI", resultTable$OTHERTYPE);
	resultTable$OTHERTYPE <- gsub("LIPID MAPS", "LIPIDMAPS", resultTable$OTHERTYPE);
	#HMDB is as is
	resultTable$OTHERTYPE <- gsub("PubChem Compound", "PUBCHEM_CID", resultTable$OTHERTYPE);
	write.table(resultTable,file=paste0(TABLENAME,".txt"),row.names=FALSE,quote=FALSE,sep="\t");
	result <- list(resultTable); names(result) <- TABLENAME;
} else {
	result <- readTables(tableNames=c(TABLENAME),tableType="txt",extdataFolder="METABCONVERSIONs");
}
return(result);
}

#' Function: makeMETABOANALYSTcompoundTable
#'
#' Gets data from Metaboanalyst program (MetaboAnalyst-4.93.war resources libs compound_db.rds) from https://www.metaboanalyst.ca/docs/Resources.xhtml
#' @param updateDesired TRUE or FALSE. If TRUE, will regenerate the table from a new web GET; each new tsv file will be written to the cwd (it is up to the developer to move this to ./metabolomics/inst/extdata and overwrite the prior table if desired). If FALSE, will use existing table in that folder.
#' @return Named list of single dataframe "METABOANALYST_TABLE" with the following columns: KEGG, OTHERID, and OTHERTYPE. OTHERTYPE is: HMDB, NAME, PUBCHEM_CID, CHEBI, METLIN, SMILES.
#' @export
makeMETABOANALYSTcompoundTable <- function(updateDesired=FALSE)
{
TABLENAME <- "METABOANALYST_TABLE";
if(updateDesired) {
	metaboanalystfile <- system.file("extdata","METABCONVERSIONs","compound_db.rds",package="multiomics");
	compounds <- readRDS(file=metaboanalystfile); # Dataframe with headers: hmdb_id, name, pubchem_id, chebi_id, kegg_id, metlin_id, lipid, smiles
	colnames(compounds) <- c("HMDB","NAME","PUBCHEM_CID","CHEBI","KEGG","METLIN","LIPID","SMILES");
	compounds$LIPID <- NULL; # Not useful column (just contains 0 or 1)
	keep <- startsWith(compounds$KEGG, "C"); # Omit entries with glycans, blanks, and NAs in KEGG IDs.
	compounds <- compounds[keep,];
	resultTable <- setNames(data.frame(matrix(ncol = 3, nrow = 0)), c("KEGG", "OTHERID", "OTHERTYPE"));
	for(i in 1:dim(compounds)[1]) {
		for(j in 1:dim(compounds)[2]) {
			if((colnames(compounds)[j]!="KEGG") && (!is.na(compounds[i,j]))) {
				resultTable[dim(resultTable)[1]+1,] <- list(compounds[i,"KEGG"],compounds[i,j],colnames(compounds)[j]);
			}
		}
	}
	write.table(resultTable,file="METABOANALYST_TABLE.txt",row.names=FALSE,quote=FALSE,sep="\t"); # Not using csv since some compounds names contain commas
	result <- list(resultTable); names(result) <- TABLENAME;
} else {
	result <- readTables(tableNames=c(TABLENAME),tableType="txt",extdataFolder="METABCONVERSIONs");
}
return(result);
}

#' Function: makeMetabolitesConversionTable
#'
#' Produces a final dataframe for metabolite identifier conversions. Uses the KEGG REST API, BioCyc Web Services, MBROLE, and MetaboAnalyst data.
#' @param updateDesired TRUE or FALSE. If TRUE, will regenerate the table from a new web GET; each new tab-delimited txt file will be written to the cwd (it is up to the developer to move this to ./metabolomics/inst/extdata and overwrite the prior tables if desired). If FALSE, will use existing table in that folder.
#' @return Named list of single dataframe "METABCONVERSIONTABLE" with the following columns: KEGGID, OTHERID, OTHERTYPE. The OTHERTYPE can be BIOCYC, LIPIDMAPS, HMDB, PUBCHEM_SID, PUBCHEM_CID, CHEBI, METLIN, NAME, COMPOUND (semicolon list of multiple).
#' @export
makeMetabolitesConversionTable <- function(updateDesired=FALSE)
{
TABLENAME <- "METABCONVERSIONTABLE";
if(updateDesired) {
	# KEGG REST API
	kegg <- makeKEGGcompoundTable(updateDesired=FALSE)[[1]];
	# BioCyc Web Services
	biocyc <- makeBIOCYCcompoundTable(updateDesired=FALSE)[[1]];
	# MetaboAnalyst
	metaboanalyst <- makeMETABOANALYSTcompoundTable(updateDesired=FALSE)[[1]];
	# MBROLE
	mbrole <- makeMBROLEcompoundTable(updateDesired=FALSE)[[1]];
	# Combine
	resultTable <- rbind(kegg,biocyc,metaboanalyst,mbrole);
	# Fix ChEBI:NNN to be just NNN
	resultTable$OTHERID <- sub("CHEBI:","",resultTable$OTHERID);
	# Deduplicate entries
	resultTable <- resultTable[!duplicated(resultTable),];
	write.table(resultTable,file=paste0(TABLENAME,".txt"),row.names=FALSE,quote=FALSE,sep="\t");
	result <- list(resultTable); names(result) <- TABLENAME;
} else {
	result <- readTables(tableNames=c(TABLENAME),tableType="txt",extdataFolder="METABCONVERSIONs");
}
return(result);
}

#' Function: makePathwaysInfoTable
#'
#' Produces a final dataframe for obtaining pathway names and size from pathway IDs.
#' @param updateDesired TRUE or FALSE. If TRUE, will regenerate the table from a new web GET; each new tab-delimited txt file will be written to the cwd (it is up to the developer to move this to ./metabolomics/inst/extdata and overwrite the prior tables if desired). If FALSE, will use existing table in that folder.
#' @return Named list of single dataframe "PATHWAYSINFOTABLE" with the following columns: PATHWAY, NAME, NGENES, NCOMPOUNDS, TYPE. The TYPE can be KEGGPATHWAY, KEGGMODULE, REACTOME, WIKIPATHWAYS, BIOCYC.
#' @export
makePathwaysInfoTable <- function(updateDesired=FALSE)
{
TABLENAME <- "PATHWAYSINFOTABLE";
if(updateDesired) {
	resultTable <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("PATHWAY", "NAME", "NGENES", "NCOMPOUNDS", "TYPE"));
	# KEGG Pathways
	toadd <- makeKEGGpathwaynameTables(updateDesired=FALSE)[["KEGGPATHWAYPATHWAYNAMES"]];
	pathwaytableGenes <- makeKEGGpathwayTables(updateDesired=FALSE)[["GENE_KEGG-KEGGPATHWAY"]];
	sizes <- table(pathwaytableGenes$PATHWAY);
	toadd$NGENES <- sizes[toadd$PATHWAY];
	pathwaytableCompounds <- makeKEGGpathwayTables(updateDesired=FALSE)[["METABOLITE_KEGG-KEGGPATHWAY"]];
	sizes <- table(pathwaytableCompounds$PATHWAY);
	toadd$NCOMPOUNDS <- sizes[toadd$PATHWAY];
	toadd$TYPE <- "KEGGPATHWAY";
	resultTable <- rbind(resultTable,toadd);
	# KEGG Modules
	toadd <- makeKEGGpathwaynameTables(updateDesired=FALSE)[["KEGGMODULEPATHWAYNAMES"]];
	pathwaytableGenes <- makeKEGGpathwayTables(updateDesired=FALSE)[["GENE_KEGG-KEGGMODULE"]];
	sizes <- table(pathwaytableGenes$PATHWAY);
	toadd$NGENES <- sizes[toadd$PATHWAY];
	pathwaytableCompounds <- makeKEGGpathwayTables(updateDesired=FALSE)[["METABOLITE_KEGG-KEGGMODULE"]];
	sizes <- table(pathwaytableCompounds$PATHWAY);
	toadd$NCOMPOUNDS <- sizes[toadd$PATHWAY];
	toadd$TYPE <- "KEGGMODULE";
	resultTable <- rbind(resultTable,toadd);
	# BioCyc
	toadd <- makeBIOCYCpathwaynameTable(updateDesired=FALSE)[["BIOCYCPATHWAYNAMES"]];
	pathwaytableGenes <- makeBIOCYCpathwayTables(updateDesired=FALSE)[["GENE_MGI-BIOCYCPATHWAY"]];
	sizes <- table(pathwaytableGenes$PATHWAY);
	toadd$NGENES <- sizes[toadd$PATHWAY];
	pathwaytableCompounds <- makeBIOCYCpathwayTables(updateDesired=FALSE)[["METABOLITE_BIOCYC-BIOCYCPATHWAY"]];
	sizes <- table(pathwaytableCompounds$PATHWAY);
	toadd$NCOMPOUNDS <- sizes[toadd$PATHWAY];
	toadd$TYPE <- "BIOCYC";
	resultTable <- rbind(resultTable,toadd);
	# Reactome
	toadd <- makeREACTOMEpathwaynameTable(updateDesired=FALSE)[["REACTOMEPATHWAYNAMES"]];
	pathwaytableGenes <- makeREACTOMEpathwayTables(updateDesired=FALSE)[["GENE_ENTREZID-REACTOMEPATHWAY"]];
	sizes <- table(pathwaytableGenes$PATHWAY);
	toadd$NGENES <- sizes[toadd$PATHWAY];
	pathwaytableCompounds <- makeREACTOMEpathwayTables(updateDesired=FALSE)[["METABOLITE_CHEBI-REACTOMEPATHWAY"]];
	sizes <- table(pathwaytableCompounds$PATHWAY);
	toadd$NCOMPOUNDS <- sizes[toadd$PATHWAY];
	toadd$TYPE <- "REACTOME";
	resultTable <- rbind(resultTable,toadd);
	# Wikipathways
	toadd <- makeWIKIPATHWAYSpathwaynameTable(updateDesired=FALSE)[["WIKIPATHWAYSPATHWAYNAMES"]];
	pathwaytableGenes <- makeWIKIPATHWAYSpathwayTables(updateDesired=FALSE)[["GENE_ENTREZID-WIKIPATHWAYS"]];
	sizes <- table(pathwaytableGenes$PATHWAY);
	toadd$NGENES <- sizes[toadd$PATHWAY];
	pathwaytableCompounds <- makeWIKIPATHWAYSpathwayTables(updateDesired=FALSE)[["METABOLITE_KEGG-WIKIPATHWAYS"]];
	sizes <- table(pathwaytableCompounds$PATHWAY);
	toadd$NCOMPOUNDS <- sizes[toadd$PATHWAY];
	toadd$TYPE <- "WIKIPATHWAYS";
	resultTable <- rbind(resultTable,toadd);
	# Deduplicate entries if any are duplicated
	resultTable <- resultTable[!duplicated(resultTable),];
	# Replace NA values with 0
	resultTable[is.na(resultTable)] = 0;
	# Remove any rows where NGENES and NCOMPOUNDS are both 0
	resultTable <- resultTable[(resultTable$NGENES + resultTable$NCOMPOUNDS) > 0,];
	# Write output
	write.table(resultTable,file=paste0(TABLENAME,".txt"),row.names=FALSE,quote=FALSE,sep="\t");
	result <- list(resultTable); names(result) <- TABLENAME;
} else {
	result <- readTables(tableNames=c(TABLENAME),tableType="txt",extdataFolder="PATHWAYNAMEs");
}
return(result);
}