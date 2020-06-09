#' Function: aggregateRanks
#'
#' Uses the RankAggreg package to find a consensus ranked list for multiple input ranked lists. Uses the Cross-Entropy Monte Carlo algorithm and the Spearman footrule distance.
#' @param enrichedpathwayslist List of dataframes, each of the format returned by indivEnrichment.
#' @param pathwaysDb KEGGPATHWAY, KEGGMODULE, BIOCYC, WIKIPATHWAYS, REACTOME
#' @param rankTypeScore To be passed to prepareRankLists: POS, NEG, ABS, NONE.
#' @param rankTypeQ To be passed to prepareRankLists: POS, NEG, ABS, NONE.
#' @param rankTypeP To be passed to prepareRankLists: POS, NEG, ABS, NONE.
#' @param cutoffType To be passed to prepareRankLists: TOPSCORE, TOPQ, TOPP, TOPN, TOPFRAC, NONE.
#' @param cutoffValue To be passed to prepareRankLists: value or NULL.
#' @param maxPathwaysToReturn Number of pathways desired in the consensus list. Any more than 20 requires a large amount of time and CPU use. Must be no larger than the largest input list.
#' @param aggregationScoreType SCORE, QVAL, PVAL, or NONE. Used to refine the distance function if not NONE. Any non-NONE significantly increases computational load.
#' @param mySeed A number or NULL.
#' @param importances Vector in the same order as the ordered lists with a number indicating the importance of each to be used in distance metric calculation. Higher number means more important. Default (NULL) is equal importance to all.
#' @return Dataframe with columns "PATHID", "PATHNAME", "PVAL", "QVAL" along with attributes "iterations" and "objectiveFunctionScore". The p-value is the chance of observing the observed set of ranks for this pathway under the null hypothesis that each the ranks of each input list are random. The q-value is the BH-corrected FDR.
#' @export
aggregateRanks <- function(enrichedpathwayslist=NULL,pathwaysDb="KEGGPATHWAY",rankTypeScore="NONE",rankTypeQ="NONE",rankTypeP="NONE",cutoffType="NONE",cutoffValue="NONE",maxPathwaysToReturn=20,aggregationScoreType="NONE",mySeed=NULL,importances=NULL)
{
# Prepare the pathways and parameters
enrichedpathwayslist <- prepareRankLists(enrichedpathwayslist=enrichedpathwayslist,rankTypeScore=rankTypeScore,rankTypeQ=rankTypeQ,rankTypeP=rankTypeP,cutoffType=cutoffType,cutoffValue=cutoffValue);
maxPathways <- nrow(enrichedpathwayslist[[which.max(lapply(enrichedpathwayslist,nrow))]]);
pathwayMatrix <- matrix(nrow=length(enrichedpathwayslist),ncol=maxPathways);
rownames(pathwayMatrix) <- rep(NA,nrow(pathwayMatrix));
if(aggregationScoreType=="NONE") {
	weightsMatrix <- NULL;
} else {
	weightsMatrix <- matrix(nrow=nrow(pathwayMatrix),ncol=ncol(pathwayMatrix));
}
for(i in 1:length(enrichedpathwayslist)) {
	rowVec <- enrichedpathwayslist[[i]][,"PATHID"];
	placeholderVec <- rep("PLACEHOLDER",maxPathways-length(rowVec));
	if(length(placeholderVec)>0) {
		placeholderVec <- paste(placeholderVec,1:length(placeholderVec),sep="");
	}
	rowVec <- c(rowVec,placeholderVec);
	pathwayMatrix[i,] <- rowVec;
	rownames(pathwayMatrix)[i] <- attributes(enrichedpathwayslist[[i]])$datacategory;
	if(aggregationScoreType!="NONE") {
		rowVec <- enrichedpathwayslist[[i]][,aggregationScoreType];
		placeholder <- rev(enrichedpathwayslist[[i]][,aggregationScoreType])[1]; #the bottom element
		rowVec <- c(rowVec,rep(placeholder,maxPathways-length(rowVec)));
		weightsMatrix[i,] <- rowVec;
	}
}
if(is.null(importances)) { importances <- rep(1,nrow(pathwayMatrix)); }
result <- RankAggreg::RankAggreg(pathwayMatrix, k=maxPathwaysToReturn, weights=weightsMatrix, seed=mySeed, method="CE", distance="Spearman", maxIter=1000, convIn=7, importance=importances, rho=0.01, verbose=FALSE);
print(paste0("RankAggreg converged in ",result$num.iter," iterations with objective function value ",result$optimal.value));
resultNames <- convertPathwayInfo(pathwayIDs=result$top.list,type=pathwaysDb,outputInfo="NAME");
resultTable <- data.frame("PATHID"=result$top.list,"PATHNAME"=resultNames);
# Get p-values
allPathways <- c();
for (i in 1:nrow(pathwayMatrix)) {
	allPathways <- c(allPathways,enrichedpathwayslist[[i]][,"PATHID"]);
}
allPathways <- unique(allPathways);
maxPathwayCounts <- unlist(lapply(enrichedpathwayslist,nrow));
normalizedRanksDf <- setNames(data.frame(matrix(ncol=nrow(pathwayMatrix),nrow=length(allPathways))), rownames(pathwayMatrix));
rownames(normalizedRanksDf) <- allPathways;
for (i in 1:nrow(normalizedRanksDf)) {
for (j in 1:ncol(normalizedRanksDf)) {
	thisRank <- match(rownames(normalizedRanksDf)[i],enrichedpathwayslist[[j]]$PATHID);
	if(is.na(thisRank)) {
		normalizedRanksDf[i,j] <- NA;
	}
	else {
		normalizedRanksDf[i,j] <- thisRank/maxPathwayCounts[j];
	}
}}
pvalues <- getRankVectorPvals(normalizedRanksDf);
resultTable$PVAL <- pvalues[match(resultTable$PATHID,names(pvalues))];
resultTable$QVAL <- p.adjust(resultTable$PVAL,method="BH");
# Finalize and return
attributes(resultTable)$iterations <- result$num.iter;
attributes(resultTable)$objectiveFunctionScore <- result$optimal.value;
return(resultTable);
}

#' Function: getRankVectorPvals
#' Uses the RobustRankAggreg package (rhoScores function) to generate p-values. For computational efficiency and stability, returns the upper limit of the p-value (conservative) rather than the exact p-value.
#' @param normalizedRanksDf Dataframe in which each pathway is a named row, and each column is a dataset category. Each element is the normalized rank of that pathway in that dataset category. Normalized means divided by the total number of ranks: values are (0,1]. Any pathway with a missing rank should have its normalized rank set to 1.
#' @return Named vector of p-values in the same order as the rows of the input dataframe. The p-value is the chance of observing the observed set of ranks for this pathway under the null hypothesis that each the ranks of each input list are random. Note that this is not corrected for multiple hypothesis testing.
#' @export
getRankVectorPvals <- function(normalizedRanksDf=NULL)
{
pvals <- apply(normalizedRanksDf,1,RobustRankAggreg::rhoScores,exact=TRUE);
return(pvals);
}

#' Function: robustAggregateRanks
#'
#' Uses the RobustRankAggreg package to generate a consensus ranked list for multiple input ranked lists.
#' @param enrichedpathwayslist List of dataframes, each of the format returned by indivEnrichment.
#' @param pathwaysDb KEGGPATHWAY, KEGGMODULE, BIOCYC, WIKIPATHWAYS, REACTOME
#' @param rankTypeScore To be passed to prepareRankLists: POS, NEG, ABS, NONE.
#' @param rankTypeQ To be passed to prepareRankLists: POS, NEG, ABS, NONE.
#' @param rankTypeP To be passed to prepareRankLists: POS, NEG, ABS, NONE.
#' @param cutoffType  To be passed to prepareRankLists: TOPSCORE, TOPSIG, TOPN, TOPFRAC, NONE.
#' @param cutoffValue To be passed to prepareRankLists: value or NULL.
#' @return Dataframe with columns "PATHID","PATHNAME","PVAL","QVAL". The p-value is the chance of observing the observed set of ranks for this pathway under the null hypothesis that each the ranks of each input list are random. The q-value is the BH-corrected FDR.
#' @export
robustAggregateRanks <- function(enrichedpathwayslist=NULL,pathwaysDb="KEGGPATHWAY",rankTypeScore="NONE",rankTypeQ="NONE",rankTypeP="NONE",cutoffType="NONE",cutoffValue=NULL)
{
# Prepare the pathways and parameters
enrichedpathwayslist <- prepareRankLists(enrichedpathwayslist=enrichedpathwayslist,rankTypeScore=rankTypeScore,rankTypeQ=rankTypeQ,rankTypeP=rankTypeP,cutoffType=cutoffType,cutoffValue=cutoffValue);
if(cutoffType=="NONE") {
	fullParameter <- TRUE; # Missing pathways will be set to rank NA
	nParameter <- NA; # Will be autocalculated as the total unique pathways (since not cutoffs were used)
	topCutoff <- NA;
} else {
	fullParameter <- FALSE; # Missing pathways will be set to bottommost rank
	nParameter <- attributes(enrichedpathwayslist)$originalpathwaycount; # Use the total unique pathways before cutoffs
	topCutoff <- c();
	for(i in 1:length(enrichedpathwayslist)) {
		topCutoff <- c(topCutoff,attributes(enrichedpathwayslist[[i]])$currentpathwaycount/attributes(enrichedpathwayslist[[i]])$originalpathwaycount);
	}
}
# Make the glist input
glist <- list();
for(i in 1:length(enrichedpathwayslist)) {
	glist[[attributes(enrichedpathwayslist[[i]])$datacategory]] <- enrichedpathwayslist[[i]][,"PATHID"];
}
# Make the rank matrix input
rmatN <- c();
for(i in 1:length(enrichedpathwayslist)) {
	rmatN <- c(rmatN,attributes(enrichedpathwayslist[[i]])$originalpathwaycount);
}
rmat <- RobustRankAggreg::rankMatrix(glist,N=rmatN,full=fullParameter);
# Run the rank aggregation
result <- RobustRankAggreg::aggregateRanks(glist, rmat=rmat, N=nParameter, method="RRA", full=fullParameter, exact=TRUE, topCutoff=topCutoff); #exact FALSE to use distributions for speed of p-value calculation
resultNames <- convertPathwayInfo(pathwayIDs=result$Name,type=pathwaysDb,outputInfo="NAME");
resultTable <- data.frame("PATHID"=result$Name,"PATHNAME"=resultNames,"PVAL"=result$Score);
resultTable$QVAL <- p.adjust(resultTable$PVAL,method="BH");
return(resultTable);
}

#' Function: prepareRankLists
#'
#' Helper function for list aggregation. Transforms, filters, and re-ranks lists based on specified criteria. Note that for p-value and q-value log-transformation, any 0 value is first imputed to be the half-min of non-zero value of all the p-values or q-values, respectively.
#' @param enrichedpathwayslist List of dataframes for filtering, each of the format returned by indivEnrichment.
#' @param rankTypeScore POS, NEG, ABS, NONE. \cr
#'  POS means the most positively scoring pathways will be ranked at the top (ties broken by QVAL). Negatively scoring pathways' p-values and q-values will be set to 1 and then transformed by -log2. Scores will be multiplied by 1 for use in aggregation. \cr
#'  NEG means the most negatively scoring pathways will be ranked at the top (ties broken by QVAL). Positively scoring pathways' p-values and q-values will be set to 1 and then transformed by -log2. Scores will be multiplied by -1 for use in aggregation. \cr
#'  ABS means scores will have their absolute value taken and ranked based on this (ties broken by QVAL). P-values and q-values will be transformed by -log2. \cr
#'  NONE means rankTypeP or rankTypeQ will be used instead.
#' @param rankTypeQ POS, NEG, ABS, NONE. \cr
#'  POS means the most significant positively scoring pathways will be ranked at the top (ties broken by SCORE). Negatively scoring pathways' p-values and q-values will be set to 1 and then transformed by -log2. Scores remain unchanged for future reference. \cr
#'  NEG means the most significant negatively scoring pathways will be ranked at the top (ties broken by SCORE). Positively scoring pathways' p-values and q-values will be set to 1 and then transformed by -log2. Scores remain unchanged for future reference. \cr
#'  ABS means pathways will be ranked based on q-value (ties broken by absolute value of SCORE). Before ranking, p-values and q-values will be transformed by -log2. Scores remain unchanged for future reference. \cr
#'  NONE means rankTypeScore will be used instead.
#' @param rankTypeP POS, NEG, ABS, NONE. See rankTypeQ for an explanation; p-values are used instead.
#' @param cutoffType After re-ranking and score/q-value transformation, type of cutoff to implement: TOPSCORE, TOPSIG, TOPN, TOPFRAC, NONE. For example, if choosing rankTypeScore="NEG", the selected cutoff value should be positive. Or if choosing cutoffType=TOPSIG, the selected cutoff value should be of the form -log2(q-value). \cr
#'  TOPSCORE means pathways must have at least this score to be retained. Note that this will be applied to the post-transformation data (if rankTypeScore is NEG, provide a positive value in cutoffValue). \cr
#'  TOPSIG means pathways must have at most this significance value to be retained. Note that this will be applied to the post-transformation data. \cr
#'  TOPN means up to this many pathways based on the post-transformation ranking will be retained. \cr
#'  TOPFRAC means this fraction of pathways (rounded down) based on the post-transformation ranking will be retained. \cr
#'  NONE means no filtering will be performed.
#' @param cutoffValue Value to be used in the cutoff, or NULL.
#' @return List of same format as enrichedpathwayslist except dataframes are post-transformation and post-cutoff application, each dataframe in the list now has the attributes "originalpathwaycount" and "currentpathwaycount", and the entire list now has the attributes total unique "originalpathwaycount" and "currentpathwaycount".
#' @export
prepareRankLists <- function(enrichedpathwayslist=NULL,rankTypeScore="NONE",rankTypeQ="NONE",rankTypeP="NONE",cutoffType="NONE",cutoffValue=NULL)
{
# Get the total number of pathways before the cutoffs
allPathways <- c();
for(i in 1:length(enrichedpathwayslist)) {
	allPathways <- c(allPathways,enrichedpathwayslist[[i]]$PATHID);
	attributes(enrichedpathwayslist[[i]])$originalpathwaycount <- length(unique(enrichedpathwayslist[[i]]$PATHID));
}
attributes(enrichedpathwayslist)$originalpathwaycount <- length(unique(allPathways));
# Before doing log-transformations on significance values, impute any p-value or q-value of 0 to be the half-min of the smallest non-zero p-value or q-value, respectively.
for(i in 1:length(enrichedpathwayslist)) {
	zeroP=enrichedpathwayslist[[i]]$PVAL==0;
	zeroQ=enrichedpathwayslist[[i]]$QVAL==0;
	halfminP=min(enrichedpathwayslist[[i]][!zeroP,"PVAL"])/2;
	halfminQ=min(enrichedpathwayslist[[i]][!zeroP,"QVAL"])/2;
	enrichedpathwayslist[[i]][zeroP,"PVAL"] <- halfminP;
	enrichedpathwayslist[[i]][zeroP,"QVAL"] <- halfminQ;
}
# Reorder and transform
if(rankTypeScore=="POS") {
	for(i in 1:length(enrichedpathwayslist)) {
		enrichedpathwayslist[[i]][enrichedpathwayslist[[i]]$SCORE<0,c("PVAL","QVAL")] <- 1;
		enrichedpathwayslist[[i]] <- enrichedpathwayslist[[i]][order(enrichedpathwayslist[[i]]$SCORE,(-1)*enrichedpathwayslist[[i]]$QVAL,decreasing=TRUE),];
	}
} else if(rankTypeScore=="NEG") {
	for(i in 1:length(enrichedpathwayslist)) {
		enrichedpathwayslist[[i]][enrichedpathwayslist[[i]]$SCORE>0,c("PVAL","QVAL")] <- 1;
		enrichedpathwayslist[[i]]$SCORE <- enrichedpathwayslist[[i]]$SCORE*(-1);
		enrichedpathwayslist[[i]] <- enrichedpathwayslist[[i]][order(enrichedpathwayslist[[i]]$SCORE,(-1)*enrichedpathwayslist[[i]]$QVAL,decreasing=TRUE),];
	}
} else if(rankTypeScore=="ABS") {
	for(i in 1:length(enrichedpathwayslist)) {
		enrichedpathwayslist[[i]]$SCORE <- abs(enrichedpathwayslist[[i]]$SCORE);
		enrichedpathwayslist[[i]] <- enrichedpathwayslist[[i]][order(enrichedpathwayslist[[i]]$SCORE,(-1)*enrichedpathwayslist[[i]]$QVAL,decreasing=TRUE),];
	}
}
if(rankTypeQ=="POS") {
	for(i in 1:length(enrichedpathwayslist)) {
		enrichedpathwayslist[[i]][enrichedpathwayslist[[i]]$SCORE<0,c("PVAL","QVAL")] <- 1;
		enrichedpathwayslist[[i]] <- enrichedpathwayslist[[i]][order(enrichedpathwayslist[[i]]$QVAL,(-1)*enrichedpathwayslist[[i]]$SCORE,decreasing=FALSE),];
	}
} else if(rankTypeQ=="NEG") {
	for(i in 1:length(enrichedpathwayslist)) {
		enrichedpathwayslist[[i]][enrichedpathwayslist[[i]]$SCORE>0,c("PVAL","QVAL")] <- 1;
		enrichedpathwayslist[[i]] <- enrichedpathwayslist[[i]][order(enrichedpathwayslist[[i]]$QVAL,enrichedpathwayslist[[i]]$SCORE,decreasing=FALSE),];
	}
} else if(rankTypeQ=="ABS") {
	for(i in 1:length(enrichedpathwayslist)) {
		enrichedpathwayslist[[i]] <- enrichedpathwayslist[[i]][order(enrichedpathwayslist[[i]]$QVAL,(-1)*abs(enrichedpathwayslist[[i]]$SCORE),decreasing=FALSE),];
	}
}
if(rankTypeP=="POS") {
	for(i in 1:length(enrichedpathwayslist)) {
		enrichedpathwayslist[[i]][enrichedpathwayslist[[i]]$SCORE<0,c("PVAL","QVAL")] <- 1;
		enrichedpathwayslist[[i]] <- enrichedpathwayslist[[i]][order(enrichedpathwayslist[[i]]$PVAL,(-1)*enrichedpathwayslist[[i]]$SCORE,decreasing=FALSE),];
	}
} else if(rankTypeP=="NEG") {
	for(i in 1:length(enrichedpathwayslist)) {
		enrichedpathwayslist[[i]][enrichedpathwayslist[[i]]$SCORE>0,c("PVAL","QVAL")] <- 1;
		enrichedpathwayslist[[i]] <- enrichedpathwayslist[[i]][order(enrichedpathwayslist[[i]]$PVAL,enrichedpathwayslist[[i]]$SCORE,decreasing=FALSE),];
	}
} else if(rankTypeP=="ABS") {
	for(i in 1:length(enrichedpathwayslist)) {
		enrichedpathwayslist[[i]] <- enrichedpathwayslist[[i]][order(enrichedpathwayslist[[i]]$PVAL,(-1)*abs(enrichedpathwayslist[[i]]$SCORE),decreasing=FALSE),];
	}
}
# Transform all the p-values and q-values by -log2.
for(i in 1:length(enrichedpathwayslist)) {
	enrichedpathwayslist[[i]]$PVAL <- -log2(enrichedpathwayslist[[i]]$PVAL);
	enrichedpathwayslist[[i]]$PVAL <- -log2(enrichedpathwayslist[[i]]$QVAL);
}
# Cutoffs
if(cutoffType=="TOPSCORE") {
	for(i in 1:length(enrichedpathwayslist)) {
		enrichedpathwayslist[[i]] <- enrichedpathwayslist[[i]][enrichedpathwayslist[[i]]$SCORE>=cutoffValue,];
	}
} else if(cutoffType=="TOPSIG") {
	for(i in 1:length(enrichedpathwayslist)) {
		enrichedpathwayslist[[i]] <- enrichedpathwayslist[[i]][enrichedpathwayslist[[i]]$QVAL<=cutoffValue,];
	}
} else if(cutoffType=="TOPN") {
	for(i in 1:length(enrichedpathwayslist)) {
		if(cutoffValue < attributes(enrichedpathwayslist[[i]])$originalpathwaycount) {
			enrichedpathwayslist[[i]] <- enrichedpathwayslist[[i]][seq(1,cutoffValue),];
		}
	}
} else if(cutoffType=="TOPFRAC") {
	for(i in 1:length(enrichedpathwayslist)) {
		myN <- floor(cutoffValue*attributes(enrichedpathwayslist[[i]])$originalpathwaycount);
		enrichedpathwayslist[[i]] <- enrichedpathwayslist[[i]][seq(1,myN),];
	}
}
# Get the total number of pathways after the cutoffs
allPathways <- c();
for(i in 1:length(enrichedpathwayslist)) {
	allPathways <- c(allPathways,enrichedpathwayslist[[i]]$PATHID);
	attributes(enrichedpathwayslist[[i]])$currentpathwaycount <- length(unique(enrichedpathwayslist[[i]]$PATHID));
}
attributes(enrichedpathwayslist)$currentpathwaycount <- length(unique(allPathways));
# Return
return(enrichedpathwayslist);
}
