#' Function: imputeZeros
#'
#' Imputes zeros to be half of the non-zero minimum of each row, or if the whole row is 0, half of the minimum of the entire matrix.
#' @param myMatrix Matrix to be used for imputation.
#' @return Matrix after imputation.
#' @export
imputeZeros <- function(myMatrix=NULL)
{ 
grandHalfMin <- min(myMatrix[myMatrix>0])/2;
for (i in 1:dim(myMatrix)[1]) {
	zeroCols <- (myMatrix[i,]==0);
	if (sum(zeroCols) < dim(myMatrix)[2]) {
		myMatrix[i,zeroCols] <- min(myMatrix[i,!zeroCols])/2;
	} else {
		myMatrix[i,zeroCols] <- grandHalfMin;
	}
}
return(myMatrix);
}

