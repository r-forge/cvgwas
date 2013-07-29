#' Genome-wide association scan based on cross-validated predictive ability
#' 
#' The function performs genome-wide association scan based on predictive
#' ability evaluated via cross validation.
#' 
#' @param formula Formula describing fixed effects to be used in analysis, e.g. 
#' y ~ a + b means that outcome (y) depends on two covariates, a and b. 
#' If no covariates used in analysis, please skip the right-hand side of the 
#' equation.
#' @param data An (optional) object of \code{\link{gwaa.data-class}} or a data frame with 
#' outcome and covariates. To read data into the require object format or to know more,
#' please refer to the \code{\link{GenABEL-package}}.
#' @param nfolds Number of folds to be established during the cross validation.
#' @param seed A positive integer to set seed for randomization so that the results
#' can be exactly repeated. If \code{NA}, no seed is used.
#' @param verbose A logic value to set whether more screen displayed output is shown
#' during the procedure.
#' @param ... Other parameters to be passed.
#' 
#' @export 
#' 
#' @details
#' The predictive ability is assessed by an $R^2$ value, which is the squared correlation
#' coefficient between the predicted phenotypes and their true values for each variant. 
#' In each fold of the cross validation, such an $R^2$ value is calculated, then the final
#' estimate of $R^2$ is given by the mean of all such values from the cross validation.
#' 
#' @return
#' A list containing \code{R2.est}, \code{R2.raw} and \code{R.raw} which provide the estimated
#' $R^2$ for each variant, the orginal $R^2$ values from cross validation, and the corresponding
#' $R$ values. 
#' 
#' @author Xia Shen
#' 
#' @references 
#' Shen X (2013). The missing heritability revealed in \emph{Arabidopsis thaliana}. 
#' \emph{Submitted}.
#' 
#' @seealso 
#' \code{\link{GenABEL-package}},
#' \code{\link{qtscore}}
#' 
#' @examples 
#' \dontrun{
#' require(GenABEL)
#' data(srdta)
#' cvres <- cvscore(qt3 ~ sex, data = srdta)
#' r2 <- cvres$R2.est
#' r2[which(is.na(r2))] <- 0
#' plot(r2, xlab = 'SNP index', ylab = 'Predictive ability')
#' }
#' 
#' @keywords cross validation, genome-wide association study, predictive ability
#' 
"cvscore" <- function(formula, data, nfolds = 5, seed = NA, verbose = TRUE, ...)
{
	options(warn = -1)
	y <- as.numeric(model.frame(formula, data = data, na.action = NULL)[,1])
	if (any(is.na(y))) {
		naidx <- which(is.na(y))
		id <- data@phdata$id[-naidx]
		n <- nrow(data@phdata) - length(naidx)
	} else {
		id <- data@phdata$id
		n <- nrow(data@phdata)
	}
	groups <- c(rep(1:(nfolds - 1), each = floor(n/nfolds)), rep(nfolds, n - (nfolds - 1)*floor(n/nfolds)))
	if (!is.na(seed)) set.seed(seed)
	groups  <- sample(groups, n)
	if (verbose) cat('Cross validation:\n')
	R.raw <- matrix(0, nfolds, ncol(data@gtdata))
	genomat <- as.double.gwaa.data(data)
	for (i in 1:nfolds) {
		gwa <- qtscore(formula, data = data, idsubset = id[groups != i])
		beta <- gwa@results$effB
		testmat <- genomat[data@phdata$id %in% id[groups == i],]
		testy <- y[data@phdata$id %in% id[groups == i]]
		R.raw[i,] <- as.numeric(cor(testy, t(t(testmat)*beta), use = "pairwise.complete.obs"))
		if (verbose) cat('fold', i, '\n')
	}
	R2.raw <- R.raw**2
	options(warn = 0)
	res <- list(R2.est = colMeans(R2.raw, na.rm = TRUE), R2.raw = R2.raw, R.raw = R.raw)
}