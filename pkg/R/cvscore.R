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
"cvscore" <- function(formula, data, nfolds = 5, seed = NA, verbose = TRUE, memory.saving = FALSE, ...)
{
	options(warn = -1)
	y <- as.numeric(model.frame(formula, data = data, na.action = NULL)[,1])
	if (any(is.na(y))) {
		naidx <- which(is.na(y))
		id <- data@phdata$id[-naidx]
		n <- nrow(data@phdata) - length(naidx)
		m0 <- lm(formula, data = data@phdata[-naidx,])
	} else {
		id <- data@phdata$id
		n <- nrow(data@phdata)
		m0 <- lm(formula, data = data@phdata)
	}
	y <- m0$resid
	
	groups <- c(rep(1:(nfolds - 1), each = floor(n/nfolds)), rep(nfolds, n - (nfolds - 1)*floor(n/nfolds)))
	if (!is.na(seed)) set.seed(seed)
	groups  <- sample(groups, n)
	
	R.raw <- MSE <- matrix(0, nfolds, ncol(data@gtdata))
	N <- as.numeric(nrow(data@gtdata))
	p <- as.numeric(ncol(data@gtdata))
	if (!memory.saving) {
		genomat <- c()
		if (verbose) cat('Allocating genotype matrix ...\n')
		if (verbose) pb <- txtProgressBar(style = 3)
		if (N*p < 9e8) {
			genomat <- as.double.gwaa.data(data)
			if (verbose) setTxtProgressBar(pb, 1)
		} else {
			npiece <- ceiling(N*p/9e8)
			nc <- c(rep(floor(p/npiece), npiece - 1), p - (npiece - 1)*floor(p/npiece))
			cumnc <- c(0, cumsum(nc))
			for (i in 1:npiece) {
				genomat <- cbind(genomat, as.double.gwaa.data(data[,(cumnc[i] + 1):(cumnc[i] + nc[i])]))
				if (verbose) setTxtProgressBar(pb, i/npiece)
			}
		}
		if (verbose) cat('\nCross validation:\n')
		for (i in 1:nfolds) {
			gwa <- qtscore(y[groups != i], data = data, idsubset = id[groups != i])
			beta <- gwa@results$effB
			testmat <- genomat[data@phdata$id %in% id[groups == i],]
			testy <- y[data@phdata$id %in% id[groups == i]]
			#R.raw[i,] <- as.numeric(cor(testy, t(t(testmat)*beta), use = "pairwise.complete.obs"))
			R.raw[i,] <- as.numeric(cor(testy, testmat, use = "pairwise.complete.obs"))
			MSE[i,] <- colMeans((testy - t(t(testmat)*beta))**2)
			if (verbose) cat('fold', i, '\n')
		}
	} else {
		if (verbose) cat('\nCross validation:\n')
		for (i in 1:nfolds) {
			gwa <- qtscore(y[groups != i], data = data, idsubset = id[groups != i])
			beta <- gwa@results$effB
			if (N*p < 9e8) {
				genomat <- as.double.gwaa.data(data)
				testmat <- genomat[data@phdata$id %in% id[groups == i],]
				testy <- y[data@phdata$id %in% id[groups == i]]
				#R.raw[i,] <- as.numeric(cor(testy, t(t(testmat)*beta), use = "pairwise.complete.obs"))
				R.raw[i,] <- as.numeric(cor(testy, testmat, use = "pairwise.complete.obs"))
				MSE[i,] <- colMeans((testy - t(t(testmat)*beta))**2)
			} else {
				npiece <- ceiling(N*p/9e8)
				nc <- c(rep(floor(p/npiece), npiece - 1), p - (npiece - 1)*floor(p/npiece))
				cumnc <- c(0, cumsum(nc))
				testy <- y[data@phdata$id %in% id[groups == i]]
				for (j in 1:npiece) {
					testmat <- as.double.gwaa.data(data[,(cumnc[j] + 1):(cumnc[j] + nc[j])])
					testmat <- testmat[data@phdata$id %in% id[groups == i],]
					#R.raw[i,(cumnc[j] + 1):(cumnc[j] + nc[j])] <- as.numeric(cor(testy, t(t(testmat)*beta[(cumnc[j] + 1):(cumnc[j] + nc[j])]), use = "pairwise.complete.obs"))
					R.raw[i,(cumnc[j] + 1):(cumnc[j] + nc[j])] <- as.numeric(cor(testy, testmat, use = "pairwise.complete.obs"))
					MSE[i,(cumnc[j] + 1):(cumnc[j] + nc[j])] <- colMeans((testy - t(t(testmat)*beta[(cumnc[j] + 1):(cumnc[j] + nc[j])]))**2)
				}
			}
			if (verbose) cat('fold', i, '\n')
		}
	}

	R2.raw <- R.raw**2
	R2.est <- colMeans(R2.raw, na.rm = TRUE)
	RMSE <- sqrt(MSE)
	R2.se <- sqrt(colVars(R2.raw)/nfolds)
	options(warn = 0)
	res <- list(R2.est = R2.est, R2.se = R2.se, R2.pt = pt(R2.est/R2.se, nfolds - 1, lower.tail = FALSE),
			    R2.raw = R2.raw, R.raw = R.raw, MRMSE = colMeans(RMSE), RMSE = RMSE)
}

'colVars' <- function(x, na.rm = TRUE) colSums((t(t(x) - colMeans(x, na.rm = na.rm)))**2, na.rm = na.rm)/(nrow(x) - 1)