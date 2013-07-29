\name{cvscore}
\alias{cvscore}
\title{Genome-wide association scan based on cross-validated predictive ability}
\usage{
  cvscore(formula, data, nfolds = 5, seed = NA,
    verbose = TRUE, ...)
}
\arguments{
  \item{formula}{Formula describing fixed effects to be
  used in analysis, e.g. y ~ a + b means that outcome (y)
  depends on two covariates, a and b. If no covariates used
  in analysis, please skip the right-hand side of the
  equation.}

  \item{data}{An (optional) object of
  \code{\link{gwaa.data-class}} or a data frame with
  outcome and covariates. To read data into the require
  object format or to know more, please refer to the
  \code{\link{GenABEL-package}}.}

  \item{nfolds}{Number of folds to be established during
  the cross validation.}

  \item{seed}{A positive integer to set seed for
  randomization so that the results can be exactly
  repeated. If \code{NA}, no seed is used.}

  \item{verbose}{A logic value to set whether more screen
  displayed output is shown during the procedure.}

  \item{...}{Other parameters to be passed.}
}
\value{
  A list containing \code{R2.est}, \code{R2.raw} and
  \code{R.raw} which provide the estimated $R^2$ for each
  variant, the orginal $R^2$ values from cross validation,
  and the corresponding $R$ values.
}
\description{
  The function performs genome-wide association scan based
  on predictive ability evaluated via cross validation.
}
\details{
  The predictive ability is assessed by an $R^2$ value,
  which is the squared correlation coefficient between the
  predicted phenotypes and their true values for each
  variant. In each fold of the cross validation, such an
  $R^2$ value is calculated, then the final estimate of
  $R^2$ is given by the mean of all such values from the
  cross validation.
}
\examples{
\dontrun{
require(GenABEL)
data(srdta)
cvres <- cvscore(qt3 ~ sex, data = srdta)
r2 <- cvres$R2.est
r2[which(is.na(r2))] <- 0
plot(r2, xlab = 'SNP index', ylab = 'Predictive ability')
}
}
\author{
  Xia Shen
}
\references{
  Shen X (2013). The missing heritability revealed in
  \emph{Arabidopsis thaliana}. \emph{Submitted}.
}
\seealso{
  \code{\link{GenABEL-package}}, \code{\link{qtscore}}
}
\keyword{ability}
\keyword{association}
\keyword{cross}
\keyword{genome-wide}
\keyword{predictive}
\keyword{study,}
\keyword{validation,}

