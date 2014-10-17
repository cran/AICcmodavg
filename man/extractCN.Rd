\name{extractCN}
\Rdversion{1.1}
\alias{extractCN}
\alias{extractCN.default}
\alias{extractCN.unmarkedFit}
\alias{print.extractCN}
\title{
Compute Condition Number
}
\description{
This function computes the condition number for models of
\code{unmarkedFit} classes as the ratio of the largest eigenvalue of the
Hessian matrix to the smallest eigenvalue of the Hessian matrix.
}
\usage{
extractCN(mod, method = "svd", \dots)

\method{extractCN}{unmarkedFit}(mod, method = "svd", \dots)
}
\arguments{
  \item{mod}{
    a model of one the \code{unmarkedFit} classes for which a condition
    number is requested.
}
\item{method}{
  specifies the method used to extract the singular values or
  eigenvalues from the Hessian matrix using singular value decomposition
  (\code{method = "svd"}) or eigenvalue decomposition
  (\code{method = "eigen"}).
}
\item{\dots}{
  additional arguments passed to the function.
  }
}
\details{
The condition number (\eqn{\kappa}) is a measure of the transfer of
error to the solution in response to small changes in the input (Cheney
and Kincaid 2008).  In this implementation, the condition number is
computed on the Hessian matrix of models of \code{unmarkedFit} classes
from the \code{optim} results stored in the model object.  The condition
number is defined as the ratio of the largest to the smallest
non-negative singular values of a given matrix (Cline et al. 1979, Dixon
1983).  In the special case of positive semi-definite matrices, the
singular values are equal to the eigenvalues (Ruhe 1975). 

Large values of the condition number may indicate problems in estimating
parameters or their variance (ill-conditioning), possibly due to a model
having too many parameters for the given data set.  Cheney and Ward
(2008) suggest using the \eqn{\log_{10}(\kappa)} of the condition number as
a crude estimate of the number of digits of precision lost.
}
\value{
  \code{extractCN} returns a list of class \code{extractCN} with the
  following components:
  
  \item{CN}{the condition number (\eqn{\kappa}) of the model.
  }
  \item{log10}{the log base 10 of the condition number.
  }
  \item{method}{the method used to extract the singular values or eigenvalues.
  }
}
\references{

Cheney, W., Kincaid, D. (2008) \emph{Numerical mathematics and
  computing}. Sixth edition. Thomson Brooks/Cole: Belmont.

Cline, A. K., Moler, C. B., Stewart, G. W., Wilkinson, J. H. (1979) An
estimate for the condition number of a matrix. \emph{SIAM Journal on
  Numerical Analysis} \bold{16}, 368--375.

Dixon, J. D. (1983) Estimating extremal eigenvalues and condition
numbers of matrices. \emph{SIAM Journal on Numerical Analysis}
\bold{20}, 812--814. 

Ruhe, A. (1975) On the closeness of eigenvalues and singular values for
almost normal matrices. \emph{Linear Algebra and its Applications}
\bold{11}, 87--94. 
}

\author{
Marc J. Mazerolle
}

\seealso{
\code{\link{c_hat}}, \code{\link{mb.gof.test}},
\code{\link{Nmix.gof.test}}, \code{\link[unmarked]{parboot}},
\code{\link{kappa}}, \code{\link{rcond}}
}

\examples{
##N-mixture model example modified from ?pcount
\dontrun{
require(unmarked)
##single season
data(mallard)
mallardUMF <- unmarkedFramePCount(mallard.y, siteCovs = mallard.site,
                                  obsCovs = mallard.obs)
##run model
fm.mallard <- pcount(~ ivel+ date + I(date^2) ~ length + elev + forest,
                     mallardUMF, K=30)

##compute condition number
extractCN(fm.mallard)

##compare against 'kappa'
kappa(fm.mallard@opt$hessian, exact = TRUE)
detach(package:unmarked)
}
}


\keyword{models}
