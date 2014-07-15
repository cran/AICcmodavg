\name{bullfrog}
\Rdversion{1.1}
\alias{bullfrog}
\docType{data}
\title{
Bullfrog occupancy data and common reed invasion
}
\description{
This is a data set from Mazerolle et al. (2014) on the occupancy of
Bullfrogs (\emph{Lithobates catesbeianus}) in 50 wetlands sampled in
2009 in the area of Montreal, QC.
}
\usage{data(bullfrog)}
\format{
  A data frame with 50 observations on the following 23 variables.
  \describe{
    \item{\code{Location}}{a factor with a unique identifier for each
      wetland.}
    \item{\code{Reed.presence}}{a binary variable, either 1 (reed
      present) or 0 (reed absent).}
    \item{\code{V1}}{a binary variable for detection (1) or non
      detection (0) of bullfrogs during the first survey.} 
    \item{\code{V2}}{a binary variable for detection (1) or non
      detection (0) of bullfrogs during the second survey.} 
    \item{\code{V3}}{a binary variable for detection (1) or non
      detection (0) of bullfrogs during the third survey.}
    \item{\code{V4}}{a binary variable for detection (1) or non
      detection (0) of bullfrogs during the fourth survey.}
    \item{\code{V5}}{a binary variable for detection (1) or non
      detection (0) of bullfrogs during the fifth survey.}
    \item{\code{V6}}{a binary variable for detection (1) or non
      detection (0) of bullfrogs during the sixth survey.}
    \item{\code{V7}}{a binary variable for detection (1) or non
      detection (0) of bullfrogs during the seventh survey.}
    \item{\code{Effort1}}{a numeric variable for the centered
      number of sampling stations during the first survey.}
    \item{\code{Effort2}}{a numeric variable for the centered
      number of sampling stations during the second survey.}
    \item{\code{Effort3}}{a numeric variable for the centered
      number of sampling stations during the third survey.}
    \item{\code{Effort4}}{a numeric variable for the centered
      number of sampling stations during the fourth survey.}
    \item{\code{Effort5}}{a numeric variable for the centered
      number of sampling stations during the fifth survey.}
    \item{\code{Effort6}}{a numeric variable for the centered
      number of sampling stations during the sixth survey.}
    \item{\code{Effort7}}{a numeric variable for the centered
      number of sampling stations during the seventh survey.}
    \item{\code{Type1}}{a binary variable to identify the survey type,
      either minnow trap (1) or call survey (0) during the first
      sampling occasion.}
    \item{\code{Type2}}{a binary variable to identify the survey type,
      either minnow trap (1) or call survey (0) during the second
      sampling occasion.}
    \item{\code{Type3}}{a binary variable to identify the survey type,
      either minnow trap (1) or call survey (0) during the third
      sampling occasion.}
    \item{\code{Type4}}{a binary variable to identify the survey type,
      either minnow trap (1) or call survey (0) during the fourth
      sampling occasion.}
    \item{\code{Type5}}{a binary variable to identify the survey type,
      either minnow trap (1) or call survey (0) during the fifth
      sampling occasion.}
    \item{\code{Type6}}{a binary variable to identify the survey type,
      either minnow trap (1) or call survey (0) during the sixth
      sampling occasion.}
    \item{\code{Type7}}{a binary variable to identify the survey type,
      either minnow trap (1) or call survey (0) during the seventh
      sampling occasion.}
  }
}
\details{
This data set is used to illustrate single-species single-season
occupancy models (MacKenzie et al. 2002). 
}
\source{
  MacKenzie, D. I., Nichols, J. D., Lachman, G. B., Droege, S., Royle,
  J. A., Langtimm, C. A. (2002). Estimating site occupancy rates when
  detection probabilities are less than one. \emph{Ecology} \bold{83},
  2248--2255. 

  Mazerolle, M. J., Perez, A., Brisson, J. (2014) Common reed
  (\emph{Phragmites australis}) invasion and amphibian distribution in
  freshwater wetlands. \emph{Wetlands Ecology and Management} \bold{22},
  325--340. 
}
\examples{
data(bullfrog)
str(bullfrog)
}
\keyword{datasets}
