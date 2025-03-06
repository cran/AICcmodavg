\name{fat}
\Rdversion{1.1}
\alias{fat}
\docType{data}
\title{
Fat Data and Body Measurements
}
\description{
This data set illustrates the relationship between body measurements and
body fat in 252 males aged between 21 and 81 years. 
}
\usage{data(fat)}
\format{
  A data frame with 252 rows and 26 variables.
  \describe{
    \item{\code{Obs}}{observation number.}
    \item{\code{Perc.body.fat.Brozek}}{percent body fat using Brozek's
      equation, i.e., \eqn{457/Density - 414.2}.}
    \item{\code{Perc.body.fat.Siri}}{percent body fat using Siri's
      equation, i.e., \eqn{495/Density - 450}.}
    \item{\code{Density}}{density (\eqn{\frac{g}{cm^3}}).}
    \item{\code{Age}}{age (years).}
    \item{\code{Weight}}{weight (lbs).}
    \item{\code{Height}}{height (inches).}
    \item{\code{Adiposity.index}}{adiposity index computed as \eqn{Weight/Height^2} (\eqn{\frac{kg}{m^2}}).}
    \item{\code{Fat.free.weight}}{fat free weight computed as \eqn{(1 -
	Brozek's percent body fat) * Weight} (lbs).}
    \item{\code{Neck.circ}}{neck circumference (cm).}
    \item{\code{Chest.circ}}{chest circumference (cm).}
    \item{\code{Abdomen.circ}}{abdomen circumference (cm) measured at
      the umbilicus and level with the iliac crest.} 
    \item{\code{Hip.circ}}{hip circumference (cm).}
    \item{\code{Thigh.circ}}{thigh circumference (cm).}
    \item{\code{Knee.circ}}{knee circumference (cm).}
    \item{\code{Ankle.circ}}{ankle circumference (cm).}
    \item{\code{Biceps.circ}}{extended biceps circumference (cm).}
    \item{\code{Forearm.circ}}{forearm circumference (cm).}
    \item{\code{Wrist.circ}}{wrist circumference (cm).}
    \item{\code{inv.Density}}{inverse of density (\eqn{\frac{cm^3}{g}}).}
    \item{\code{z1}}{log of weight divided by log of height (allometric measure).}
    \item{\code{z2}}{abdomen circumference divided by chest
      circumference (beer gut factor).}
    \item{\code{z3}}{index based on knee, wrist, and ankle circumference
      relative to height (\eqn{\frac{(Knee.circ * Wrist.circ * Ankle.circ)^(1/3)}{Height}}).}
    \item{\code{z4}}{fleshiness index based on biceps, thigh, forearm,
      knee, wrist, and ankle circumference (\eqn{\frac{Biceps.circ *
	Thigh.circ * Forearm.circ}{Knee.circ * Wrist.circ *
	Ankle.circ}}).}
    \item{\code{z5}}{age standardized to zero mean and unit variance.}
    \item{\code{z6}}{square of standardized age.}
  }
}
\details{
Burnham and Anderson (2002, p. 268) use this data set to show model
selection uncertainty in the context of all possible combinations of
explanatory variables.  The data are originally from Penrose et
al. (1985) who used only the first 143 cases of the 252 observations in
the data set.  Johnson (1996) later used these data as an example of
multiple regression.  Note that observation number 42 originally had an
erroneous height of 29.5 inches and that this value was changed to 69.5
inches.

Burnham and Anderson (2002, p. 274) created six indices based on the
original measurements (i.e., z1 -- z6).  Although Burnham and Anderson
(2002) indicate that the fleshiness index (\code{z4}) involved the cubic
root in the equation, the result table for the full model on p. 276
suggests that the index did not include the cubic root for \code{z4}.
The latter is the version of \code{z4} used in the data set here.
}
\source{
  Burnham, K. P., Anderson, D. R. (2002) \emph{Model Selection and
    Multimodel Inference: a practical information-theoretic
    approach}. Second edition. Springer: New York.

  Johnson, J. W. (1996). Fitting percentage of body fat to simple body
  measurements. \emph{Journal of Statistics Education} \bold{4} [online].

  Penrose, K., Nelson, A., Fisher, A. (1985) Generalized body composition
  prediction equation for men using simple measurement techniques
  \emph{Medicine and Science in Sports and Exercise} \bold{17}, 189.
}
\examples{
data(fat)
str(fat)
}
\keyword{datasets}
