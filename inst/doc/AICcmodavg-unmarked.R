### R code from vignette source 'AICcmodavg-unmarked.Rnw'

###################################################
### code chunk number 1: AICcmodavg-unmarked.Rnw:37-38
###################################################
options(width=70, continue = "  ")


###################################################
### code chunk number 2: loadPackage
###################################################
##load package
library(AICcmodavg)
##load data frame
data(bullfrog)


###################################################
### code chunk number 3: checkBullfrog
###################################################
##check data structure
str(bullfrog)
##first rows
head(bullfrog)


###################################################
### code chunk number 4: formatData
###################################################
##extract detections
yObs <- bullfrog[, c("V1", "V2", "V3", "V4", "V5", "V6", "V7")]

##extract site variables
siteVars <- bullfrog[, c("Location", "Reed.presence")]

##extract observation variables
##centered sampling effort on each visit
effort <- bullfrog[, c("Effort1", "Effort2", "Effort3", "Effort4",
                       "Effort5", "Effort6", "Effort7")]
##survey type (0 = call survey, 1 = minnow trap)
type <- bullfrog[, c("Type1", "Type2", "Type3", "Type4", "Type5",
                     "Type6", "Type7")]


###################################################
### code chunk number 5: loadFormat
###################################################
##load package
library(unmarked)
##format data
bfrogData <- unmarkedFrameOccu(y = yObs,
                               siteCovs = siteVars,
                               obsCovs = list(Type = type, Effort = effort))


###################################################
### code chunk number 6: summary1
###################################################
summary(bfrogData)


###################################################
### code chunk number 7: detHist
###################################################
detHist(bfrogData)


###################################################
### code chunk number 8: fitOccu
###################################################
##null model
m1 <- occu(~ 1 ~ 1, data = bfrogData)
##p varies with survey type and effort, occupancy is constant
m2 <- occu(~ Type + Effort ~ 1, data = bfrogData)
##p constant, occupancy varies with reed presence
m3 <- occu(~ 1 ~ Reed.presence, data = bfrogData)
##global model
m4 <- occu(~ Type + Effort ~ Reed.presence, data = bfrogData)


###################################################
### code chunk number 9: checkOut
###################################################
summary(m4)
summaryOD(m4, out.type = "confint")
summaryOD(m4, out.type = "nhst")


###################################################
### code chunk number 10: createList
###################################################
bfrogMods <- list("null" = m1, "psidot.pTypeEffort" = m2,
                  "psiReed.pdot" = m3, 
                  "psiReed.pTypeEffort" = m4)


###################################################
### code chunk number 11: checkConv
###################################################
##check convergence for a single model
checkConv(m1)
##extract values across all models
sapply(bfrogMods, checkConv)


###################################################
### code chunk number 12: extractCN
###################################################
##extract condition number of single model
extractCN(m1)
##extract condition across all models
sapply(bfrogMods, extractCN)


###################################################
### code chunk number 13: largeSE
###################################################
##check highest SE in single model
checkParms(m1)
##check highest SE across all models
lapply(bfrogMods, checkParms)


###################################################
### code chunk number 14: LRT
###################################################
##compare global model vs null
anovaOD(mod.simple = m1, mod.complex = m3)


###################################################
### code chunk number 15: gof (eval = FALSE)
###################################################
## ##this takes 226 min. using 2 cores
## gof <- mb.gof.test(mod = m4, nsim = 10000, parallel = TRUE, ncores = 2)
## gof
## save(gof, file = "gofMod3.Rdata")


###################################################
### code chunk number 16: gof2
###################################################
load("gofMod3.Rdata")
gof
p.value <- sum(gof$t.star >= gof$chi.square)/gof$nsim
if (p.value == 0) {
    p.display <- paste("<", round(1/gof$nsim, digits = 4))
} else {
    p.display <- paste("=", round(p.value, digits = 4))
}
hist(gof$t.star, 
     main = "Bootstrapped MacKenzie and Bailey fit statistic (10 000 samples)", 
     xlim = range(c(gof$t.star, gof$chi.square)), 
     xlab = paste("Simulated statistic (observed = ", 
                  round(gof$chi.square, digits = 2), ")", sep = ""), 
     cex.axis = 1.2, cex.lab = 1.2, cex.main = 1.2)
title(main = bquote(paste(italic(P), " ", .(p.display))), 
      line = 0.5, cex.main = 1.2)
abline(v = gof$chi.square, lty = "dashed", 
       col = "red")


###################################################
### code chunk number 17: summaryOD2
###################################################
##compare inferences 
summaryOD(m3)
summaryOD(m3, c.hat = 1.08)


###################################################
### code chunk number 18: aic
###################################################
##when no overdispersion is present
outTab <- aictab(cand.set = bfrogMods) 
##accounting for overdispersion
outTabC <- aictab(cand.set = bfrogMods, c.hat = 1.08)
outTab
outTabC


###################################################
### code chunk number 19: evidenceRatio
###################################################
##evidence ratio between top-ranked model vs second-ranked model
evidence(aic.table = outTabC)


###################################################
### code chunk number 20: modavgShrink
###################################################
##model-averaged estimate of reed presence - shrinkage estimator
estReed <- modavgShrink(cand.set = bfrogMods, 
                        parm = "Reed.presence", parm.type = "psi",
                        c.hat = 1.08)
estReed


###################################################
### code chunk number 21: modavgShrink2
###################################################
estType <- modavgShrink(cand.set = bfrogMods, 
                        parm = "Type", parm.type = "detect",
                        c.hat = 1.08)
estType
estEffort <- modavgShrink(cand.set = bfrogMods, 
                          parm = "Effort", parm.type = "detect",
                          c.hat = 1.08)
estEffort


###################################################
### code chunk number 22: checkXpsi
###################################################
##variables on psi
extractX(cand.set = bfrogMods, parm.type = "psi")
##variables on p
extractX(cand.set = bfrogMods, parm.type = "detect")


###################################################
### code chunk number 23: predReed
###################################################
reedFrame <- data.frame(Reed.presence = c(0, 1))


###################################################
### code chunk number 24: predReed2
###################################################
outReed <- modavgPred(cand.set = bfrogMods, newdata = reedFrame, 
                      parm.type = "psi", c.hat = 1.08)
outReed


###################################################
### code chunk number 25: storeReed
###################################################
##store predictions and confidence intervals in data frame
reedFrame$fit <- outReed$mod.avg.pred
reedFrame$low95 <- outReed$lower.CL
reedFrame$upp95 <- outReed$upper.CL


###################################################
### code chunk number 26: plotReed
###################################################
##create plot
xvals <- c(0.2, 0.4)
plot(fit ~ xvals, data = reedFrame,
     ylab = "Probability of occupancy",
     xlab = "Presence of reed",
     ylim = c(0, 1),
     cex = 1.2, cex.axis = 1.2, cex.lab = 1.2,
     xlim = c(0, 0.6),
     xaxt = "n")
#add x axis
axis(side = 1, at = xvals, 
     labels = c("absent", "present"), 
     cex.axis = 1.2)
##add error bars
segments(x0 = xvals, y0 = reedFrame$low95,
         x1 = xvals, y1 = reedFrame$upp95)


###################################################
### code chunk number 27: predType
###################################################
##vary Type, hold Effort constant at its mean 
typeFrame <- data.frame(Type = c(0, 1), Effort = 0)
##model-averaged predictions
outType <- modavgPred(cand.set = bfrogMods, newdata = typeFrame, 
                      parm.type = "detect", c.hat = 1.08)
outType


###################################################
### code chunk number 28: plotType
###################################################
##store predictions and confidence intervals in data frame
typeFrame$fit <- outType$mod.avg.pred
typeFrame$low95 <- outType$lower.CL
typeFrame$upp95 <- outType$upper.CL
##create plot
xvals <- c(0.2, 0.4)
plot(fit ~ xvals, data = typeFrame,
     ylab = "Detection probability",
     xlab = "Survey type",
     ylim = c(0, 1),
     cex = 1.2, cex.axis = 1.2, cex.lab = 1.2,
     xlim = c(0, 0.6),
     xaxt = "n")
#add x axis
axis(side = 1, at = xvals, 
     labels = c("call survey", "minnow trapping"), 
     cex.axis = 1.2)
##add error bars
segments(x0 = xvals, y0 = typeFrame$low95,
         x1 = xvals, y1 = typeFrame$upp95)


###################################################
### code chunk number 29: extractEffort
###################################################
##extract centered values of sampling effort
effort <- bfrogData@obsCovs$Effort

##create a series of 30 values to plot
Effort.cent <- seq(from = min(effort), to = max(effort),
                   length.out = 30)
##back-transform values to original scale of variable
Effort.mean <- 8.67 #mean of original variable see ?bullfrog
Effort.orig <- Effort.cent + Effort.mean


###################################################
### code chunk number 30: predEffort
###################################################
##note that all variables on the parameter must appear here
pred.dataEffort <- data.frame(Effort.orig = Effort.orig,
                              Effort = Effort.cent, #centered variable
                              Type = 1) 
#Recall that \texttt{Type} was coded 1 (minnow trap) or 0 (call survey)

##compute model-averaged predictions with modavgPred on probability scale
out.predsEffort <- modavgPred(cand.set = bfrogMods,
                              newdata = pred.dataEffort, parm.type = "detect",
                              type = "response", c.hat = 1.08)


###################################################
### code chunk number 31: addPreds
###################################################
##add predictions to data set to keep everything in the same place
pred.dataEffort$fit <- out.predsEffort$mod.avg.pred
pred.dataEffort$se.fit <- out.predsEffort$uncond.se
pred.dataEffort$low95 <- out.predsEffort$lower.CL
pred.dataEffort$upp95 <- out.predsEffort$upper.CL


###################################################
### code chunk number 32: plotEffort
###################################################
##create plot

##plot
plot(fit ~ Effort.orig,
     ylab = "Detection probability",
     xlab = "Sampling effort",
     ylim = c(0, 1),
     type = "l", 
     cex = 1.2, cex.lab = 1.2, cex.axis = 1.2,
     data = pred.dataEffort)

##add 95% CI around predictions
lines(low95 ~ Effort.orig, data = pred.dataEffort,
      lty = "dashed")
lines(upp95 ~ Effort.orig, data = pred.dataEffort,
      lty = "dashed")


###################################################
### code chunk number 33: xtable0
###################################################
library(xtable)
xtable(summaryOD(m4))


###################################################
### code chunk number 34: xtable1
###################################################
xtable(outTabC)


###################################################
### code chunk number 35: table2
###################################################
xtable(estReed)


###################################################
### code chunk number 36: table3
###################################################
xtable(detHist(m3))


###################################################
### code chunk number 37: table4
###################################################
xtable(mb.chisq(m3))


###################################################
### code chunk number 38: xtableOptions
###################################################
#add caption, suppress log-likelihood, and include cumulative Akaike weight
print(xtable(outTabC, 
             caption = "Model selection accounting for overdispersion in the bullfrog data.", 
             include.LL = FALSE, include.Cum.Wt = TRUE),
      caption.placement = "top", include.rownames = FALSE)


