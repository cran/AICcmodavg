### R code from vignette source 'AICcmodavg.Rnw'

###################################################
### code chunk number 1: AICcmodavg.Rnw:37-38
###################################################
options(width=70, continue = "  ")


###################################################
### code chunk number 2: import
###################################################
library(AICcmodavg)
data(dry.frog)


###################################################
### code chunk number 3: subData
###################################################
##extract only first 7 columns
frog <- dry.frog[, 1:7]
##first lines
head(frog)
##structure of data frame
str(frog)


###################################################
### code chunk number 4: na
###################################################
any(is.na(frog))


###################################################
### code chunk number 5: centInitialMass
###################################################
##center initial mass
frog$InitMass_cent <- frog$Initial_mass - mean(frog$Initial_mass)


###################################################
### code chunk number 6: InitialMass2
###################################################
frog$InitMass2 <- frog$InitMass_cent^2


###################################################
### code chunk number 7: checkDiag
###################################################
##run global model
global <- lm(Mass_lost ~ InitMass_cent + InitMass2 + Substrate + Shade, 
             data = frog)
par(mfrow = c(2, 2))
plot(global)


###################################################
### code chunk number 8: logMass
###################################################
frog$logMass_lost <- log(frog$Mass_lost + 1) #adding 1 due to presence of 0's


###################################################
### code chunk number 9: checkDiag2
###################################################
##run global model
global.log <- lm(logMass_lost ~ InitMass_cent + InitMass2 + Substrate + Shade, 
                 data = frog)
par(mfrow = c(2, 2))
plot(global.log)


###################################################
### code chunk number 10: fitCands
###################################################
m.null <- lm(logMass_lost ~ 1,
             data = frog)
m.shade <- lm(logMass_lost ~ Shade, 
              data = frog)
m.substrate <- lm(logMass_lost ~ Substrate, 
                  data = frog)
m.shade.substrate <- lm(logMass_lost ~ Shade + Substrate, 
                        data = frog)
m.null.mass <- lm(logMass_lost ~ InitMass_cent + InitMass2,
                  data = frog)
m.shade.mass <- lm(logMass_lost ~ InitMass_cent + InitMass2 + Shade, 
                   data = frog)
m.substrate.mass <- lm(logMass_lost ~ InitMass_cent + InitMass2 + Substrate, 
                       data = frog)
m.global.mass <- global.log


###################################################
### code chunk number 11: storeList
###################################################
##store models in named list
Cand.models <- list("null" = m.null, "shade" = m.shade, 
                    "substrate" = m.substrate, 
                    "shade + substrate" = m.shade.substrate, 
                    "mass" = m.null.mass, "mass + shade" = m.shade.mass, 
                    "mass + substrate" = m.substrate.mass, 
                    "global" = m.global.mass)


###################################################
### code chunk number 12: modTableAICc
###################################################
selectionTable <- aictab(cand.set = Cand.models)
selectionTable


###################################################
### code chunk number 13: modTableAIC
###################################################
aictab(Cand.models, second.ord = FALSE)


###################################################
### code chunk number 14: exportTable3 (eval = FALSE)
###################################################
## library(xtable)
## print(xtable(selectionTable, caption = "Model selection table on frog mass lost.", 
##              label = "tab:selection"),
##       include.rownames = FALSE, caption.placement = "top")


###################################################
### code chunk number 15: exportTable4
###################################################
library(xtable)
print(xtable(selectionTable, caption = "Model selection table on frog mass lost.", 
             label = "tab:selection"),
      include.rownames = FALSE, caption.placement = "top", )


###################################################
### code chunk number 16: confSet
###################################################
##confidence set of models
confset(cand.set = Cand.models)


###################################################
### code chunk number 17: evidence
###################################################
##evidence ratios
evidence(aic.table = selectionTable)


###################################################
### code chunk number 18: evidenceSilent
###################################################
evRatio <- evidence(selectionTable)


###################################################
### code chunk number 19: evidence2
###################################################
##compare "substrate" vs "shade"
evidence(selectionTable, model.high = "substrate",
         model.low = "shade")


###################################################
### code chunk number 20: evidence2Silent
###################################################
##compare "substrate" vs "shade"
evRatio2 <- evidence(selectionTable, model.high = "substrate",
                     model.low = "shade")


###################################################
### code chunk number 21: evidenceNull
###################################################
evidence(selectionTable, model.high = "global",
         model.low = "null")


###################################################
### code chunk number 22: confint
###################################################
confint(m.global.mass)


###################################################
### code chunk number 23: modavg
###################################################
modavg(cand.set = Cand.models, parm = "Shade")


###################################################
### code chunk number 24: modavg2
###################################################
modavgShade <- modavg(cand.set = Cand.models, parm = "Shade")


###################################################
### code chunk number 25: coef
###################################################
coef(m.global.mass)


###################################################
### code chunk number 26: substrateSPHAG
###################################################
modavg(Cand.models, parm = "SubstrateSPHAGNUM")


###################################################
### code chunk number 27: substrateSPHAG2
###################################################
modavgSphag <- modavg(Cand.models, parm = "SubstrateSPHAGNUM")


###################################################
### code chunk number 28: modavg
###################################################
modavgShrink(cand.set = Cand.models, parm = "Shade")


###################################################
### code chunk number 29: substrateSPHAGShrink
###################################################
modavgShrink(Cand.models, parm = "SubstrateSPHAGNUM")


###################################################
### code chunk number 30: shadePred
###################################################
##data frame to make predictions
##all variables are held constant, except Shade
predData <- data.frame(InitMass_cent = c(0, 0),
                       InitMass2 = c(0, 0),
                       Substrate = factor("SOIL", 
                                          levels = levels(frog$Substrate)),
                       Shade = c(0, 1))
##predictions from global model
predict(m.global.mass, newdata = predData, se.fit = TRUE)
##predictions from null model
predict(m.null, newdata = predData, se.fit = TRUE)


###################################################
### code chunk number 31: extractX
###################################################
extractX(cand.set = Cand.models)


###################################################
### code chunk number 32: modavgPred
###################################################
modavgPred(cand.set = Cand.models, newdata = predData)


###################################################
### code chunk number 33: modavgPredSub
###################################################
##data frame holding all variables constant, except Substrate
predSub <- data.frame(InitMass_cent = c(0, 0, 0),
                      InitMass2 = c(0, 0, 0),
                      Substrate = factor(c("PEAT", "SOIL", "SPHAGNUM"),
                                         levels = levels(frog$Substrate)),
                      Shade = c(1, 1, 1))
##model-average predictions
predsMod <- modavgPred(Cand.models, newdata = predSub)
predsMod


###################################################
### code chunk number 34: checkContent
###################################################
##check content of object
str(predsMod)


###################################################
### code chunk number 35: savePreds
###################################################
##add predictions, lower CL, and upper CL
predSub$fit <- predsMod$mod.avg.pred
predSub$low95 <- predsMod$lower.CL
predSub$upp95 <- predsMod$upper.CL


###################################################
### code chunk number 36: plotPreds
###################################################
##create vector for X axis
predSub$xvals <- c(0.25, 0.5, 0.75)
##create empty box
plot(fit ~ xvals, 
     data = predSub,
     xlim = c(0, 1),
     ylim = range(low95, upp95),
     xlab = "Substrate type",
     ylab = "Predicted mass lost (log of mass in g)",
     xaxt = "n",
     cex.axis = 1.2, 
     cex.lab = 1.2)
##add x axis
axis(side = 1, at = predSub$xvals,
     labels = c("Peat", "Soil", "Sphagnum"),
     cex.axis = 1.2)
##add CI's
segments(x0 = predSub$xvals, x1 = predSub$xvals,
         y0 = predSub$low95, predSub$upp95)


###################################################
### code chunk number 37: compGroups
###################################################
predComp <- data.frame(InitMass_cent = c(0, 0),
                       InitMass2 = c(0, 0),
                       Substrate = factor(c("PEAT", "SPHAGNUM"),
                                          levels = levels(frog$Substrate)),
                       Shade = c(1, 1))
##model-average predictions
modavgEffect(Cand.models, newdata = predComp)


###################################################
### code chunk number 38: customAICc
###################################################
##log-likelihoods
modL <- c(-225.4180, -224.0697, -225.4161)
##number of parameters
modK <- c(2, 3, 3)
##model selection
outTab <- aictabCustom(logL = modL,
                       K = modK,
                       modnames = c("null", "phi(SVL)p(.)",
                                    "phi(Road)p(.)"),
                       nobs = 621)


###################################################
### code chunk number 39: evRatioCustom (eval = FALSE)
###################################################
## evidence(outTab, model.high = "phi(SVL)p(.)",
##          model.low = "phi(Road)p(.)")


###################################################
### code chunk number 40: evRatioCustom
###################################################
evRatioCust <- evidence(outTab, model.high = "phi(SVL)p(.)",
                        model.low = "phi(Road)p(.)")


###################################################
### code chunk number 41: estSE
###################################################
##survival estimates with road mitigation
modEst <- c(0.1384450, 0.1266030, 0.1378745)
##SE's of survival estimates with road mitigation
modSE <- c(0.03670327, 0.03347475, 0.03862634)


###################################################
### code chunk number 42: customModavg
###################################################
##model-averaged survival with road mitigation
modavgCustom(logL = modL,
             K = modK,
             modnames = c("null", "phi(SVL)p(.)",
                          "phi(Road)p(.)"),
             estimate = modEst,
             se = modSE,
             nobs = 621)


###################################################
### code chunk number 43: customModavg2
###################################################
##survival estimates without road mitigation
modEst2 <- c(0.1384450, 0.1266030, 0.1399727)

##SE's of survival estimates without road mitigation
modSE2 <- c(0.03670327, 0.03347475, 0.04981635)

##model-averaged survival
modavgCustom(logL = modL,
             K = modK,
             modnames = c("null", "phi(SVL)p(.)",
                          "phi(Road)p(.)"),
             estimate = modEst2,
             se = modSE2,
             nobs = 621)


