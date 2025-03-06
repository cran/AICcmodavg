##summarize time to detection data
detTime <- function(object, plot.time = TRUE, plot.seasons = FALSE,
                    cex.axis = 1, cex.lab = 1, cex.main = 1, ...){
    UseMethod("detTime", object)
}



detTime.default <- function(object, plot.time = TRUE, plot.seasons = FALSE,
                            cex.axis = 1, cex.lab = 1, cex.main = 1, ...){
  stop("\nFunction not yet defined for this object class\n")
}



##for ummarkedFrameOccuTTD
detTime.unmarkedFrameOccuTTD <- function(object, plot.time = TRUE, plot.seasons = FALSE,
                                         cex.axis = 1, cex.lab = 1, cex.main = 1, ...) {

    ##extract data
    yMat <- object@y
    nsites <- nrow(yMat)
    n.seasons <- object@numPrimary
    nvisits <- ncol(yMat)
    ##visits per season
    n.visits.season <- nvisits/n.seasons
    n.seasons.adj <- n.seasons #total number of plots fixed to 11 or 12, depending on plots requested
    seasonNames <- paste("season", 1:n.seasons, sep = "")
    surveyLength <- object@surveyLength
    

    if(plot.seasons && n.seasons == 1) {
        warning("\nCannot plot data across seasons with only 1 season of data: reset to plot.seasons = FALSE\n")
        plot.seasons <- FALSE
    }
    
    if(plot.time && !plot.seasons) {
        nRows <- 1
        nCols <- 1
        
        ##reset graphics parameters and save in object
        oldpar <- par(mfrow = c(nRows, nCols))
    }

    
    if(!plot.time && plot.seasons) {
        ##determine arrangement of plots in matrix
        if(plot.seasons && n.seasons >= 12) {
            n.seasons.adj <- 12
            warning("\nOnly first 12 seasons are plotted\n")
        }
    
        if(plot.seasons && n.seasons.adj <= 12) {

            ##if n.seasons < 12
            ##if 12, 11, 10 <- 4 x 3
            ##if 9, 8, 7 <- 3 x 3
            ##if 6, 5 <- 3 x 2
            ##if 4 <- 2 x 2
            ##if 3 <- 3 x 1
            ##if 2 <- 2 x 1
        
            if(n.seasons.adj >= 10) {
                nRows <- 4
                nCols <- 3
            } else {
            
                if(n.seasons.adj >= 7) {
                    nRows <- 3
                    nCols <- 3
                } else {
                
                    if(n.seasons.adj >= 5) {
                        nRows <- 3
                        nCols <- 2
                    } else {
                        if(n.seasons.adj == 4) {
                            nRows <- 2
                            nCols <- 2
                        } else {
                            if(n.seasons.adj == 3) {
                                nRows <- 3
                                nCols <- 1
                            } else {
                                nRows <- 2
                                nCols <- 1
                            }
                        }
                    }
                }
            }
            ##reset graphics parameters and save in object
            oldpar <- par(mfrow = c(nRows, nCols))
        }
    }
    



    ##if both plots for seasons and combined are requested
    if(plot.time && plot.seasons) {
        
        if(plot.seasons && n.seasons >= 12) {
            n.seasons.adj <- 11
            warning("\nOnly first 11 seasons are plotted\n")
        }

        if(plot.seasons && n.seasons.adj <= 11) {
            
            if(n.seasons.adj >= 9) {
                nRows <- 4
                nCols <- 3
            } else {
                
                if(n.seasons.adj >= 6) {
                    nRows <- 3
                    nCols <- 3
                } else {
                    
                    if(n.seasons.adj >= 4) {
                        nRows <- 3
                        nCols <- 2
                    } else {
                        if(n.seasons.adj == 3) {
                            nRows <- 2
                            nCols <- 2
                        } else {
                            if(n.seasons.adj == 2) {
                                nRows <- 3
                                nCols <- 1
                            }
                        }
                    }
                }
            }
        }
        
        ##reset graphics parameters and save in object
        oldpar <- par(mfrow = c(nRows, nCols))
    }
        
        
    ##combine all seasons
    ##Censoring distance
    censoredDist.full <- surveyLength
    uniqueDist.full <- unique(as.vector(censoredDist.full))
    
    ##determine data that were censored
    uncensoredIndex.full <- yMat < censoredDist.full
    uncensoredData.full <- yMat[uncensoredIndex.full]
    ncensored.full <- nsites - sum(uncensoredIndex.full)

    
    if(plot.time) {
    
        ##check that maximum times are the same for all sites
        if(n.seasons == 1) {
            
            if(length(uniqueDist.full) == 1) {
                main.title <- paste("Distribution of time to detection (survey length: ", uniqueDist.full, " min.)",
                                    sep = "")
            } else {
                main.title  <- paste("Distribution of time to detection (survey length: ", min(uniqueDist.full), "-",
                                     max(uniqueDist.full), " min.)",
                                     sep = "")
            }
        }
        
        if(n.seasons > 1) {
            
            if(length(uniqueDist.full) == 1) {
                main.title <- paste("Distribution of time to detection (", n.seasons, " seasons, survey length: ",
                                    uniqueDist.full, " min.)",
                                    sep = "")
            } else {
                main.title  <- paste("Distribution of time to detection (", n.seasons, " seasons, survey length: ",
                                     min(uniqueDist.full), "-",
                                     max(uniqueDist.full), " min.)",
                                     sep = "")
            }
        }
        
        hist(uncensoredData.full, xlim = c(0, uniqueDist.full),
             xlab = "Time to detection (min.)",
             main = main.title, cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
    }
        

    ##quantiles for entire data
    time.table.full <- quantile(uncensoredData.full, na.rm = TRUE)

    ##store data for each season
    
    ##columns for each season
    col.seasons <- seq(1, nvisits, by = n.visits.season)

    ##list to store raw data
    yMat.seasons <- vector(mode = "list", length = n.seasons)
    names(yMat.seasons) <- seasonNames
    minusOne <- n.visits.season - 1
        
    ##list to store quantiles excluding censored times
    time.table.seasons <- vector("list", n.seasons)
    names(time.table.seasons) <- seasonNames

    ##list of uncensored observations
    uncensoredData.seasons <- vector("list", n.seasons)
    names(uncensoredData.seasons) <- seasonNames
    
    ##vector of censored observations
    censored.seasons <- vector("numeric", n.seasons)
    names(censored.seasons) <- seasonNames

    ##list of unique values of maximum effort
    uniqueDist.seasons <- vector("list", n.seasons)
    names(uniqueDist.seasons) <- seasonNames

    ##list of maximum effort
    censoredDist.seasons <- vector("list", n.seasons)
    names(censoredDist.seasons) <- seasonNames
    
    ##iterate over each season
    for(i in 1:n.seasons) {
        ##extract yMat for each season
        yMat1 <- yMat[, col.seasons[i]:(col.seasons[i]+minusOne), drop = FALSE]
        yMat.seasons[[i]] <- yMat1
            
        ##Censoring values
        censoredDist <- surveyLength[, col.seasons[i]:(col.seasons[i]+minusOne), drop = FALSE]
        censoredDist.seasons[[i]] <- censoredDist
        uniqueDist.seasons[[i]] <- unique(as.vector(censoredDist))
        
        ##determine data that were censored
        uncensoredIndex <- yMat1 < censoredDist
        uncensoredData <- yMat1[uncensoredIndex]
        uncensoredData.seasons[[i]] <- uncensoredData
        ncensored <- nsites - sum(uncensoredIndex)
        
        ##summarize times per season
        time.quantiles <- quantile(uncensoredData, na.rm = TRUE)
        time.table.seasons[[i]] <- time.quantiles
        censored.seasons[i] <- ncensored
        
    }

    ##check if any seasons were not sampled
    y.seasonsNA <- sapply(yMat.seasons, FUN = function(i) all(is.na(i)))
    
    if(plot.seasons) {
        
        for(i in 1:n.seasons) {
            
            if(length(uniqueDist.seasons[[i]]) == 1) {
                main.title <- paste("Distribution of time to detection (season ", i,
                                    ", survey length: ", uniqueDist.seasons[[i]],
                                    " min.)", sep = "")
            } else {
                main.title <- paste("Distribution of time to detection (season ", i,
                                    ", survey length: ", min(uniqueDist.seasons[[i]]), "-",
                                    max(uniqueDist.seasons[[i]]), " min.)",
                                    sep = "")
            }
            
                
            ##create histogram
            ##check for missing season
            if(y.seasonsNA[i]) {next} #skip to next season if current season not sampled
            hist(uncensoredData.seasons[[i]], xlim = c(0, uniqueDist.seasons[[i]]),
                 xlab = "Time to detection (min.)",
                 main = main.title, cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
        }
    }
        
    if(n.seasons == 1) {
        ##for each season, determine frequencies
        out.freqs <- matrix(data = NA, ncol = 2, nrow = n.seasons)
        colnames(out.freqs) <- c("sampled", "detected")
        rownames(out.freqs) <- "Season-1"
        
        ##sequence of visits
        for(i in 1:n.seasons) {
            ySeason <- yMat.seasons[[i]]
            censored <- censoredDist.seasons[[i]]
            
            uncensoredObs <- matrix(NA, ncol = n.visits.season,
                                    nrow = nsites)
            
            for(j in 1:ncol(ySeason)){
                ##observations
                uncensoredObs[, j] <- ySeason[, j] < censored[, j]
            }
            ##determine proportion of sites with at least 1 detection
            det.sum <- rowSums(uncensoredObs, na.rm = TRUE)
            
            ##check sites with observed detections and deal with NA's
            sum.rows <- rowSums(uncensoredObs, na.rm = TRUE)
            is.na(sum.rows) <- rowSums(is.na(uncensoredObs)) == ncol(ySeason)
            
            ##number of sites sampled
            out.freqs[i, 1] <- sum(!is.na(sum.rows))
            out.freqs[i, 2] <- sum(det.sum)
            
        }
              
    
        ##create a matrix with proportion of sites with colonizations
        ##and extinctions based on raw data
        out.props <- matrix(NA, nrow = nrow(out.freqs), ncol = 1)
        colnames(out.props) <- "naive.occ"
        rownames(out.props) <- rownames(out.freqs)
        out.props[, 1] <- out.freqs[, 2]/out.freqs[, 1]
    }       

        
    if(n.seasons > 1) {
        out.seasons <- vector("list", n.seasons)
        names(out.seasons) <- seasonNames
        
        ##for each season, determine frequencies
        out.freqs <- matrix(data = NA, ncol = 6, nrow = n.seasons)
        colnames(out.freqs) <- c("sampled", "detected", "colonized",
                                 "extinct", "static", "common")
        rownames(out.freqs) <- paste("Season-", 1:n.seasons, sep = "")
        
        ##sequence of visits
        for(i in 1:n.seasons) {
            ySeason <- yMat.seasons[[i]]
            censored <- censoredDist.seasons[[i]]
            
            uncensoredObs <- matrix(NA, ncol = n.visits.season,
                                    nrow = nsites)
            
            for(j in 1:ncol(ySeason)){
                ##observations
                uncensoredObs[, j] <- ySeason[, j] < censored[, j]
            }
            ##determine proportion of sites with at least 1 detection
            det.sum <- rowSums(uncensoredObs, na.rm = TRUE)
            
            ##check sites with observed detections and deal with NA's
            sum.rows <- rowSums(uncensoredObs, na.rm = TRUE)
            is.na(sum.rows) <- rowSums(is.na(uncensoredObs)) == ncol(ySeason)
                
            ##number of sites sampled
            out.freqs[i, 1] <- sum(!is.na(sum.rows))
            out.freqs[i, 2] <- sum(det.sum)
            
            
            ##sites without detections
            none <- which(sum.rows == 0)
            ##sites with at least one detection
            some <- which(sum.rows != 0) 
            out.seasons[[i]] <- list("none" = none, "some" = some)
        }
  
        ##populate out.freqs with freqs of extinctions and colonizations
        for(j in 2:n.seasons) {
            none1 <- out.seasons[[j-1]]$none
            some1 <- out.seasons[[j-1]]$some
            none2 <- out.seasons[[j]]$none
            some2 <- out.seasons[[j]]$some

            ##add check for seasons without sampling or previous season without sampling
            if(y.seasonsNA[j] || y.seasonsNA[j-1]) {
                if(y.seasonsNA[j]) {
                    out.freqs[j, 2:6] <- NA
                }
                if(y.seasonsNA[j-1]) {
                    out.freqs[j, 3:6] <- NA
                }
            } else {
                ##colonizations
                out.freqs[j, 3] <- sum(duplicated(c(some2, none1)))
                ##extinctions
                out.freqs[j, 4] <- sum(duplicated(c(some1, none2)))
                ##no change
                out.freqs[j, 5] <- sum(duplicated(c(some1, some2))) + sum(duplicated(c(none1, none2)))
                ##sites both sampled in t and t-1
                year1 <- c(none1, some1)
                year2 <- c(none2, some2)
                out.freqs[j, 6] <- sum(duplicated(c(year1, year2)))
            }
        }
        
        ##create a matrix with proportion of sites with colonizations
        ##and extinctions based on raw data
        out.props <- matrix(NA, nrow = nrow(out.freqs), ncol = 4)
        colnames(out.props) <- c("naive.occ", "naive.colonization", "naive.extinction", "naive.static")
        rownames(out.props) <- rownames(out.freqs)
        
        for(k in 1:n.seasons) {
            ##proportion of sites with detections
            out.props[k, 1] <- out.freqs[k, 2]/out.freqs[k, 1]
            ##add check for seasons without sampling
            if(y.seasonsNA[k]) {
                out.props[k, 2:4] <- NA
            } else {
                ##proportion colonized
                out.props[k, 2] <- out.freqs[k, 3]/out.freqs[k, 6]
                ##proportion extinct
                out.props[k, 3] <- out.freqs[k, 4]/out.freqs[k, 6]
                ##proportion static
                out.props[k, 4] <- out.freqs[k, 5]/out.freqs[k, 6]
            }
        }
    }

    ##reset to original values
    if(any(plot.time || plot.seasons)) {
        on.exit(par(oldpar))
    }
    
    out.det <- list("time.table.full" = time.table.full,
                    "time.table.seasons" = time.table.seasons,
                    "out.freqs" = out.freqs, "out.props" = out.props,
                    "n.seasons" = n.seasons,
                    "n.visits.season" = n.visits.season,
                    "missing.seasons" = y.seasonsNA)
    class(out.det) <- "detTime"
    return(out.det)
}
    


##for occuTTD
detTime.unmarkedFitOccuTTD <- function(object, plot.time = TRUE, plot.seasons = FALSE,
                                       cex.axis = 1, cex.lab = 1, cex.main = 1, ...) {

    ##extract data
    yMat <- object@data@y
    nsites <- nrow(yMat)
    n.seasons <- object@data@numPrimary
    nvisits <- ncol(yMat)
    ##visits per season
    n.visits.season <- nvisits/n.seasons
    n.seasons.adj <- n.seasons #total number of plots fixed to 11 or 12, depending on plots requested
    seasonNames <- paste("season", 1:n.seasons, sep = "")
    surveyLength <- object@data@surveyLength
    
        if(plot.seasons && n.seasons == 1) {
        warning("\nCannot plot data across seasons with only 1 season of data: reset to plot.seasons = FALSE\n")
        plot.seasons <- FALSE
    }
    
    if(plot.time && !plot.seasons) {
        nRows <- 1
        nCols <- 1
        ##reset graphics parameters and save in object
        oldpar <- par(mfrow = c(nRows, nCols))
    }

    
    if(!plot.time && plot.seasons) {
        ##determine arrangement of plots in matrix
        if(plot.seasons && n.seasons >= 12) {
            n.seasons.adj <- 12
            warning("\nOnly first 12 seasons are plotted\n")
        }
    
        if(plot.seasons && n.seasons.adj <= 12) {

            ##if n.seasons < 12
            ##if 12, 11, 10 <- 4 x 3
            ##if 9, 8, 7 <- 3 x 3
            ##if 6, 5 <- 3 x 2
            ##if 4 <- 2 x 2
            ##if 3 <- 3 x 1
            ##if 2 <- 2 x 1
            
            if(n.seasons.adj >= 10) {
                nRows <- 4
                nCols <- 3
            } else {
            
                if(n.seasons.adj >= 7) {
                    nRows <- 3
                    nCols <- 3
                } else {
                    
                    if(n.seasons.adj >= 5) {
                        nRows <- 3
                        nCols <- 2
                    } else {
                        if(n.seasons.adj == 4) {
                            nRows <- 2
                            nCols <- 2
                        } else {
                            if(n.seasons.adj == 3) {
                                nRows <- 3
                                nCols <- 1
                            } else {
                                nRows <- 2
                                nCols <- 1
                            }
                        }
                    }
                }
            }
            ##reset graphics parameters and save in object
            oldpar <- par(mfrow = c(nRows, nCols))
        }
    }



    ##if both plots for seasons and combined are requested
    if(plot.time && plot.seasons) {
        
        if(plot.seasons && n.seasons >= 12) {
            n.seasons.adj <- 11
            warning("\nOnly first 11 seasons are plotted\n")
        }

        if(plot.seasons && n.seasons.adj <= 11) {
            
            if(n.seasons.adj >= 9) {
                nRows <- 4
                nCols <- 3
            } else {
                
                if(n.seasons.adj >= 6) {
                    nRows <- 3
                    nCols <- 3
                } else {
                    
                    if(n.seasons.adj >= 4) {
                        nRows <- 3
                        nCols <- 2
                    } else {
                        if(n.seasons.adj == 3) {
                            nRows <- 2
                            nCols <- 2
                        } else {
                            if(n.seasons.adj == 2) {
                                nRows <- 3
                                nCols <- 1
                            }
                        }
                    }
                }
            }
        }
        
        ##reset graphics parameters and save in object
        oldpar <- par(mfrow = c(nRows, nCols))
    }

        
            
    ##combine all seasons
    ##Censoring distance
    censoredDist.full <- surveyLength
    uniqueDist.full <- unique(as.vector(censoredDist.full))
        
    ##determine data that were censored
    uncensoredIndex.full <- yMat < censoredDist.full
    uncensoredData.full <- yMat[uncensoredIndex.full]
    ncensored.full <- nsites - sum(uncensoredIndex.full)

    
    if(plot.time) {
    
        ##check that maximum times are the same for all sites
        if(n.seasons == 1) {
            
            if(length(uniqueDist.full) == 1) {
                main.title <- paste("Distribution of time to detection (survey length: ", uniqueDist.full, " min.)",
                                    sep = "")
            } else {
                main.title  <- paste("Distribution of time to detection (survey length: ", min(uniqueDist.full), "-",
                                     max(uniqueDist.full), " min.)",
                                     sep = "")
            }
        }
        
        if(n.seasons > 1) {
            
            if(length(uniqueDist.full) == 1) {
                main.title <- paste("Distribution of time to detection (", n.seasons, " seasons, survey length: ",
                                    uniqueDist.full, " min.)",
                                    sep = "")
            } else {
                main.title  <- paste("Distribution of time to detection (", n.seasons, " seasons, survey length: ",
                                     min(uniqueDist.full), "-",
                                     max(uniqueDist.full), " min.)",
                                     sep = "")
            }
        }
        
        hist(uncensoredData.full, xlim = c(0, uniqueDist.full),
             xlab = "Time to detection (min.)",
             main = main.title,
             cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
    }


    ##quantiles for entire data
    time.table.full <- quantile(uncensoredData.full, na.rm = TRUE)

    ##store data for each season

    ##columns for each season
    col.seasons <- seq(1, nvisits, by = n.visits.season)

    ##list to store raw data
    yMat.seasons <- vector(mode = "list", length = n.seasons)
    names(yMat.seasons) <- seasonNames
    minusOne <- n.visits.season - 1

    ##list to store quantiles excluding censored times
    time.table.seasons <- vector("list", n.seasons)
    names(time.table.seasons) <- seasonNames

    ##list of uncensored observations
    uncensoredData.seasons <- vector("list", n.seasons)
    names(uncensoredData.seasons) <- seasonNames
    
    ##vector of censored observations
    censored.seasons <- vector("numeric", n.seasons)
    names(censored.seasons) <- seasonNames

    ##list of unique values of maximum effort
    uniqueDist.seasons <- vector("list", n.seasons)
    names(uniqueDist.seasons) <- seasonNames

    ##list of maximum effort
    censoredDist.seasons <- vector("list", n.seasons)
    names(censoredDist.seasons) <- seasonNames
    
    ##iterate over each season
    for(i in 1:n.seasons) {
        ##extract yMat for each season
        yMat1 <- yMat[, col.seasons[i]:(col.seasons[i]+minusOne), drop = FALSE]
        yMat.seasons[[i]] <- yMat1
            
        ##Censoring values
        censoredDist <- surveyLength[, col.seasons[i]:(col.seasons[i]+minusOne), drop = FALSE]
        censoredDist.seasons[[i]] <- censoredDist
        uniqueDist.seasons[[i]] <- unique(as.vector(censoredDist))
        
        ##determine data that were censored
        uncensoredIndex <- yMat1 < censoredDist
        uncensoredData <- yMat1[uncensoredIndex]
        uncensoredData.seasons[[i]] <- uncensoredData
        ncensored <- nsites - sum(uncensoredIndex)

        ##summarize times per season
        time.quantiles <- quantile(uncensoredData, na.rm = TRUE)
        time.table.seasons[[i]] <- time.quantiles
        censored.seasons[i] <- ncensored

    }
        
    ##check if any seasons were not sampled
    y.seasonsNA <- sapply(yMat.seasons, FUN = function(i) all(is.na(i)))

    if(plot.seasons) {

        for(i in 1:n.seasons) {
            
            if(length(uniqueDist.seasons[[i]]) == 1) {
                main.title <- paste("Distribution of time to detection (season ", i,
                                    ", survey length: ", uniqueDist.seasons[[i]],
                                    " min.)", sep = "")
            } else {
                main.title <- paste("Distribution of time to detection (season ", i,
                                    ", survey length: ", min(uniqueDist.seasons[[i]]), "-",
                                    max(uniqueDist.seasons[[i]]), " min.)",
                                    sep = "")
            }
            
            
            ##create histogram
            ##check for missing season
            if(y.seasonsNA[i]) {next} #skip to next season if current season not sampled
            hist(uncensoredData.seasons[[i]], xlim = c(0, uniqueDist.seasons[[i]]),
                 xlab = "Time to detection (min.)",
                 main = main.title,
                 cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
        }
    }
        
    if(n.seasons == 1) {
        ##for each season, determine frequencies
        out.freqs <- matrix(data = NA, ncol = 2, nrow = n.seasons)
        colnames(out.freqs) <- c("sampled", "detected")
        rownames(out.freqs) <- "Season-1"
        
        ##sequence of visits
        for(i in 1:n.seasons) {
            ySeason <- yMat.seasons[[i]]
            censored <- censoredDist.seasons[[i]]

            uncensoredObs <- matrix(NA, ncol = n.visits.season,
                                    nrow = nsites)
            
            for(j in 1:ncol(ySeason)){
                ##observations
                uncensoredObs[, j] <- ySeason[, j] < censored[, j]
            }
            ##determine proportion of sites with at least 1 detection
            det.sum <- rowSums(uncensoredObs, na.rm = TRUE)
            
            ##check sites with observed detections and deal with NA's
            sum.rows <- rowSums(uncensoredObs, na.rm = TRUE)
            is.na(sum.rows) <- rowSums(is.na(uncensoredObs)) == ncol(ySeason)
                
            ##number of sites sampled
            out.freqs[i, 1] <- sum(!is.na(sum.rows))
            out.freqs[i, 2] <- sum(det.sum)
            
        }
              
    
        ##create a matrix with proportion of sites with colonizations
        ##and extinctions based on raw data
        out.props <- matrix(NA, nrow = nrow(out.freqs), ncol = 1)
        colnames(out.props) <- "naive.occ"
        rownames(out.props) <- rownames(out.freqs)
        out.props[, 1] <- out.freqs[, 2]/out.freqs[, 1]
    }       

        
    if(n.seasons > 1) {
        out.seasons <- vector("list", n.seasons)
        names(out.seasons) <- seasonNames
        
        ##for each season, determine frequencies
        out.freqs <- matrix(data = NA, ncol = 6, nrow = n.seasons)
        colnames(out.freqs) <- c("sampled", "detected", "colonized",
                                 "extinct", "static", "common")
        rownames(out.freqs) <- paste("Season-", 1:n.seasons, sep = "")
        
        ##sequence of visits
        for(i in 1:n.seasons) {
            ySeason <- yMat.seasons[[i]]
            censored <- censoredDist.seasons[[i]]
            
            uncensoredObs <- matrix(NA, ncol = n.visits.season,
                                    nrow = nsites)
            
            for(j in 1:ncol(ySeason)){
                ##observations
                uncensoredObs[, j] <- ySeason[, j] < censored[, j]
            }
            ##determine proportion of sites with at least 1 detection
            det.sum <- rowSums(uncensoredObs, na.rm = TRUE)
            
            ##check sites with observed detections and deal with NA's
            sum.rows <- rowSums(uncensoredObs, na.rm = TRUE)
            is.na(sum.rows) <- rowSums(is.na(uncensoredObs)) == ncol(ySeason)
            
            ##number of sites sampled
            out.freqs[i, 1] <- sum(!is.na(sum.rows))
            out.freqs[i, 2] <- sum(det.sum)
            
            
            ##sites without detections
            none <- which(sum.rows == 0)
            ##sites with at least one detection
            some <- which(sum.rows != 0) 
            out.seasons[[i]] <- list("none" = none, "some" = some)
        }
  
        ##populate out.freqs with freqs of extinctions and colonizations
        for(j in 2:n.seasons) {
            none1 <- out.seasons[[j-1]]$none
            some1 <- out.seasons[[j-1]]$some
            none2 <- out.seasons[[j]]$none
            some2 <- out.seasons[[j]]$some

            ##add check for seasons without sampling or previous season without sampling
            if(y.seasonsNA[j] || y.seasonsNA[j-1]) {
                if(y.seasonsNA[j]) {
                    out.freqs[j, 2:6] <- NA
                }
                if(y.seasonsNA[j-1]) {
                    out.freqs[j, 3:6] <- NA
                }
            } else {

                ##colonizations
                out.freqs[j, 3] <- sum(duplicated(c(some2, none1)))
                ##extinctions
                out.freqs[j, 4] <- sum(duplicated(c(some1, none2)))
                ##no change
                out.freqs[j, 5] <- sum(duplicated(c(some1, some2))) + sum(duplicated(c(none1, none2)))
                ##sites both sampled in t and t-1
                year1 <- c(none1, some1)
                year2 <- c(none2, some2)
                out.freqs[j, 6] <- sum(duplicated(c(year1, year2)))
            }
        }
                
        ##create a matrix with proportion of sites with colonizations
        ##and extinctions based on raw data
        out.props <- matrix(NA, nrow = nrow(out.freqs), ncol = 4)
        colnames(out.props) <- c("naive.occ", "naive.colonization", "naive.extinction", "naive.static")
        rownames(out.props) <- rownames(out.freqs)

        for(k in 1:n.seasons) {
            ##proportion of sites with detections
            out.props[k, 1] <- out.freqs[k, 2]/out.freqs[k, 1]
            ##add check for seasons without sampling
            if(y.seasonsNA[k]) {
                out.props[k, 2:4] <- NA
            } else {
                ##proportion colonized
                out.props[k, 2] <- out.freqs[k, 3]/out.freqs[k, 6]
                ##proportion extinct
                out.props[k, 3] <- out.freqs[k, 4]/out.freqs[k, 6]
                ##proportion static
                out.props[k, 4] <- out.freqs[k, 5]/out.freqs[k, 6]
            }
        }
    }

    ##reset to original values
    if(any(plot.time || plot.seasons)) {
        on.exit(par(oldpar))
    }

    out.det <- list("time.table.full" = time.table.full,
                    "time.table.seasons" = time.table.seasons,
                    "out.freqs" = out.freqs, "out.props" = out.props,
                    "n.seasons" = n.seasons,
                    "n.visits.season" = n.visits.season,
                    "missing.seasons" = y.seasonsNA)
    class(out.det) <- "detTime"
    return(out.det)
}



##print method
print.detTime <- function(x, digits = 2, ...) {
    if(identical(x$n.seasons, 1)) {

        cat("\nSummary of time to detection:\n")
        time.mat <- matrix(x$time.table.full, nrow = 1)
        colnames(time.mat) <- names(x$time.table.full)
        rownames(time.mat) <- "Times"
        print(round(time.mat, digits = digits))
        
    cat("\nProportion of sites with at least one detection:\n", round(x$out.props[, "naive.occ"], digits), "\n\n")
    
    cat("Frequencies of sites with detections:\n")
    ##add matrix of frequencies
    print(x$out.freqs)

  } else {
    cat("\nSummary of time to detection (", x$n.seasons, " seasons combined): \n", sep ="")
    time.mat <- matrix(x$time.table.full, nrow = 1)
    colnames(time.mat) <- names(x$time.table.full)
    rownames(time.mat) <- "Time"
    print(round(time.mat, digits = digits))

    ##if some seasons have not been sampled
    if(any(x$missing.seasons)) {
        if(sum(x$missing.seasons) == 1) {
            cat("\nNote: season", which(x$missing.seasons), "was not sampled\n")
        } else {
            cat("\nNote: seasons",
                paste(which(x$missing.seasons), sep = ", "),
                "were not sampled\n")
        }
    }            
    cat("\nSeason-specific time to detection: \n")
    cat("\n")
    for(i in 1:x$n.seasons) {
        if(!x$missing.seasons[i]) {
            cat("Season", i, "\n")
        } else {
            cat("Season", i, "(no sites sampled)", "\n")
        }
        
        temp.tab <- x$time.table.seasons[[i]]
        out.mat <- matrix(temp.tab, nrow = 1)
        colnames(out.mat) <- names(temp.tab)
        rownames(out.mat) <- "Time"
        print(round(out.mat, digits = digits))
        cat("--------\n\n")
    }
    
    cat("Frequencies of sites with detections, extinctions, and colonizations:\n")
    ##add matrix of frequencies
    print(x$out.freqs)
  }
}



