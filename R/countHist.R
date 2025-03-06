##summarize detection histories and count data
countHist <- function(object, plot.freq = TRUE,
                      cex.axis = 1, cex.lab = 1, cex.main = 1, ...){
  UseMethod("countHist", object)
}



countHist.default <- function(object, plot.freq = TRUE,
                              cex.axis = 1, cex.lab = 1, cex.main = 1, ...){
  stop("\nFunction not yet defined for this object class\n")
}



##for unmarkedFramePCount 
countHist.unmarkedFramePCount <- function(object, plot.freq = TRUE,
                                          cex.axis = 1, cex.lab = 1, cex.main = 1, ...) {

    ##extract data
    yMat <- object@y
    nsites <- nrow(yMat)
    n.seasons <- 1
    nvisits <- ncol(yMat)
    ##visits per season
    n.visits.season <- nvisits/n.seasons
    seasonNames <- paste("season", 1:n.seasons, sep = "")
    
    ##collapse yMat into a single vector
    yVec <- as.vector(yMat)
  
    ##summarize detection histories
    if(plot.freq) {

        ##create histogram
        barplot(table(yVec), ylab = "Frequency", xlab = "Counts of individuals",
                main = "Distribution of raw counts",
                cex.axis = cex.axis, cex.names = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
    }

    ##raw counts
    count.table.full <- table(yVec, exclude = NULL, deparse.level = 0)
    count.table.seasons <- list(count.table.full)
    names(count.table.seasons) <- seasonNames

    ##summarize count histories
    hist.full <- apply(X = yMat, MARGIN = 1, FUN = function(i) paste(i, collapse = "|"))
    hist.table.full <- table(hist.full, deparse.level = 0)

    ##for each season, determine frequencies
    hist.table.seasons <- vector(mode = "list", length = n.seasons)
    names(hist.table.seasons) <- seasonNames

    out.freqs <- matrix(data = NA, ncol = 2, nrow = n.seasons)
    colnames(out.freqs) <- c("sampled", "detected")
    rownames(out.freqs) <- "Season-1"

    hist.table.seasons[[1]] <- hist.table.full
    
    ##determine proportion of sites with at least 1 detection
    det.sum <- apply(X = yMat, MARGIN = 1, FUN = function(i) ifelse(sum(i, na.rm = TRUE) > 0, 1, 0))

    ##check sites with observed detections and deal with NA's
    sum.rows <- rowSums(yMat, na.rm = TRUE)
    is.na(sum.rows) <- rowSums(is.na(yMat)) == ncol(yMat)
    
    ##number of sites sampled
    out.freqs[1, 1] <- sum(!is.na(sum.rows))
    ##number of sites with at least 1 detection
    out.freqs[1, 2] <- sum(det.sum)

    ##create a matrix with proportion of sites with colonizations
    ##and extinctions based on raw data
    out.props <- matrix(NA, nrow = nrow(out.freqs), ncol = 1)
    colnames(out.props) <- "naive.occ"
    rownames(out.props) <- rownames(out.freqs)
    out.props[, 1] <- out.freqs[, 2]/out.freqs[, 1]
    
    out.count <- list("count.table.full" = count.table.full,
                      "count.table.seasons" = count.table.seasons,
                      "hist.table.full" = hist.table.full,
                      "hist.table.seasons" = hist.table.seasons,
                      "out.freqs" = out.freqs, "out.props" = out.props,
                      "n.seasons" = n.seasons,
                      "n.visits.season" = n.visits.season,
                      "missing.seasons" = FALSE)
    class(out.count) <- "countHist"
    return(out.count)
}



##for unmarkedFitPCount 
countHist.unmarkedFitPCount <- function(object, plot.freq = TRUE,
                                        cex.axis = 1, cex.lab = 1, cex.main = 1, ...) {

    ##extract data
    yMat <- object@data@y
    nsites <- nrow(yMat)
    n.seasons <- 1
    nvisits <- ncol(yMat)
    ##visits per season
    n.visits.season <- nvisits/n.seasons
    seasonNames <- paste("season", 1:n.seasons, sep = "")

    ##collapse yMat into a single vector
    yVec <- as.vector(yMat)
  
    ##summarize detection histories
    if(plot.freq) {

        ##create histogram
        barplot(table(yVec), ylab = "Frequency", xlab = "Counts of individuals",
                main = "Distribution of raw counts",
                cex.axis = cex.axis, cex.names = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
    }

    ##raw counts
    count.table.full <- table(yVec, exclude = NULL, deparse.level = 0)
    count.table.seasons <- list(count.table.full)
    names(count.table.seasons) <- seasonNames
    
    ##count histories
    hist.full <- apply(X = yMat, MARGIN = 1, FUN = function(i) paste(i, collapse = "|"))
    hist.table.full <- table(hist.full, deparse.level = 0)

    ##for each season, determine frequencies
    hist.table.seasons <- vector(mode = "list", length = n.seasons)
    names(hist.table.seasons) <- seasonNames
    
    out.freqs <- matrix(data = NA, ncol = 2, nrow = n.seasons)
    colnames(out.freqs) <- c("sampled", "detected")
    rownames(out.freqs) <- "Season-1"

    hist.table.seasons[[1]] <- hist.table.full
    
    ##determine proportion of sites with at least 1 detection
    det.sum <- apply(X = yMat, MARGIN = 1, FUN = function(i) ifelse(sum(i, na.rm = TRUE) > 0, 1, 0))

    ##check sites with observed detections and deal with NA's
    sum.rows <- rowSums(yMat, na.rm = TRUE)
    is.na(sum.rows) <- rowSums(is.na(yMat)) == ncol(yMat)
    
    ##number of sites sampled
    out.freqs[1, 1] <- sum(!is.na(sum.rows))
    ##number of sites with at least 1 detection
    out.freqs[1, 2] <- sum(det.sum)

    ##create a matrix with proportion of sites with colonizations
    ##and extinctions based on raw data
    out.props <- matrix(NA, nrow = nrow(out.freqs), ncol = 1)
    colnames(out.props) <- "naive.occ"
    rownames(out.props) <- rownames(out.freqs)
    out.props[, 1] <- out.freqs[, 2]/out.freqs[, 1]

    out.count <- list("count.table.full" = count.table.full,
                      "count.table.seasons" = count.table.seasons,
                      "hist.table.full" = hist.table.full,
                      "hist.table.seasons" = hist.table.seasons,
                      "out.freqs" = out.freqs, "out.props" = out.props,
                      "n.seasons" = n.seasons,
                      "n.visits.season" = n.visits.season,
                      "missing.seasons" = FALSE)
    class(out.count) <- "countHist"
    return(out.count)
}



##for unmarkedFrameGPC 
countHist.unmarkedFrameGPC <- function(object, plot.freq = TRUE,
                                       cex.axis = 1, cex.lab = 1, cex.main = 1, ...) {

    ##extract data
    yMat <- object@y
    nsites <- nrow(yMat)
    n.seasons <- 1
    nvisits <- ncol(yMat)
    ##visits per season
    n.visits.season <- nvisits/n.seasons
    seasonNames <- paste("season", 1:n.seasons, sep = "")

    
    ##collapse yMat into a single vector
    yVec <- as.vector(yMat)
  
    ##summarize detection histories
    if(plot.freq) {

        ##create histogram
        barplot(table(yVec), ylab = "Frequency", xlab = "Counts of individuals",
                main = "Distribution of raw counts",
                cex.axis = cex.axis, cex.names = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
    }

    ##raw counts
    count.table.full <- table(yVec, exclude = NULL, deparse.level = 0)
    count.table.seasons <- list(count.table.full)
    names(count.table.seasons) <- seasonNames

    ##summarize count histories
    hist.full <- apply(X = yMat, MARGIN = 1, FUN = function(i) paste(i, collapse = "|"))
    hist.table.full <- table(hist.full, deparse.level = 0)

    ##for each season, determine frequencies
    hist.table.seasons <- vector(mode = "list", length = n.seasons)
    names(hist.table.seasons) <- seasonNames
    
    out.freqs <- matrix(data = NA, ncol = 2, nrow = n.seasons)
    colnames(out.freqs) <- c("sampled", "detected")
    rownames(out.freqs) <- "Season-1"

    hist.table.seasons[[1]] <- hist.table.full
    
    ##determine proportion of sites with at least 1 detection
    det.sum <- apply(X = yMat, MARGIN = 1, FUN = function(i) ifelse(sum(i, na.rm = TRUE) > 0, 1, 0))

    ##check sites with observed detections and deal with NA's
    sum.rows <- rowSums(yMat, na.rm = TRUE)
    is.na(sum.rows) <- rowSums(is.na(yMat)) == ncol(yMat)
    
    ##number of sites sampled
    out.freqs[1, 1] <- sum(!is.na(sum.rows))
    ##number of sites with at least 1 detection
    out.freqs[1, 2] <- sum(det.sum)
    
    ##create a matrix with proportion of sites with colonizations
    ##and extinctions based on raw data
    out.props <- matrix(NA, nrow = nrow(out.freqs), ncol = 1)
    colnames(out.props) <- "naive.occ"
    rownames(out.props) <- rownames(out.freqs)
    out.props[, 1] <- out.freqs[, 2]/out.freqs[, 1]
    
    out.count <- list("count.table.full" = count.table.full,
                      "count.table.seasons" = count.table.seasons,
                      "hist.table.full" = hist.table.full,
                      "hist.table.seasons" = hist.table.seasons,
                      "out.freqs" = out.freqs, "out.props" = out.props,
                      "n.seasons" = n.seasons,
                      "n.visits.season" = n.visits.season,
                      "missing.seasons" = FALSE)
    class(out.count) <- "countHist"
    return(out.count)
}



##for unmarkedFitGPC 
countHist.unmarkedFitGPC <- function(object, plot.freq = TRUE,
                                     cex.axis = 1, cex.lab = 1, cex.main = 1, ...) {

    ##extract data
    yMat <- object@data@y
    nsites <- nrow(yMat)
    n.seasons <- 1
    nvisits <- ncol(yMat)
    ##visits per season
    n.visits.season <- nvisits/n.seasons
    seasonNames <- paste("season", 1:n.seasons, sep = "")

    
    ##collapse yMat into a single vector
    yVec <- as.vector(yMat)
    
    ##summarize detection histories
    if(plot.freq) {

        ##create histogram
        barplot(table(yVec), ylab = "Frequency", xlab = "Counts of individuals",
                main = "Distribution of raw counts",
                cex.axis = cex.axis, cex.names = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
    }

    ##raw counts
    count.table.full <- table(yVec, exclude = NULL, deparse.level = 0)
    count.table.seasons <- list(count.table.full)
    names(count.table.seasons) <- seasonNames
    
    ##count histories
    hist.full <- apply(X = yMat, MARGIN = 1, FUN = function(i) paste(i, collapse = "|"))
    hist.table.full <- table(hist.full, deparse.level = 0)

    ##for each season, determine frequencies
    hist.table.seasons <- vector(mode = "list", length = n.seasons)
    names(hist.table.seasons) <- seasonNames
        
    out.freqs <- matrix(data = NA, ncol = 2, nrow = n.seasons)
    colnames(out.freqs) <- c("sampled", "detected")
    rownames(out.freqs) <- "Season-1"

    hist.table.seasons[[1]] <- hist.table.full
    
    ##determine proportion of sites with at least 1 detection
    det.sum <- apply(X = yMat, MARGIN = 1, FUN = function(i) ifelse(sum(i, na.rm = TRUE) > 0, 1, 0))

    ##check sites with observed detections and deal with NA's
    sum.rows <- rowSums(yMat, na.rm = TRUE)
    is.na(sum.rows) <- rowSums(is.na(yMat)) == ncol(yMat)
    
    ##number of sites sampled
    out.freqs[1, 1] <- sum(!is.na(sum.rows))
    ##number of sites with at least 1 detection
    out.freqs[1, 2] <- sum(det.sum)

    ##create a matrix with proportion of sites with colonizations
    ##and extinctions based on raw data
    out.props <- matrix(NA, nrow = nrow(out.freqs), ncol = 1)
    colnames(out.props) <- "naive.occ"
    rownames(out.props) <- rownames(out.freqs)
    out.props[, 1] <- out.freqs[, 2]/out.freqs[, 1]

    out.count <- list("count.table.full" = count.table.full,
                      "count.table.seasons" = count.table.seasons,
                      "hist.table.full" = hist.table.full,
                      "hist.table.seasons" = hist.table.seasons,
                      "out.freqs" = out.freqs, "out.props" = out.props,
                      "n.seasons" = n.seasons,
                      "n.visits.season" = n.visits.season,
                      "missing.seasons" = FALSE)
    class(out.count) <- "countHist"
    return(out.count)
}



##for unmarkedFrameMPois 
countHist.unmarkedFrameMPois <- function(object, plot.freq = TRUE,
                                         cex.axis = 1, cex.lab = 1, cex.main = 1, ...) {

    ##extract data
    yMat <- object@y
    nsites <- nrow(yMat)
    n.seasons <- 1
    nvisits <- ncol(yMat)
    ##visits per season
    n.visits.season <- nvisits/n.seasons
    seasonNames <- paste("season", 1:n.seasons, sep = "")
    
    ##collapse yMat into a single vector
    yVec <- as.vector(yMat)
    
    ##summarize detection histories
    if(plot.freq) {

        ##create histogram
        barplot(table(yVec), ylab = "Frequency", xlab = "Counts of individuals",
                main = "Distribution of raw counts",
                cex.axis = cex.axis, cex.names = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
    }

    ##raw counts
    count.table.full <- table(yVec, exclude = NULL, deparse.level = 0)
    count.table.seasons <- list(count.table.full)
    names(count.table.seasons) <- seasonNames
    
    ##summarize count histories
    hist.full <- apply(X = yMat, MARGIN = 1, FUN = function(i) paste(i, collapse = "|"))
    hist.table.full <- table(hist.full, deparse.level = 0)

    ##for each season, determine frequencies
    hist.table.seasons <- vector(mode = "list", length = n.seasons)
    names(hist.table.seasons) <- seasonNames
    
    out.freqs <- matrix(data = NA, ncol = 2, nrow = n.seasons)
    colnames(out.freqs) <- c("sampled", "detected")
    rownames(out.freqs) <- "Season-1"

    hist.table.seasons[[1]] <- hist.table.full
    
    ##determine proportion of sites with at least 1 detection
    det.sum <- apply(X = yMat, MARGIN = 1, FUN = function(i) ifelse(sum(i, na.rm = TRUE) > 0, 1, 0))

    ##check sites with observed detections and deal with NA's
    sum.rows <- rowSums(yMat, na.rm = TRUE)
    is.na(sum.rows) <- rowSums(is.na(yMat)) == ncol(yMat)
    
    ##number of sites sampled
    out.freqs[1, 1] <- sum(!is.na(sum.rows))
    ##number of sites with at least 1 detection
    out.freqs[1, 2] <- sum(det.sum)

    ##create a matrix with proportion of sites with colonizations
    ##and extinctions based on raw data
    out.props <- matrix(NA, nrow = nrow(out.freqs), ncol = 1)
    colnames(out.props) <- "naive.occ"
    rownames(out.props) <- rownames(out.freqs)
    out.props[, 1] <- out.freqs[, 2]/out.freqs[, 1]

    out.count <- list("count.table.full" = count.table.full,
                      "count.table.seasons" = count.table.seasons,
                      "hist.table.full" = hist.table.full,
                      "hist.table.seasons" = hist.table.seasons,
                      "out.freqs" = out.freqs, "out.props" = out.props,
                      "n.seasons" = n.seasons,
                      "n.visits.season" = n.visits.season,
                      "missing.seasons" = FALSE)
    class(out.count) <- "countHist"
    return(out.count)
}



##for unmarkedFitMPois 
countHist.unmarkedFitMPois <- function(object, plot.freq = TRUE,
                                       cex.axis = 1, cex.lab = 1, cex.main = 1, ...) {

    ##extract data
    yMat <- object@data@y
    nsites <- nrow(yMat)
    n.seasons <- 1
    nvisits <- ncol(yMat)
    ##visits per season
    n.visits.season <- nvisits/n.seasons
    seasonNames <- paste("season", 1:n.seasons, sep = "")
    
    ##collapse yMat into a single vector
    yVec <- as.vector(yMat)
    
    ##summarize detection histories
    if(plot.freq) {

        ##create histogram
        barplot(table(yVec), ylab = "Frequency", xlab = "Counts of individuals",
                main = "Distribution of raw counts",
                cex.axis = cex.axis, cex.names = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
    }
    
    ##raw counts
    count.table.full <- table(yVec, exclude = NULL, deparse.level = 0)
    count.table.seasons <- list(count.table.full)
    names(count.table.seasons) <- seasonNames
    
    ##summarize count histories
    hist.full <- apply(X = yMat, MARGIN = 1, FUN = function(i) paste(i, collapse = "|"))
    hist.table.full <- table(hist.full, deparse.level = 0)

    ##for each season, determine frequencies
    hist.table.seasons <- vector(mode = "list", length = n.seasons)
    names(hist.table.seasons) <- seasonNames
    
    out.freqs <- matrix(data = NA, ncol = 2, nrow = n.seasons)
    colnames(out.freqs) <- c("sampled", "detected")
    rownames(out.freqs) <- "Season-1"

    hist.table.seasons[[1]] <- hist.table.full
    
    ##determine proportion of sites with at least 1 detection
    det.sum <- apply(X = yMat, MARGIN = 1, FUN = function(i) ifelse(sum(i, na.rm = TRUE) > 0, 1, 0))

    ##check sites with observed detections and deal with NA's
    sum.rows <- rowSums(yMat, na.rm = TRUE)
    is.na(sum.rows) <- rowSums(is.na(yMat)) == ncol(yMat)
    
    ##number of sites sampled
    out.freqs[1, 1] <- sum(!is.na(sum.rows))
    ##number of sites with at least 1 detection
    out.freqs[1, 2] <- sum(det.sum)

    ##create a matrix with proportion of sites with colonizations
    ##and extinctions based on raw data
    out.props <- matrix(NA, nrow = nrow(out.freqs), ncol = 1)
    colnames(out.props) <- "naive.occ"
    rownames(out.props) <- rownames(out.freqs)
    out.props[, 1] <- out.freqs[, 2]/out.freqs[, 1]

    out.count <- list("count.table.full" = count.table.full,
                      "count.table.seasons" = count.table.seasons,
                      "hist.table.full" = hist.table.full,
                      "hist.table.seasons" = hist.table.seasons,
                      "out.freqs" = out.freqs, "out.props" = out.props,
                      "n.seasons" = n.seasons,
                      "n.visits.season" = n.visits.season,
                      "missing.seasons" = FALSE)
    class(out.count) <- "countHist"
    return(out.count)
}



##for unmarkedFramePCO
countHist.unmarkedFramePCO <- function(object, plot.freq = TRUE,
                                       cex.axis = 1, cex.lab = 1, cex.main = 1,
                                       plot.seasons = FALSE, ...) {

    ##extract data
    yMat <- object@y
    nsites <- nrow(yMat)
    n.seasons <- object@numPrimary
    n.seasons.adj <- n.seasons #total number of plots fixed to 11 or 12, depending on plots requested
    nvisits <- ncol(yMat)
    ##visits per season
    n.visits.season <- nvisits/n.seasons
    seasonNames <- paste("season", 1:n.seasons, sep = "")
    
    ##collapse yMat into a single vector
    yVec.full <- as.vector(yMat)

    ##raw counts
    count.table.full <- table(yVec.full, exclude = NULL, deparse.level = 0)
  
    ##summarize count histories
    hist.full <- apply(X = yMat, MARGIN = 1, FUN = function(i) paste(i, collapse = "|"))
    hist.table.full <- table(hist.full, deparse.level = 0)

    ##for each season, determine frequencies
    yVectors <- vector(mode = "list", length = n.seasons)
    out.seasons <- vector(mode = "list", length = n.seasons)
    count.table.seasons <- vector(mode = "list", length = n.seasons)
    names(count.table.seasons) <- seasonNames

    hist.table.seasons <- vector(mode = "list", length = n.seasons)
    names(hist.table.seasons) <- seasonNames

    
    out.freqs <- matrix(data = NA, ncol = 6, nrow = n.seasons)
    colnames(out.freqs) <- c("sampled", "detected", "colonized",
                             "extinct", "static", "common")
    rownames(out.freqs) <- paste("Season-", 1:n.seasons, sep = "")

    ##sequence of visits
    vis.seq <- seq(from = 1, to = nvisits, by = n.visits.season)
    for(i in 1:n.seasons) {
        col.start <- vis.seq[i]
        col.end <- col.start + (n.visits.season - 1)
        ySeason <- yMat[, col.start:col.end]
        ##summarize count histories
        if(is.null(ncol(ySeason))){
            ySeason <- as.matrix(ySeason)
        }
        yVec.season <- as.vector(ySeason)
        yVectors[[i]] <- yVec.season
        
        det.hist <- apply(X = ySeason, MARGIN = 1, FUN = function(i) paste(i, collapse = "|"))
        hist.table.seasons[[i]] <- table(det.hist, deparse.level = 0)
        count.table.seasons[[i]] <- table(yVec.season, exclude = NULL)

        ##determine proportion of sites with at least 1 detection
        det.sum <- apply(X = ySeason, MARGIN = 1, FUN = function(i) ifelse(sum(i, na.rm = TRUE) > 0, 1, 0))

        ##check sites with observed detections and deal with NA's
        sum.rows <- rowSums(ySeason, na.rm = TRUE)
        is.na(sum.rows) <- rowSums(is.na(ySeason)) == ncol(ySeason)
    
        ##number of sites sampled
        out.freqs[i, 1] <- sum(!is.na(sum.rows))
        out.freqs[i, 2] <- sum(det.sum)


        ##sites without detections
        none <- which(sum.rows == 0)
        ##sites with at least one detection
        some <- which(sum.rows != 0) 
        out.seasons[[i]] <- list("none" = none, "some" = some)
    }

    ##check if any seasons were not sampled
    y.seasonsNA <- sapply(yVectors, FUN = function(i) all(is.na(i)))

    
    ##if only season-specific plots are requested
    if(!plot.freq && plot.seasons) {
        ##determine arrangement of plots in matrix
        if(n.seasons >= 12) {
            n.seasons.adj <- 12
            warning("\nOnly first 12 seasons are plotted\n")
        }
        
        if(n.seasons.adj <= 12) {

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
        }
        ##reset graphics parameters and save in object
        oldpar <- par(mfrow = c(nRows, nCols))
    }


    ##if both plots for seasons and combined are requested
    ##summarize detection histories
    if(plot.freq) {

        if(!plot.seasons) {
            nRows <- 1
            nCols <- 1
        }
        
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
        
        ##histogram for data combined across seasons
        barplot(table(yVec.full), ylab = "Frequency", xlab = "Counts of individuals",
                main = paste("Distribution of raw counts (", n.seasons, " seasons combined)", sep = ""),
                cex.axis = cex.axis, cex.names = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
    }

    ##iterate over each season
    if(plot.seasons) {
        for(k in 1:n.seasons.adj) {
            ##check for missing season
            if(y.seasonsNA[k]) {next} #skip to next season if current season not sampled
            ##histogram for data combined across seasons
            barplot(table(yVectors[[k]]), ylab = "Frequency", xlab = "Counts of individuals",
                    main = paste("Distribution of raw counts (season ", k, ")", sep = ""),
                    cex.axis = cex.axis, cex.names = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
        }
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

    ##reset to original values
    if(any(plot.freq || plot.seasons)) {
        on.exit(par(oldpar))
    }
  
    out.count <- list("count.table.full" = count.table.full,
                      "count.table.seasons" = count.table.seasons,
                      "hist.table.full" = hist.table.full,
                      "hist.table.seasons" = hist.table.seasons,
                      "out.freqs" = out.freqs, "out.props" = out.props,
                      "n.seasons" = n.seasons,
                      "n.visits.season" = n.visits.season,
                      "missing.seasons" = y.seasonsNA)
    class(out.count) <- "countHist"
    return(out.count)
}



##for unmarkedFitPCO
countHist.unmarkedFitPCO <- function(object, plot.freq = TRUE, 
                                     cex.axis = 1, cex.lab = 1, cex.main = 1,
                                     plot.seasons = FALSE, ...) {

    ##extract data
    yMat <- object@data@y
    nsites <- nrow(yMat)
    n.seasons <- object@data@numPrimary
    n.seasons.adj <- n.seasons #total number of plots fixed to 11 or 12, depending on plots requested
    nvisits <- ncol(yMat)
    ##visits per season
    n.visits.season <- nvisits/n.seasons
    seasonNames <- paste("season", 1:n.seasons, sep = "")

    ##collapse yMat into a single vector
    yVec.full <- as.vector(yMat)
  
    ##raw counts
    count.table.full <- table(yVec.full, exclude = NULL, deparse.level = 0)
  
    ##summarize count histories
    hist.full <- apply(X = yMat, MARGIN = 1, FUN = function(i) paste(i, collapse = "|"))
    hist.table.full <- table(hist.full, deparse.level = 0)

    ##for each season, determine frequencies
    yVectors <- vector(mode = "list", length = n.seasons)
    out.seasons <- vector(mode = "list", length = n.seasons)
    count.table.seasons <- vector(mode = "list", length = n.seasons)
    names(count.table.seasons) <- seasonNames
    
    hist.table.seasons <- vector(mode = "list", length = n.seasons)
    names(hist.table.seasons) <- seasonNames
    
    out.freqs <- matrix(data = NA, ncol = 6, nrow = n.seasons)
    colnames(out.freqs) <- c("sampled", "detected", "colonized",
                             "extinct", "static", "common")
    rownames(out.freqs) <- paste("Season-", 1:n.seasons, sep = "")

    ##sequence of visits
    vis.seq <- seq(from = 1, to = nvisits, by = n.visits.season)
    for(i in 1:n.seasons) {
        col.start <- vis.seq[i]
        col.end <- col.start + (n.visits.season - 1)
        ySeason <- yMat[, col.start:col.end]
        ##summarize count histories
        if(is.null(ncol(ySeason))){
            ySeason <- as.matrix(ySeason)
        }
        yVec.season <- as.vector(ySeason)
        yVectors[[i]] <- yVec.season
        
        det.hist <- apply(X = ySeason, MARGIN = 1, FUN = function(i) paste(i, collapse = "|"))
        hist.table.seasons[[i]] <- table(det.hist, deparse.level = 0)
        count.table.seasons[[i]] <- table(yVec.season, exclude = NULL)
        
        ##determine proportion of sites with at least 1 detection
        det.sum <- apply(X = ySeason, MARGIN = 1, FUN = function(i) ifelse(sum(i, na.rm = TRUE) > 0, 1, 0))
        
        ##check sites with observed detections and deal with NA's
        sum.rows <- rowSums(ySeason, na.rm = TRUE)
        is.na(sum.rows) <- rowSums(is.na(ySeason)) == ncol(ySeason)
        
        ##number of sites sampled
        out.freqs[i, 1] <- sum(!is.na(sum.rows))
        out.freqs[i, 2] <- sum(det.sum)
        
        
        ##sites without detections
        none <- which(sum.rows == 0)
        ##sites with at least one detection
        some <- which(sum.rows != 0) 
        out.seasons[[i]] <- list("none" = none, "some" = some)
    }
    
    ##check if any seasons were not sampled
    y.seasonsNA <- sapply(yVectors, FUN = function(i) all(is.na(i)))
    

    ##if only season-specific plots are requested
    if(!plot.freq && plot.seasons) {
        ##determine arrangement of plots in matrix
        if(n.seasons >= 12) {
            n.seasons.adj <- 12
            warning("\nOnly first 12 seasons are plotted\n")
        }
    
        if(n.seasons.adj <= 12) {

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
        }
        ##reset graphics parameters and save in object
        oldpar <- par(mfrow = c(nRows, nCols))
    }


    ##if both plots for seasons and combined are requested
    ##summarize detection histories
    if(plot.freq) {

        if(!plot.seasons) {
            nRows <- 1
            nCols <- 1
        }
      
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
      
        ##histogram for data combined across seasons
        barplot(table(yVec.full), ylab = "Frequency", xlab = "Counts of individuals",
                main = paste("Distribution of raw counts (", n.seasons, " seasons combined)", sep = ""),
                cex.axis = cex.axis, cex.names = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
    }

    
    ##iterate over each season
    if(plot.seasons) {
        for(k in 1:n.seasons.adj) {
            ##check for missing season
            if(y.seasonsNA[k]) {next} #skip to next season if current season not sampled
            ##histogram for data combined across seasons
            barplot(table(yVectors[[k]]), ylab = "Frequency", xlab = "Counts of individuals",
                    main = paste("Distribution of raw counts (season ", k, ")", sep = ""),
                    cex.axis = cex.axis, cex.names = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
        }
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
    
    ##reset to original values
    if(any(plot.freq || plot.seasons)) {
        on.exit(par(oldpar))
    }
  
    out.count <- list("count.table.full" = count.table.full,
                      "count.table.seasons" = count.table.seasons,
                      "hist.table.full" = hist.table.full,
                      "hist.table.seasons" = hist.table.seasons,
                      "out.freqs" = out.freqs, "out.props" = out.props,
                      "n.seasons" = n.seasons,
                      "n.visits.season" = n.visits.season,
                      "missing.seasons" = y.seasonsNA)
    class(out.count) <- "countHist"
    return(out.count)
}



##for unmarkedFrameGMM
countHist.unmarkedFrameGMM <- function(object, plot.freq = TRUE,
                                       cex.axis = 1, cex.lab = 1, cex.main = 1,
                                       plot.seasons = FALSE, ...) {

    ##extract data
    yMat <- object@y
    nsites <- nrow(yMat)
    n.seasons <- object@numPrimary
    n.seasons.adj <- n.seasons #total number of plots fixed to 11 or 12, depending on plots requested
    nvisits <- ncol(yMat)
    ##visits per season
    n.visits.season <- nvisits/n.seasons
    seasonNames <- paste("season", 1:n.seasons, sep = "")
    
    ##collapse yMat into a single vector
    yVec.full <- as.vector(yMat)

    ##raw counts
    count.table.full <- table(yVec.full, exclude = NULL, deparse.level = 0)
  
    ##summarize count histories
    hist.full <- apply(X = yMat, MARGIN = 1, FUN = function(i) paste(i, collapse = "|"))
    hist.table.full <- table(hist.full, deparse.level = 0)

    ##for each season, determine frequencies
    yVectors <- vector(mode = "list", length = n.seasons)
    out.seasons <- vector(mode = "list", length = n.seasons)
    count.table.seasons <- vector(mode = "list", length = n.seasons)
    names(count.table.seasons) <- seasonNames
    
    hist.table.seasons <- vector(mode = "list", length = n.seasons)
    names(hist.table.seasons) <- seasonNames
    
    out.freqs <- matrix(data = NA, ncol = 6, nrow = n.seasons)
    colnames(out.freqs) <- c("sampled", "detected", "colonized",
                             "extinct", "static", "common")
    rownames(out.freqs) <- paste("Season-", 1:n.seasons, sep = "")

    ##sequence of visits
    vis.seq <- seq(from = 1, to = nvisits, by = n.visits.season)
    for(i in 1:n.seasons) {
        col.start <- vis.seq[i]
        col.end <- col.start + (n.visits.season - 1)
        ySeason <- yMat[, col.start:col.end]
        ##summarize count histories
        if(is.null(ncol(ySeason))){
            ySeason <- as.matrix(ySeason)
        }
        yVec.season <- as.vector(ySeason)
        yVectors[[i]] <- yVec.season
        
        det.hist <- apply(X = ySeason, MARGIN = 1, FUN = function(i) paste(i, collapse = "|"))
        hist.table.seasons[[i]] <- table(det.hist, deparse.level = 0)
        count.table.seasons[[i]] <- table(yVec.season, exclude = NULL)

        ##determine proportion of sites with at least 1 detection
        det.sum <- apply(X = ySeason, MARGIN = 1, FUN = function(i) ifelse(sum(i, na.rm = TRUE) > 0, 1, 0))
        
        ##check sites with observed detections and deal with NA's
        sum.rows <- rowSums(ySeason, na.rm = TRUE)
        is.na(sum.rows) <- rowSums(is.na(ySeason)) == ncol(ySeason)
    
        ##number of sites sampled
        out.freqs[i, 1] <- sum(!is.na(sum.rows))
        out.freqs[i, 2] <- sum(det.sum)


        ##sites without detections
        none <- which(sum.rows == 0)
        ##sites with at least one detection
        some <- which(sum.rows != 0) 
        out.seasons[[i]] <- list("none" = none, "some" = some)
    }

    ##check if any seasons were not sampled
    y.seasonsNA <- sapply(yVectors, FUN = function(i) all(is.na(i)))

    ##if only season-specific plots are requested
    if(!plot.freq && plot.seasons) {
        ##determine arrangement of plots in matrix
        if(n.seasons >= 12) {
            n.seasons.adj <- 12
            warning("\nOnly first 12 seasons are plotted\n")
        }
    
        if(n.seasons.adj <= 12) {

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
        }
        ##reset graphics parameters and save in object
        oldpar <- par(mfrow = c(nRows, nCols))
    }


    ##if both plots for seasons and combined are requested
    ##summarize detection histories
    if(plot.freq) {

        if(!plot.seasons) {
            nRows <- 1
            nCols <- 1
        }
      
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
    
        ##histogram for data combined across seasons
        barplot(table(yVec.full), ylab = "Frequency", xlab = "Counts of individuals",
                main = paste("Distribution of raw counts (", n.seasons, " seasons combined)", sep = ""),
                cex.axis = cex.axis, cex.names = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
    }


    ##iterate over each season
    if(plot.seasons) {
        for(k in 1:n.seasons.adj) {
            ##check for missing season
            if(y.seasonsNA[k]) {next} #skip to next season if current season not sampled
            ##histogram for data combined across seasons
            barplot(table(yVectors[[k]]), ylab = "Frequency", xlab = "Counts of individuals",
                    main = paste("Distribution of raw counts (season ", k, ")", sep = ""),
                    cex.axis = cex.axis, cex.names = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
        }
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

    ##reset to original values
    if(any(plot.freq || plot.seasons)) {
        on.exit(par(oldpar))
    }
  
    out.count <- list("count.table.full" = count.table.full,
                      "count.table.seasons" = count.table.seasons,
                      "hist.table.full" = hist.table.full,
                      "hist.table.seasons" = hist.table.seasons,
                      "out.freqs" = out.freqs, "out.props" = out.props,
                      "n.seasons" = n.seasons,
                      "n.visits.season" = n.visits.season,
                      "missing.seasons" = y.seasonsNA)
    class(out.count) <- "countHist"
    return(out.count)
}



##for unmarkedFitGMM
countHist.unmarkedFitGMM <- function(object, plot.freq = TRUE, 
                                     cex.axis = 1, cex.lab = 1, cex.main = 1,
                                     plot.seasons = FALSE, ...) {

    ##extract data
    yMat <- object@data@y
    nsites <- nrow(yMat)
    n.seasons <- object@data@numPrimary
    n.seasons.adj <- n.seasons #total number of plots fixed to 11 or 12, depending on plots requested
    nvisits <- ncol(yMat)
    ##visits per season
    n.visits.season <- nvisits/n.seasons
    seasonNames <- paste("season", 1:n.seasons, sep = "")
    
    ##collapse yMat into a single vector
    yVec.full <- as.vector(yMat)
  
    ##raw counts
    count.table.full <- table(yVec.full, exclude = NULL, deparse.level = 0)
  
    ##summarize count histories
    hist.full <- apply(X = yMat, MARGIN = 1, FUN = function(i) paste(i, collapse = "|"))
    hist.table.full <- table(hist.full, deparse.level = 0)

    ##for each season, determine frequencies
    yVectors <- vector(mode = "list", length = n.seasons)
    out.seasons <- vector(mode = "list", length = n.seasons)
    count.table.seasons <- vector(mode = "list", length = n.seasons)
    names(count.table.seasons) <- seasonNames
    
    hist.table.seasons <- vector(mode = "list", length = n.seasons)
    names(hist.table.seasons) <- seasonNames

    out.freqs <- matrix(data = NA, ncol = 6, nrow = n.seasons)
    colnames(out.freqs) <- c("sampled", "detected", "colonized",
                             "extinct", "static", "common")
    rownames(out.freqs) <- paste("Season-", 1:n.seasons, sep = "")

    ##sequence of visits
    vis.seq <- seq(from = 1, to = nvisits, by = n.visits.season)
    for(i in 1:n.seasons) {
        col.start <- vis.seq[i]
        col.end <- col.start + (n.visits.season - 1)
        ySeason <- yMat[, col.start:col.end]
        ##summarize count histories
        if(is.null(ncol(ySeason))){
            ySeason <- as.matrix(ySeason)
        }
        yVec.season <- as.vector(ySeason)
        yVectors[[i]] <- yVec.season
    
        det.hist <- apply(X = ySeason, MARGIN = 1, FUN = function(i) paste(i, collapse = "|"))
        hist.table.seasons[[i]] <- table(det.hist, deparse.level = 0)
        count.table.seasons[[i]] <- table(yVec.season, exclude = NULL)

        ##determine proportion of sites with at least 1 detection
        det.sum <- apply(X = ySeason, MARGIN = 1, FUN = function(i) ifelse(sum(i, na.rm = TRUE) > 0, 1, 0))

        ##check sites with observed detections and deal with NA's
        sum.rows <- rowSums(ySeason, na.rm = TRUE)
        is.na(sum.rows) <- rowSums(is.na(ySeason)) == ncol(ySeason)
    
        ##number of sites sampled
        out.freqs[i, 1] <- sum(!is.na(sum.rows))
        out.freqs[i, 2] <- sum(det.sum)


        ##sites without detections
        none <- which(sum.rows == 0)
        ##sites with at least one detection
        some <- which(sum.rows != 0) 
        out.seasons[[i]] <- list("none" = none, "some" = some)
    }

    ##check if any seasons were not sampled
    y.seasonsNA <- sapply(yVectors, FUN = function(i) all(is.na(i)))

    ##if only season-specific plots are requested
    if(!plot.freq && plot.seasons) {
        ##determine arrangement of plots in matrix
        if(n.seasons >= 12) {
            n.seasons.adj <- 12
            warning("\nOnly first 12 seasons are plotted\n")
        }
    
        if(n.seasons.adj <= 12) {

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
        }
        ##reset graphics parameters and save in object
        oldpar <- par(mfrow = c(nRows, nCols))
  }

    
    ##if both plots for seasons and combined are requested
    ##summarize detection histories
    if(plot.freq) {

        if(!plot.seasons) {
            nRows <- 1
            nCols <- 1
        }
        
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
        
        ##histogram for data combined across seasons
        barplot(table(yVec.full), ylab = "Frequency", xlab = "Counts of individuals",
                main = paste("Distribution of raw counts (", n.seasons, " seasons combined)", sep = ""),
                cex.axis = cex.axis, cex.names = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
    }


    ##iterate over each season
    if(plot.seasons) {
        for(k in 1:n.seasons.adj) {
            ##check for missing season
            if(y.seasonsNA[k]) {next} #skip to next season if current season not sampled
            ##histogram for data combined across seasons
            barplot(table(yVectors[[k]]), ylab = "Frequency", xlab = "Counts of individuals",
                    main = paste("Distribution of raw counts (season ", k, ")", sep = ""),
                    cex.axis = cex.axis, cex.names = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
        }
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
  
    ##reset to original values
    if(any(plot.freq || plot.seasons)) {
        on.exit(par(oldpar))
    }

    out.count <- list("count.table.full" = count.table.full,
                      "count.table.seasons" = count.table.seasons,
                      "hist.table.full" = hist.table.full,
                      "hist.table.seasons" = hist.table.seasons,
                      "out.freqs" = out.freqs, "out.props" = out.props,
                      "n.seasons" = n.seasons,
                      "n.visits.season" = n.visits.season,
                      "missing.seasons" = y.seasonsNA)
    class(out.count) <- "countHist"
    return(out.count)
}



##multmixOpen
countHist.unmarkedFrameMMO <- function(object, plot.freq = TRUE, 
                                       cex.axis = 1, cex.lab = 1, cex.main = 1,
                                       plot.seasons = FALSE, ...) {

    ##extract data
    yMat <- object@y
    nsites <- nrow(yMat)
    n.seasons <- object@numPrimary
    n.seasons.adj <- n.seasons #total number of plots fixed to 11 or 12, depending on plots requested
    nvisits <- ncol(yMat)
    ##visits per season
    n.visits.season <- nvisits/n.seasons
    seasonNames <- paste("season", 1:n.seasons, sep = "")
    
    ##collapse yMat into a single vector
    yVec.full <- as.vector(yMat)
  
    ##raw counts
    count.table.full <- table(yVec.full, exclude = NULL, deparse.level = 0)
  
    ##summarize count histories
    hist.full <- apply(X = yMat, MARGIN = 1, FUN = function(i) paste(i, collapse = "|"))
    hist.table.full <- table(hist.full, deparse.level = 0)

    ##for each season, determine frequencies
    yVectors <- vector(mode = "list", length = n.seasons)
    out.seasons <- vector(mode = "list", length = n.seasons)
    count.table.seasons <- vector(mode = "list", length = n.seasons)
    names(count.table.seasons) <- seasonNames

    hist.table.seasons <- vector(mode = "list", length = n.seasons)
    names(hist.table.seasons) <- seasonNames
    
    out.freqs <- matrix(data = NA, ncol = 6, nrow = n.seasons)
    colnames(out.freqs) <- c("sampled", "detected", "colonized",
                             "extinct", "static", "common")
    rownames(out.freqs) <- paste("Season-", 1:n.seasons, sep = "")

    ##sequence of visits
    vis.seq <- seq(from = 1, to = nvisits, by = n.visits.season)
    for(i in 1:n.seasons) {
        col.start <- vis.seq[i]
        col.end <- col.start + (n.visits.season - 1)
        ySeason <- yMat[, col.start:col.end]
        ##summarize count histories
        if(is.null(ncol(ySeason))){
            ySeason <- as.matrix(ySeason)
        }
        yVec.season <- as.vector(ySeason)
        yVectors[[i]] <- yVec.season
        
        det.hist <- apply(X = ySeason, MARGIN = 1, FUN = function(i) paste(i, collapse = "|"))
        hist.table.seasons[[i]] <- table(det.hist, deparse.level = 0)
        count.table.seasons[[i]] <- table(yVec.season, exclude = NULL)

        ##determine proportion of sites with at least 1 detection
        det.sum <- apply(X = ySeason, MARGIN = 1, FUN = function(i) ifelse(sum(i, na.rm = TRUE) > 0, 1, 0))

        ##check sites with observed detections and deal with NA's
        sum.rows <- rowSums(ySeason, na.rm = TRUE)
        is.na(sum.rows) <- rowSums(is.na(ySeason)) == ncol(ySeason)
    
        ##number of sites sampled
        out.freqs[i, 1] <- sum(!is.na(sum.rows))
        out.freqs[i, 2] <- sum(det.sum)


        ##sites without detections
        none <- which(sum.rows == 0)
        ##sites with at least one detection
        some <- which(sum.rows != 0) 
        out.seasons[[i]] <- list("none" = none, "some" = some)
    }

    ##check if any seasons were not sampled
    y.seasonsNA <- sapply(yVectors, FUN = function(i) all(is.na(i)))

    ##if only season-specific plots are requested
    if(!plot.freq && plot.seasons) {
        ##determine arrangement of plots in matrix
        if(n.seasons >= 12) {
            n.seasons.adj <- 12
            warning("\nOnly first 12 seasons are plotted\n")
        }
        
        if(n.seasons.adj <= 12) {

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
        }
        ##reset graphics parameters and save in object
        oldpar <- par(mfrow = c(nRows, nCols))
    }
    

    ##if both plots for seasons and combined are requested
    ##summarize detection histories
    if(plot.freq) {

        if(!plot.seasons) {
            nRows <- 1
            nCols <- 1
        }
        
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

        ##histogram for data combined across seasons
        barplot(table(yVec.full), ylab = "Frequency", xlab = "Counts of individuals",
                main = paste("Distribution of raw counts (", n.seasons, " seasons combined)", sep = ""),
                cex.axis = cex.axis, cex.names = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
    }


    ##iterate over each season
    if(plot.seasons) {
        for(k in 1:n.seasons.adj) {
            ##check for missing season
            if(y.seasonsNA[k]) {next} #skip to next season if current season not sampled
            ##histogram for data combined across seasons
            barplot(table(yVectors[[k]]), ylab = "Frequency", xlab = "Counts of individuals",
                    main = paste("Distribution of raw counts (season ", k, ")", sep = ""),
                    cex.axis = cex.axis, cex.names = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
        }
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

    ##reset to original values
    if(any(plot.freq || plot.seasons)) {
        on.exit(par(oldpar))
    }
   
    out.count <- list("count.table.full" = count.table.full,
                      "count.table.seasons" = count.table.seasons,
                      "hist.table.full" = hist.table.full,
                      "hist.table.seasons" = hist.table.seasons,
                      "out.freqs" = out.freqs, "out.props" = out.props,
                      "n.seasons" = n.seasons,
                      "n.visits.season" = n.visits.season,
                      "missing.seasons" = y.seasonsNA)
    class(out.count) <- "countHist"
    return(out.count)
}



##for unmarkedFitMMO
countHist.unmarkedFitMMO <- function(object, plot.freq = TRUE, 
                                     cex.axis = 1, cex.lab = 1, cex.main = 1,
                                     plot.seasons = FALSE, ...) {

    ##extract data
    yMat <- object@data@y
    nsites <- nrow(yMat)
    n.seasons <- object@data@numPrimary
    n.seasons.adj <- n.seasons #total number of plots fixed to 11 or 12, depending on plots requested
    nvisits <- ncol(yMat)
    ##visits per season
    n.visits.season <- nvisits/n.seasons
    seasonNames <- paste("season", 1:n.seasons, sep = "")

    ##collapse yMat into a single vector
    yVec.full <- as.vector(yMat)
  
    ##raw counts
    count.table.full <- table(yVec.full, exclude = NULL, deparse.level = 0)
  
    ##summarize count histories
    hist.full <- apply(X = yMat, MARGIN = 1, FUN = function(i) paste(i, collapse = "|"))
    hist.table.full <- table(hist.full, deparse.level = 0)

    ##for each season, determine frequencies
    yVectors <- vector(mode = "list", length = n.seasons)
    out.seasons <- vector(mode = "list", length = n.seasons)
    count.table.seasons <- vector(mode = "list", length = n.seasons)
    names(count.table.seasons) <- seasonNames

    hist.table.seasons <- vector(mode = "list", length = n.seasons)
    names(hist.table.seasons) <- seasonNames
    
    out.freqs <- matrix(data = NA, ncol = 6, nrow = n.seasons)
    colnames(out.freqs) <- c("sampled", "detected", "colonized",
                             "extinct", "static", "common")
    rownames(out.freqs) <- paste("Season-", 1:n.seasons, sep = "")

    ##sequence of visits
    vis.seq <- seq(from = 1, to = nvisits, by = n.visits.season)
    for(i in 1:n.seasons) {
        col.start <- vis.seq[i]
        col.end <- col.start + (n.visits.season - 1)
        ySeason <- yMat[, col.start:col.end]
        ##summarize count histories
        if(is.null(ncol(ySeason))){
            ySeason <- as.matrix(ySeason)
        }
        yVec.season <- as.vector(ySeason)
        yVectors[[i]] <- yVec.season
        
        det.hist <- apply(X = ySeason, MARGIN = 1, FUN = function(i) paste(i, collapse = "|"))
        hist.table.seasons[[i]] <- table(det.hist, deparse.level = 0)
        count.table.seasons[[i]] <- table(yVec.season, exclude = NULL)
        
        ##determine proportion of sites with at least 1 detection
        det.sum <- apply(X = ySeason, MARGIN = 1, FUN = function(i) ifelse(sum(i, na.rm = TRUE) > 0, 1, 0))

        ##check sites with observed detections and deal with NA's
        sum.rows <- rowSums(ySeason, na.rm = TRUE)
        is.na(sum.rows) <- rowSums(is.na(ySeason)) == ncol(ySeason)
    
        ##number of sites sampled
        out.freqs[i, 1] <- sum(!is.na(sum.rows))
        out.freqs[i, 2] <- sum(det.sum)


        ##sites without detections
        none <- which(sum.rows == 0)
        ##sites with at least one detection
        some <- which(sum.rows != 0) 
        out.seasons[[i]] <- list("none" = none, "some" = some)
    }
    
    ##check if any seasons were not sampled
    y.seasonsNA <- sapply(yVectors, FUN = function(i) all(is.na(i)))

    ##if only season-specific plots are requested
    if(!plot.freq && plot.seasons) {
        ##determine arrangement of plots in matrix
        if(n.seasons >= 12) {
            n.seasons.adj <- 12
            warning("\nOnly first 12 seasons are plotted\n")
        }
    
        if(n.seasons.adj <= 12) {

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
        }
        ##reset graphics parameters and save in object
        oldpar <- par(mfrow = c(nRows, nCols))
    }

    
    ##if both plots for seasons and combined are requested
    ##summarize detection histories
    if(plot.freq) {

        if(!plot.seasons) {
            nRows <- 1
            nCols <- 1
        }
                
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

        ##histogram for data combined across seasons
        barplot(table(yVec.full), ylab = "Frequency", xlab = "Counts of individuals",
                main = paste("Distribution of raw counts (", n.seasons, " seasons combined)", sep = ""),
                cex.axis = cex.axis, cex.names = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
    }

    ##iterate over each season
    if(plot.seasons) {
        for(k in 1:n.seasons.adj) {
            ##check for missing season
            if(y.seasonsNA[k]) {next} #skip to next season if current season not sampled
            ##histogram for data combined across seasons
            barplot(table(yVectors[[k]]), ylab = "Frequency", xlab = "Counts of individuals",
                    main = paste("Distribution of raw counts (season ", k, ")", sep = ""),
                    cex.axis = cex.axis, cex.names = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
        }
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

    ##reset to original values
    if(any(plot.freq || plot.seasons)) {
        on.exit(par(oldpar))
    }
  
    out.count <- list("count.table.full" = count.table.full,
                      "count.table.seasons" = count.table.seasons,
                      "hist.table.full" = hist.table.full,
                      "hist.table.seasons" = hist.table.seasons,
                      "out.freqs" = out.freqs, "out.props" = out.props,
                      "n.seasons" = n.seasons,
                      "n.visits.season" = n.visits.season,
                      "missing.seasons" = y.seasonsNA)
    class(out.count) <- "countHist"
    return(out.count)
}



##print method
print.countHist <- function(x, digits = 2, ...) {
    if(x$n.seasons == 1) {

        ##convert NA to . for nicer printing
        hist.names <- names(x$hist.table.full)
        names(x$hist.table.full) <- gsub(pattern = "NA",
                                          replacement = ".",
                                          x = hist.names)
        
        cat("\nSummary of counts:\n")
        count.mat <- matrix(x$count.table.full, nrow = 1)
        colnames(count.mat) <- names(x$count.table.full)
        rownames(count.mat) <- "Frequency"
        print(count.mat)
    
        cat("\nSummary of count histories:\n")
        ##account for number of visits, number of unique histories, number of separators
        num.chars <- nchar(paste(names(x$hist.table.full), collapse = ""))
        if(num.chars >= 80) {
            cat("\nNote:  Count histories exceed 80 characters and are not displayed\n")
        } else {
            out.mat <- matrix(x$hist.table.full, nrow = 1)
            colnames(out.mat) <- names(x$hist.table.full)
            rownames(out.mat) <- "Frequency"
            print(out.mat)
        }
    
        cat("\nProportion of sites with at least one detection:\n", round(x$out.props[, "naive.occ"], digits), "\n\n")
        
        cat("Frequencies of sites with detections:\n")
        ##add matrix of frequencies
        print(x$out.freqs)

    }

    if(x$n.seasons > 1) {
        
        cat("\nSummary of counts (", x$n.seasons, " seasons combined): \n", sep ="")
        count.mat <- matrix(x$count.table.full, nrow = 1)
        colnames(count.mat) <- names(x$count.table.full)
        rownames(count.mat) <- "Frequency"
        print(count.mat)
    
        cat("\nSummary of count histories:\n")

        if(x$n.visits.season == 1) {
            visits <- 1
        } else {
            visits <- x$n.visits.season - 1
        }
        
        num.chars <- nchar(paste(names(x$hist.table.full), collapse = ""))
        
        if(num.chars >= 80) {
            cat("\nNote:  Count histories exceed 80 characters and are not displayed\n")
        } else {
            out.mat <- matrix(x$hist.table.full, nrow = 1)
            colnames(out.mat) <- names(x$hist.table.full)
            rownames(out.mat) <- "Frequency"
            print(out.mat)
        }

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
        cat("\nSeason-specific counts: \n")
        cat("\n")
        for(i in 1:x$n.seasons) {
            if(!x$missing.seasons[i]) {
                cat("Season", i, "\n")
            } else {
                cat("Season", i, "(no sites sampled)", "\n")
            }
            
            temp.tab <- x$count.table.seasons[[i]]
            out.mat <- matrix(temp.tab, nrow = 1)
            colnames(out.mat) <- names(temp.tab)
            rownames(out.mat) <- "Frequency"
            print(out.mat)
            cat("--------\n\n")
        }
    
        cat("Frequencies of sites with detections, extinctions, and colonizations:\n")
        ##add matrix of frequencies
        print(x$out.freqs)
        
    }
}

