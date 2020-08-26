##summarize detection histories and count data
countDist <- function(object, plot.freq = TRUE, plot.distance = TRUE,
                      cex.axis = 1, cex.lab = 1, cex.main = 1, ...){
  UseMethod("countDist", object)
}



countDist.default <- function(object, plot.freq = TRUE, plot.distance = TRUE,
                              cex.axis = 1, cex.lab = 1, cex.main = 1, ...){
  stop("\nFunction not yet defined for this object class\n")
}




##for unmarkedFrameDS
countDist.unmarkedFrameDS <- function(object, plot.freq = TRUE, plot.distance = TRUE,
                                      cex.axis = 1, cex.lab = 1, cex.main = 1, ...) {

  ##extract data
  yMat <- object@y
  nsites <- nrow(yMat)
  n.seasons <- 1
  nvisits <- 1
    ##visits per season
    n.visits.season <- 1
    
  ##distance classes
  dist.classes <- object@dist.breaks

  ##number of distance classes
  n.dist.classes <- length(dist.classes) - 1
    
  ##units
  unitsIn <- object@unitsIn
  
  ##create string of names
  dist.names <- rep(NA, n.dist.classes)
  for(i in 1:n.dist.classes){
    dist.names[i] <- paste(dist.classes[i], "-", dist.classes[i+1], sep = "")
  }
  
  ##collapse yMat into a single vector
  yVec <- as.vector(yMat)

  ##determine size of plot window
  ##when both types are requested
  if(plot.freq && plot.distance) {
      par(mfrow = c(1, 2))
  }
  
  ##summarize counts
  if(plot.freq) {
    
    ##create histogram
    barplot(table(yVec), ylab = "Frequency", xlab = "Counts of individuals",
            main = "Distribution of raw counts",
            cex.axis = cex.axis, cex.names = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
  }

  ##summarize counts per distance
  dist.sums.full <- colSums(yMat)
  names(dist.sums.full) <- dist.names
  dist.table.seasons <- list(dist.sums.full)
  names(dist.table.seasons)  <- "season1"
    
  if(plot.distance) {
    
    ##create histogram
    barplot(dist.sums.full, ylab = "Frequency",
            xlab = paste("Distance class (", unitsIn, ")", sep = ""),
            main = "Distribution of distance data",
            cex.axis = cex.axis, cex.names = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
  }

  ##raw counts
  count.table.full <- table(yVec, exclude = NULL, deparse.level = 0)
  count.table.seasons <- list(count.table.full)


  ##for each season, determine frequencies
  out.freqs <- matrix(data = NA, ncol = 2, nrow = n.seasons)
  colnames(out.freqs) <- c("sampled", "detected")
  rownames(out.freqs) <- "Season-1"

  ySeason <- yMat
    
  ##determine proportion of sites with at least 1 detection
  det.sum <- apply(X = ySeason, MARGIN = 1, FUN = function(i) ifelse(sum(i, na.rm = TRUE) > 0, 1, 0))

  ##check sites with observed detections and deal with NA's
  sum.rows <- rowSums(ySeason, na.rm = TRUE)
  is.na(sum.rows) <- rowSums(is.na(ySeason)) == ncol(ySeason)
    
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
                    "dist.sums.full" = dist.sums.full,
                    "dist.table.seasons" = dist.table.seasons,
                    "dist.names" = dist.names,
                    "n.dist.classes" = n.dist.classes,
                    "out.freqs" = out.freqs,
                    "out.props" = out.props,
                    "n.seasons" = n.seasons,
                    "n.visits.season" = n.visits.season)
  class(out.count) <- "countDist"
  return(out.count)
}



##for unmarkedFitDS
countDist.unmarkedFitDS <- function(object, plot.freq = TRUE, plot.distance = TRUE,
                                    cex.axis = 1, cex.lab = 1, cex.main = 1, ...) {

    ##extract data
    yMat <- object@data@y
    nsites <- nrow(yMat)
    n.seasons <- 1
    nvisits <- 1
    ##visits per season
    n.visits.season <- 1
    
    ##distance classes
    dist.classes <- object@data@dist.breaks
    
    ##number of distance classes
    n.dist.classes <- length(dist.classes) - 1

    ##units
    unitsIn <- object@data@unitsIn
  
    ##create string of names
    dist.names <- rep(NA, n.dist.classes)
    for(i in 1:n.dist.classes){
        dist.names[i] <- paste(dist.classes[i], "-", dist.classes[i+1], sep = "")
    }
  
    ##collapse yMat into a single vector
    yVec <- as.vector(yMat)

    ##determine size of plot window
    ##when both types are requested
    if(plot.freq && plot.distance) {
        par(mfrow = c(1, 2))
    }
  
    ##summarize counts
    if(plot.freq) {
        
        ##create histogram
        barplot(table(yVec), ylab = "Frequency", xlab = "Counts of individuals",
                main = "Distribution of raw counts",
                cex.axis = cex.axis, cex.names = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
    }
    
    ##summarize counts per distance
    dist.sums.full <- colSums(yMat)
    names(dist.sums.full) <- dist.names
    dist.table.seasons <- list(dist.sums.full)
    names(dist.table.seasons) <- "season1"
    
    if(plot.distance) {
    
        ##create histogram
        barplot(dist.sums.full, ylab = "Frequency",
                xlab = paste("Distance class (", unitsIn, ")", sep = ""),
                main = "Distribution of distance data",
                cex.axis = cex.axis, cex.names = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
    }

    ##raw counts
    count.table.full <- table(yVec, exclude = NULL, deparse.level = 0)
    count.table.seasons <- list(count.table.full)


    ##for each season, determine frequencies
    out.freqs <- matrix(data = NA, ncol = 2, nrow = n.seasons)
    colnames(out.freqs) <- c("sampled", "detected")
    rownames(out.freqs) <- "Season-1"

    ySeason <- yMat
    
    ##determine proportion of sites with at least 1 detection
    det.sum <- apply(X = ySeason, MARGIN = 1, FUN = function(i) ifelse(sum(i, na.rm = TRUE) > 0, 1, 0))
    
    ##check sites with observed detections and deal with NA's
    sum.rows <- rowSums(ySeason, na.rm = TRUE)
    is.na(sum.rows) <- rowSums(is.na(ySeason)) == ncol(ySeason)
    
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
                      "dist.sums.full" = dist.sums.full,
                      "dist.table.seasons" = dist.table.seasons,
                      "dist.names" = dist.names,
                      "n.dist.classes" = n.dist.classes,
                      "out.freqs" = out.freqs,
                      "out.props" = out.props,
                      "n.seasons" = n.seasons,
                      "n.visits.season" = n.visits.season)
    class(out.count) <- "countDist"
    return(out.count)
}



##for unmarkedFrameGDS
countDist.unmarkedFrameGDS <- function(object, plot.freq = TRUE, plot.distance = TRUE,
                                       cex.axis = 1, cex.lab = 1, cex.main = 1, ...) {

  ##extract data
  yMat <- object@y
  nsites <- nrow(yMat)
  n.seasons <- 1
  ##for GDS - several visits in single season
  nvisits <- object@numPrimary
  ##visits per season
  n.visits.season <- 1
  ##distance classes
  dist.classes <- object@dist.breaks

  ##number of distance classes
  n.dist.classes <- length(dist.classes) - 1

  ##units
  unitsIn <- object@unitsIn
  
  ##create string of names
  dist.names <- rep(NA, n.dist.classes)
  for(i in 1:n.dist.classes){
    dist.names[i] <- paste(dist.classes[i], "-", dist.classes[i+1], sep = "")
  }
  
  ##collapse yMat into a single vector
  yVec <- as.vector(yMat)

  ##determine size of plot window
  ##when both types are requested
  if(plot.freq && plot.distance) {
      par(mfrow = c(1, 2))
  }
  
  ##summarize counts
  if(plot.freq) {
    
    ##create histogram
    barplot(table(yVec), ylab = "Frequency", xlab = "Counts of individuals",
            main = "Distribution of raw counts",
            cex.axis = cex.axis, cex.names = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
  }

  ##summarize counts per distance
  dist.sums.full <- rep(NA, n.dist.classes)
  ##create matrix to hold indices of visits x dist.classes
  mat.dist <- matrix(1:(nvisits*n.dist.classes),
                     nrow = n.dist.classes,
                     ncol = nvisits)
  for(j in 1:n.dist.classes) {
    dist.sums.full[j] <- sum(colSums(yMat[, mat.dist[j, ]], na.rm = TRUE), na.rm = TRUE)
  }

    names(dist.sums.full) <- dist.names
    dist.table.seasons <- list(dist.sums.full)
    names(dist.table.seasons) <- "season1"
    
  if(plot.distance) {
    
    ##create histogram
    barplot(dist.sums.full, ylab = "Frequency",
            xlab = paste("Distance class (", unitsIn, ")", sep = ""),
            main = "Distribution of distance data",
            cex.axis = cex.axis, cex.names = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
  }

  ##raw counts
  count.table.full <- table(yVec, exclude = NULL, deparse.level = 0)
  count.table.seasons <- list(count.table.full)


  ##for each season, determine frequencies
  out.freqs <- matrix(data = NA, ncol = 2, nrow = n.seasons)
  colnames(out.freqs) <- c("sampled", "detected")
  rownames(out.freqs) <- "Season-1"

  ySeason <- yMat
    
  ##determine proportion of sites with at least 1 detection
  det.sum <- apply(X = ySeason, MARGIN = 1, FUN = function(i) ifelse(sum(i, na.rm = TRUE) > 0, 1, 0))

  ##check sites with observed detections and deal with NA's
  sum.rows <- rowSums(ySeason, na.rm = TRUE)
  is.na(sum.rows) <- rowSums(is.na(ySeason)) == ncol(ySeason)
    
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
                    "dist.sums.full" = dist.sums.full,
                    "dist.table.seasons" = dist.table.seasons,
                    "dist.names" = dist.names,
                    "n.dist.classes" = n.dist.classes,
                    "out.freqs" = out.freqs,
                    "out.props" = out.props,
                    "n.seasons" = n.seasons,
                    "n.visits.season" = n.visits.season)
  class(out.count) <- "countDist"
  return(out.count)
}



##for unmarkedFitGDS
countDist.unmarkedFitGDS <- function(object, plot.freq = TRUE, plot.distance = TRUE,
                                     cex.axis = 1, cex.lab = 1, cex.main = 1, ...) {

  ##extract data
  yMat <- object@data@y
  nsites <- nrow(yMat)
  n.seasons <- 1
  ##for GDS - several visits in single season
  nvisits <- object@data@numPrimary
  ##visits per season
  n.visits.season <- 1

  ##distance classes
  dist.classes <- object@data@dist.breaks

  ##number of distance classes
  n.dist.classes <- length(dist.classes) - 1

  ##units
  unitsIn <- object@data@unitsIn
  
  ##create string of names
  dist.names <- rep(NA, n.dist.classes)
  for(i in 1:n.dist.classes){
    dist.names[i] <- paste(dist.classes[i], "-", dist.classes[i+1], sep = "")
  }
  

  ##collapse yMat into a single vector
  yVec <- as.vector(yMat)

  ##determine size of plot window
  ##when both types are requested
  if(plot.freq && plot.distance) {
      par(mfrow = c(1, 2))
  }
  
  ##summarize counts
  if(plot.freq) {
    
    ##create histogram
    barplot(table(yVec), ylab = "Frequency", xlab = "Counts of individuals",
            main = "Distribution of raw counts",
            cex.axis = cex.axis, cex.names = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
  }

  ##summarize counts per distance
  dist.sums.full <- rep(NA, n.dist.classes)
  ##create matrix to hold indices of visits x dist.classes
  mat.dist <- matrix(1:(nvisits*n.dist.classes),
                     nrow = n.dist.classes,
                     ncol = nvisits)
  for(j in 1:n.dist.classes) {
      dist.sums.full[j] <- sum(colSums(yMat[, mat.dist[j, ]], na.rm = TRUE), na.rm = TRUE)
  }  

    names(dist.sums.full) <- dist.names
    dist.table.seasons <- list(dist.sums.full)
    names(dist.table.seasons) <- "season1"
    
  if(plot.distance) {
    
    ##create histogram
    barplot(dist.sums.full, ylab = "Frequency",
            xlab = paste("Distance class (", unitsIn, ")", sep = ""),
            main = "Distribution of distance data",
            cex.axis = cex.axis, cex.names = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
  }

  ##raw counts
  count.table.full <- table(yVec, exclude = NULL, deparse.level = 0)
  count.table.seasons <- list(count.table.full)


  ##for each season, determine frequencies
  out.freqs <- matrix(data = NA, ncol = 2, nrow = n.seasons)
  colnames(out.freqs) <- c("sampled", "detected")
  rownames(out.freqs) <- "Season-1"

  ySeason <- yMat
    
  ##determine proportion of sites with at least 1 detection
  det.sum <- apply(X = ySeason, MARGIN = 1, FUN = function(i) ifelse(sum(i, na.rm = TRUE) > 0, 1, 0))

  ##check sites with observed detections and deal with NA's
  sum.rows <- rowSums(ySeason, na.rm = TRUE)
  is.na(sum.rows) <- rowSums(is.na(ySeason)) == ncol(ySeason)
    
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
                    "dist.sums.full" = dist.sums.full,
                    "dist.table.seasons" = dist.table.seasons,
                    "dist.names" = dist.names,
                    "n.dist.classes" = n.dist.classes,
                    "out.freqs" = out.freqs,
                    "out.props" = out.props,
                    "n.seasons" = n.seasons,
                    "n.visits.season" = n.visits.season)
  class(out.count) <- "countDist"
  return(out.count)
}



##for unmarkedFrameDSO
countDist.unmarkedFrameDSO <- function(object, plot.freq = TRUE, plot.distance = TRUE,
                                       cex.axis = 1, cex.lab = 1, cex.main = 1,
                                       plot.seasons = FALSE, ...) {

    ##extract data
    yMat <- object@y
    nsites <- nrow(yMat)
    n.seasons <- object@numPrimary
    n.seasons.adj <- n.seasons #total number of plots fixed to 11 or 12, depending on plots requested
    n.columns <- ncol(yMat)
    seasonNames <- paste("season", 1:n.seasons, sep = "")
    ##visits per season
    n.visits.season <- 1

    ##distance breaks
    dist.breaks <- object@dist.breaks

    ##number of distance classes
    n.dist.classes <- length(dist.breaks) - 1
    
    ##units
    unitsIn <- object@unitsIn
  
    ##create string of names
    dist.names <- rep(NA, n.dist.classes)
    for(i in 1:n.dist.classes){
        dist.names[i] <- paste(dist.breaks[i], "-", dist.breaks[i+1], sep = "")
    }

    ##determine size of plot window
    ##when two types are requested
    if(plot.freq && plot.distance && !plot.seasons) {
        par(mfrow = c(1, 2))
    }

    if(plot.freq && !plot.distance && !plot.seasons) {
        par(mfrow = c(1, 1))
    }

    if(!plot.freq && plot.distance && !plot.seasons) {
        par(mfrow = c(1, 1))
    }


    ##if only season-specific plots are requested
    if(plot.seasons) {

        if(!plot.freq && !plot.distance) {
            ##determine arrangement of plots in matrix
            if(plot.seasons && n.seasons >= 12) {
                n.seasons.adj <- 12
                warning("\nOnly first 12 seasons are plotted\n")
            }
        }

        if(!plot.freq && !plot.distance) {
            ##determine arrangement of plots in matrix
            if(plot.seasons && n.seasons >= 12) {
                n.seasons.adj <- 12
                warning("\nOnly first 12 seasons are plotted\n")
            }
        }

        
        if(plot.freq && !plot.distance || !plot.freq && plot.distance) {
            ##determine arrangement of plots in matrix
            if(plot.seasons && n.seasons >= 11) {
                n.seasons.adj <- 11
                warning("\nOnly first 11 seasons are plotted\n")
            }
        }

        if(plot.freq && plot.distance) {
            ##determine arrangement of plots in matrix
            if(plot.seasons && n.seasons >= 10) {
                n.seasons.adj <- 10
                warning("\nOnly first 10 seasons are plotted\n")
            }
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
        }

        par(mfrow = c(nRows, nCols))
    }
    
    ##distances across years
    yVec <- as.vector(yMat)

    ##summarize counts
    if(plot.freq) {

        ##create histogram
        barplot(table(yVec), ylab = "Frequency", xlab = "Counts of individuals",
                main = paste("Distribution of raw counts (", n.seasons, " seasons combined)", sep = ""),
                cex.axis = cex.axis, cex.names = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
    }

 
    ##columns for each season
    col.seasons <- seq(1, n.columns, by = n.dist.classes)
    yMat.seasons <- vector(mode = "list", length = n.seasons)
    names(yMat.seasons) <- seasonNames
    minusOne <- n.dist.classes - 1
    
    ##iterate over each season
    for(i in 1:n.seasons) {
        
        ##extract yMat for each year
        yMat1 <- yMat[, col.seasons[i]:(col.seasons[i]+minusOne)]

        yMat.seasons[[i]] <- yMat1

    }


    ##collapse list into matrix
    fullData <- do.call("rbind", yMat.seasons)
    ##summarize counts per distance
    dist.sums.full <- colSums(fullData)
    names(dist.sums.full) <- dist.names
    
    if(plot.distance) {
    
            ##create histogram
            barplot(dist.sums.full, ylab = "Frequency",
                    xlab = paste("Distance class (", unitsIn, ")", sep = ""),
                    main = paste("Distribution of distance data (", n.seasons, " seasons combined)", sep = ""),
                    cex.axis = cex.axis, cex.names = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
        }


    dist.table.seasons <- vector("list", n.seasons)
    names(dist.table.seasons) <- seasonNames
    
    for(i in 1:n.seasons) {
        
        yMat2 <- yMat.seasons[[i]]
        
        ##summarize counts per distance
        dist.sums.season <- colSums(yMat2)
        names(dist.sums.season) <- dist.names
        dist.table.seasons[[i]] <- dist.sums.season
  
        if(plot.seasons) {
    
            ##create histogram
            barplot(dist.sums.season, ylab = "Frequency",
                    xlab = paste("Distance class (", unitsIn, ")", sep = ""),
                    main = paste("Distribution of distance data (season ", i, ")", sep = ""),
                    cex.axis = cex.axis, cex.names = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
        }
    }

    
  ##raw counts
  count.table.full <- table(yVec, exclude = NULL, deparse.level = 0)
  count.table.seasons <- lapply(yMat.seasons, table)

        
    ##for each season, determine frequencies
    out.seasons <- vector(mode = "list", length = n.seasons)
    out.freqs <- matrix(data = NA, ncol = 6, nrow = n.seasons)
    colnames(out.freqs) <- c("sampled", "detected", "colonized",
                             "extinct", "static", "common")
    rownames(out.freqs) <- paste("Season-", 1:n.seasons, sep = "")
    
    for(i in 1:n.seasons) {

        ySeason <- yMat.seasons[[i]]
        ##determine proportion of sites with at least 1 detection
        det.sum <- apply(X = ySeason, MARGIN = 1, FUN = function(i) ifelse(sum(i, na.rm = TRUE) > 0, 1, 0))

        ##check sites with observed detections and deal with NA's
        sum.rows <- rowSums(ySeason, na.rm = TRUE)
        is.na(sum.rows) <- rowSums(is.na(ySeason)) == ncol(ySeason)
    
        ##number of sites sampled
        out.freqs[i, 1] <- sum(!is.na(sum.rows))
        ##number of sites with at least 1 detection
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

    ##create a matrix with proportion of sites with colonizations
    ##and extinctions based on raw data
    out.props <- matrix(NA, nrow = nrow(out.freqs), ncol = 4)
    colnames(out.props) <- c("naive.occ", "naive.colonization", "naive.extinction", "naive.static")
    rownames(out.props) <- rownames(out.freqs)
    out.props[, 1] <- out.freqs[, 2]/out.freqs[, 1]
    out.props[, 2] <- out.freqs[, 3]/out.freqs[, 6]
    out.props[, 3] <- out.freqs[, 4]/out.freqs[, 6]
    out.props[, 4] <- out.freqs[, 5]/out.freqs[, 6]

    
    out.count <- list("count.table.full" = count.table.full,
                      "count.table.seasons" = count.table.seasons,
                      "dist.sums.full" = dist.sums.full,
                      "dist.table.seasons" = dist.table.seasons,
                      "dist.names" = dist.names,
                      "n.dist.classes" = n.dist.classes,
                      "out.freqs" = out.freqs,
                      "out.props" = out.props,
                      "n.seasons" = n.seasons,
                      "n.visits.season" = n.visits.season)
    class(out.count) <- "countDist"
    return(out.count)
}



##for unmarkedFitDSO
countDist.unmarkedFitDSO <- function(object, plot.freq = TRUE, plot.distance = TRUE,
                                     cex.axis = 1, cex.lab = 1, cex.main = 1,
                                     plot.seasons = FALSE, ...) {

    ##extract data
    yMat <- object@data@y
    nsites <- nrow(yMat)
    n.seasons <- object@data@numPrimary
    n.seasons.adj <- n.seasons #total number of plots fixed to 11 or 12, depending on plots requested
    n.columns <- ncol(yMat)
    seasonNames <- paste("season", 1:n.seasons, sep = "")
    ##visits per season
    n.visits.season <- 1

    ##distance breaks
    dist.breaks <- object@data@dist.breaks

    ##number of distance classes
    n.dist.classes <- length(dist.breaks) - 1
    
    ##units
    unitsIn <- object@data@unitsIn
  
    ##create string of names
    dist.names <- rep(NA, n.dist.classes)
    for(i in 1:n.dist.classes){
        dist.names[i] <- paste(dist.breaks[i], "-", dist.breaks[i+1], sep = "")
    }

    ##determine size of plot window
    ##when two types are requested
    if(plot.freq && plot.distance && !plot.seasons) {
        par(mfrow = c(1, 2))
    }

    if(plot.freq && !plot.distance && !plot.seasons) {
        par(mfrow = c(1, 1))
    }

    if(!plot.freq && plot.distance && !plot.seasons) {
        par(mfrow = c(1, 1))
    }


    ##if only season-specific plots are requested
    if(plot.seasons) {

        if(!plot.freq && !plot.distance) {
            ##determine arrangement of plots in matrix
            if(plot.seasons && n.seasons >= 12) {
                n.seasons.adj <- 12
                warning("\nOnly first 12 seasons are plotted\n")
            }
        }

        
        if(plot.freq && !plot.distance || !plot.freq && plot.distance) {
            ##determine arrangement of plots in matrix
            if(plot.seasons && n.seasons >= 11) {
                n.seasons.adj <- 11
                warning("\nOnly first 11 seasons are plotted\n")
            }
        }

        if(plot.freq && plot.distance) {
            ##determine arrangement of plots in matrix
            if(plot.seasons && n.seasons >= 10) {
                n.seasons.adj <- 10
                warning("\nOnly first 10 seasons are plotted\n")
            }
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
        }

        par(mfrow = c(nRows, nCols))
    }

    ##distances across years
    yVec <- as.vector(yMat)

    ##summarize counts
    if(plot.freq) {

        ##create histogram
        barplot(table(yVec), ylab = "Frequency", xlab = "Counts of individuals",
                main = paste("Distribution of raw counts (", n.seasons, " seasons combined)", sep = ""),
                cex.axis = cex.axis, cex.names = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
    }

 
    ##columns for each season
    col.seasons <- seq(1, n.columns, by = n.dist.classes)
    yMat.seasons <- vector(mode = "list", length = n.seasons)
    names(yMat.seasons) <- seasonNames
    minusOne <- n.dist.classes - 1
    ##iterate over each season
    for(i in 1:n.seasons) {
        
        ##extract yMat for each season
        yMat1 <- yMat[, col.seasons[i]:(col.seasons[i]+minusOne)]

        yMat.seasons[[i]] <- yMat1

    }


    ##collapse list into matrix
    fullData <- do.call("rbind", yMat.seasons)
    ##summarize counts per distance
    dist.sums.full <- colSums(fullData)
    names(dist.sums.full) <- dist.names
    
    if(plot.distance) {
    
            ##create histogram
            barplot(dist.sums.full, ylab = "Frequency",
                    xlab = paste("Distance class (", unitsIn, ")", sep = ""),
                    main = paste("Distribution of distance data (", n.seasons, " seasons combined)", sep = ""),
                    cex.axis = cex.axis, cex.names = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
        }


    dist.table.seasons <- vector("list", n.seasons)
    names(dist.table.seasons) <- seasonNames
    
    for(i in 1:n.seasons) {
        
        yMat2 <- yMat.seasons[[i]]
        
        ##summarize counts per distance
        dist.sums.season <- colSums(yMat2)
        names(dist.sums.season) <- dist.names
        dist.table.seasons[[i]] <- dist.sums.season
  
        if(plot.seasons) {
    
            ##create histogram
            barplot(dist.sums.season, ylab = "Frequency",
                    xlab = paste("Distance class (", unitsIn, ")", sep = ""),
                    main = paste("Distribution of distance data (season ", i, ")", sep = ""),
                    cex.axis = cex.axis, cex.names = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
        }
    }

    
  ##raw counts
  count.table.full <- table(yVec, exclude = NULL, deparse.level = 0)
  count.table.seasons <- lapply(yMat.seasons, table)

        
    ##for each season, determine frequencies
    out.seasons <- vector(mode = "list", length = n.seasons)
    out.freqs <- matrix(data = NA, ncol = 6, nrow = n.seasons)
    colnames(out.freqs) <- c("sampled", "detected", "colonized",
                             "extinct", "static", "common")
    rownames(out.freqs) <- paste("Season-", 1:n.seasons, sep = "")
    
    for(i in 1:n.seasons) {

        ySeason <- yMat.seasons[[i]]
        ##determine proportion of sites with at least 1 detection
        det.sum <- apply(X = ySeason, MARGIN = 1, FUN = function(i) ifelse(sum(i, na.rm = TRUE) > 0, 1, 0))

        ##check sites with observed detections and deal with NA's
        sum.rows <- rowSums(ySeason, na.rm = TRUE)
        is.na(sum.rows) <- rowSums(is.na(ySeason)) == ncol(ySeason)
    
        ##number of sites sampled
        out.freqs[i, 1] <- sum(!is.na(sum.rows))
        ##number of sites with at least 1 detection
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

    ##create a matrix with proportion of sites with colonizations
    ##and extinctions based on raw data
    out.props <- matrix(NA, nrow = nrow(out.freqs), ncol = 4)
    colnames(out.props) <- c("naive.occ", "naive.colonization", "naive.extinction", "naive.static")
    rownames(out.props) <- rownames(out.freqs)
    out.props[, 1] <- out.freqs[, 2]/out.freqs[, 1]
    out.props[, 2] <- out.freqs[, 3]/out.freqs[, 6]
    out.props[, 3] <- out.freqs[, 4]/out.freqs[, 6]
    out.props[, 4] <- out.freqs[, 5]/out.freqs[, 6]

    
    out.count <- list("count.table.full" = count.table.full,
                      "count.table.seasons" = count.table.seasons,
                      "dist.sums.full" = dist.sums.full,
                      "dist.table.seasons" = dist.table.seasons,
                      "dist.names" = dist.names,
                      "n.dist.classes" = n.dist.classes,
                      "out.freqs" = out.freqs,
                      "out.props" = out.props,
                      "n.seasons" = n.seasons,
                      "n.visits.season" = n.visits.season)
    class(out.count) <- "countDist"
    return(out.count)
}



##print method
print.countDist <- function(x, digits = 2, ...) {
  if(x$n.seasons == 1) {
    cat("\nSummary of counts:\n")
    count.mat <- matrix(x$count.table.full, nrow = 1)
    colnames(count.mat) <- names(x$count.table.full)
    rownames(count.mat) <- "Frequency"
    print(count.mat)
    
    cat("\nSummary of distance data:\n")
    out.mat <- matrix(x$dist.sums.full, nrow = 1)
    colnames(out.mat) <- names(x$dist.sums.full)
    rownames(out.mat) <- "Frequency"
    print(out.mat)

    cat("\nProportion of sites with at least one detection:\n", round(x$out.props[, "naive.occ"], digits), "\n\n")
    
    cat("Frequencies of sites with detections:\n")
    ##add matrix of frequencies
    print(x$out.freqs)
  }

  if(x$n.seasons > 1) {
      
      cat("\nSummary of counts across", x$n.seasons, "seasons:\n")
      count.mat <- matrix(x$count.table.full, nrow = 1)
      colnames(count.mat) <- names(x$count.table.full)
      rownames(count.mat) <- "Frequency"
      print(count.mat)

      cat("\nSummary of counts for each distance class across", x$n.seasons, "seasons:\n")
      dist.mat <- matrix(x$dist.sums.full, nrow = 1)
      colnames(dist.mat) <- names(x$dist.sums.full)
      rownames(dist.mat) <- "Frequency"
      print(dist.mat)
      cat("\n")
      
      cat("\nSeason-specific counts for each distance class: \n")
      cat("\n")
      for(i in 1:x$n.seasons) {
            cat("Season", i, "\n")
            temp.tab <- x$dist.table.seasons[[i]]
            out.mat <- matrix(temp.tab, nrow = 1)
            colnames(out.mat) <- names(temp.tab)
            rownames(out.mat) <- "Frequency"
            print(out.mat)
            cat("--------\n\n")
      }

      ##cat("\n")
      cat("Frequencies of sites with detections, extinctions, and colonizations:\n")
      ##add matrix of frequencies
      print(x$out.freqs)
  }
}

