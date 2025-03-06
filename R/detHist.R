##summarize detection histories and count data
detHist <- function(object, ...){
  UseMethod("detHist", object)
}


detHist.default <- function(object, ...){
  stop("\nFunction not yet defined for this object class\n")
}



##for unmarkedFrameOccu (same as data format for occuRN)
detHist.unmarkedFrameOccu <- function(object, ...) {

    ##extract data
    yMat <- object@y
    nsites <- nrow(yMat)
    n.seasons <- 1
    nvisits <- ncol(yMat)
    ##visits per season
    n.visits.season <- nvisits/n.seasons
    seasonNames <- paste("season", 1:n.seasons, sep = "")
    
    ##summarize detection histories
    hist.full <- apply(X = yMat, MARGIN = 1, FUN = function(i) paste(i, collapse = ""))
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

    out.det <- list("hist.table.full" = hist.table.full,
                    "hist.table.seasons" = hist.table.seasons,
                    "out.freqs" = out.freqs, "out.props" = out.props,
                    "n.seasons" = n.seasons,
                    "n.visits.season" = n.visits.season,
                    "n.species" = 1, "missing.seasons" = FALSE)
    class(out.det) <- "detHist"
    return(out.det)
}



##for occu
detHist.unmarkedFitOccu <- function(object, ...) {

    ##extract data
    yMat <- object@data@y
    nsites <- nrow(yMat)
    n.seasons <- 1
    nvisits <- ncol(yMat)
    ##visits per season
    n.visits.season <- nvisits/n.seasons
    seasonNames <- paste("season", 1:n.seasons, sep = "")
    
    ##summarize detection histories
    hist.full <- apply(X = yMat, MARGIN = 1, FUN = function(i) paste(i, collapse = ""))
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

    out.det <- list("hist.table.full" = hist.table.full,
                    "hist.table.seasons" = hist.table.seasons,
                    "out.freqs" = out.freqs, "out.props" = out.props,
                    "n.seasons" = n.seasons,
                    "n.visits.season" = n.visits.season,
                    "n.species" = 1, "missing.seasons" = FALSE)
  class(out.det) <- "detHist"
  return(out.det)
}



##for unmarkedFrameOccuFP
detHist.unmarkedFrameOccuFP <- function(object, ...) {

    ##extract data
    yMat <- object@y
    nsites <- nrow(yMat)
    n.seasons <- 1
    nvisits <- ncol(yMat)
    ##visits per season
    n.visits.season <- nvisits/n.seasons
    seasonNames <- paste("season", 1:n.seasons, sep = "")
    
    ##summarize detection histories
    hist.full <- apply(X = yMat, MARGIN = 1, FUN = function(i) paste(i, collapse = ""))
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

    out.det <- list("hist.table.full" = hist.table.full,
                    "hist.table.seasons" = hist.table.seasons,
                    "out.freqs" = out.freqs, "out.props" = out.props,
                    "n.seasons" = n.seasons,
                    "n.visits.season" = n.visits.season,
                    "n.species" = 1, "missing.seasons" = FALSE)
    class(out.det) <- "detHist"
    return(out.det)
}



##for occuFP
detHist.unmarkedFitOccuFP <- function(object, ...) {

    ##extract data
    yMat <- object@data@y
    nsites <- nrow(yMat)
    n.seasons <- 1
    nvisits <- ncol(yMat)
    ##visits per season
    n.visits.season <- nvisits/n.seasons
    seasonNames <- paste("season", 1:n.seasons, sep = "")
    
    ##summarize detection histories
    hist.full <- apply(X = yMat, MARGIN = 1, FUN = function(i) paste(i, collapse = ""))
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

    out.det <- list("hist.table.full" = hist.table.full,
                    "hist.table.seasons" = hist.table.seasons,
                    "out.freqs" = out.freqs, "out.props" = out.props,
                    "n.seasons" = n.seasons,
                    "n.visits.season" = n.visits.season,
                    "n.species" = 1, "missing.seasons" = FALSE)
  class(out.det) <- "detHist"
  return(out.det)
}



##for occuRN
detHist.unmarkedFitOccuRN <- function(object, ...) {

    ##extract data
    yMat <- object@data@y
    nsites <- nrow(yMat)
    n.seasons <- 1
    nvisits <- ncol(yMat)
    ##visits per season
    n.visits.season <- nvisits/n.seasons
    seasonNames <- paste("season", 1:n.seasons, sep = "")
    
    ##summarize detection histories
    hist.full <- apply(X = yMat, MARGIN = 1, FUN = function(i) paste(i, collapse = ""))
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
    
    out.det <- list("hist.table.full" = hist.table.full,
                    "hist.table.seasons" = hist.table.seasons,
                    "out.freqs" = out.freqs, "out.props" = out.props,
                    "n.seasons" = n.seasons,
                    "n.visits.season" = n.visits.season,
                    "n.species" = 1, "missing.seasons" = FALSE)
    class(out.det) <- "detHist"
    return(out.det)
}



##for unmarkedMultFrame
detHist.unmarkedMultFrame <- function(object, ...) {

    ##extract data
    yMat <- object@y
    nsites <- nrow(yMat)
    n.seasons <- object@numPrimary
    nvisits <- ncol(yMat)
    ##visits per season
    n.visits.season <- nvisits/n.seasons
    seasonNames <- paste("season", 1:n.seasons, sep = "")
    
    ##summarize detection histories
    ##starting and ending columns
    colStarts <- seq(from = 1, to = nvisits, by = n.visits.season)
    colEnds <- colStarts + (n.visits.season - 1)
    yrows <- list( )

    ##add check for seasons not sampled
    y.seasons <- list( )
    
    ##subsequent seasons
    for(i in 1:n.seasons) {
        yrows[[i]] <- apply(yMat[, colStarts[i]:colEnds[i]], MARGIN = 1, FUN = function(i) paste(i, collapse = ""))
        y.seasons[[i]] <- yMat[, colStarts[i]:colEnds[i]]
    }

    ##check if any seasons were not sampled
    y.seasonsNA <- sapply(y.seasons, FUN = function(i) all(is.na(i)))

    ##organize and paste rows
    hist.full <- do.call(what = "paste", args = c(yrows, sep = "'"))
    hist.table.full <- table(hist.full, deparse.level = 0)
    
    ##for each season, determine frequencies
    out.seasons <- vector(mode = "list", length = n.seasons)
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
        ##summarize detection histories
        det.hist <- apply(X = ySeason, MARGIN = 1, FUN = function(i) paste(i, collapse = ""))
        hist.table.seasons[[i]] <- table(det.hist, deparse.level = 0)
        
        ##determine proportion of sites with at least 1 detection
        det.sum <- apply(X = ySeason, MARGIN = 1, FUN = function(i) ifelse(sum(i, na.rm = TRUE) > 0, 1, 0))
        
        ##check sites with observed detections and deal with NA's
        sum.rows <- rowSums(ySeason, na.rm = TRUE)
        is.na(sum.rows) <- rowSums(is.na(ySeason)) == ncol(ySeason)
        
        ##number of sites sampled
        out.freqs[i, 1] <- sum(!is.na(sum.rows))
        ##detections
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
    
    out.det <- list("hist.table.full" = hist.table.full,
                    "hist.table.seasons" = hist.table.seasons,
                    "out.freqs" = out.freqs, "out.props" = out.props,
                    "n.seasons" = n.seasons,
                    "n.visits.season" = n.visits.season,
                    "n.species" = 1, "missing.seasons" = y.seasonsNA)
    class(out.det) <- "detHist"
    return(out.det)
}



##for colext
detHist.unmarkedFitColExt <- function(object, ...) {

    ##extract data
    yMat <- object@data@y
    nsites <- nrow(yMat)
    n.seasons <- object@data@numPrimary
    nvisits <- ncol(yMat)
    ##visits per season
    n.visits.season <- nvisits/n.seasons
    seasonNames <- paste("season", 1:n.seasons, sep = "")

    ##summarize detection histories
    ##starting and ending columns
    colStarts <- seq(from = 1, to = nvisits, by = n.visits.season)
    colEnds <- colStarts + (n.visits.season - 1)
    yrows <- list( )

    ##add check for seasons not sampled
    y.seasons <- list( )
       
    ##subsequent seasons
    for(i in 1:n.seasons) {
        yrows[[i]] <- apply(yMat[, colStarts[i]:colEnds[i]], MARGIN = 1, FUN = function(i) paste(i, collapse = ""))
        y.seasons[[i]] <- yMat[, colStarts[i]:colEnds[i]]
    }

    ##check if any seasons were not sampled
    y.seasonsNA <- sapply(y.seasons, FUN = function(i) all(is.na(i)))

    ##organize and paste rows
    hist.full <- do.call(what = "paste", args = c(yrows, sep = "'"))
    hist.table.full <- table(hist.full, deparse.level = 0)

    ##for each season, determine frequencies
    out.seasons <- vector(mode = "list", length = n.seasons)
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
        ##summarize detection histories
        det.hist <- apply(X = ySeason, MARGIN = 1, FUN = function(i) paste(i, collapse = ""))
        hist.table.seasons[[i]] <- table(det.hist, deparse.level = 0)

        ##determine proportion of sites with at least 1 detection
        det.sum <- apply(X = ySeason, MARGIN = 1, FUN = function(i) ifelse(sum(i, na.rm = TRUE) > 0, 1, 0))

        ##check sites with observed detections and deal with NA's
        sum.rows <- rowSums(ySeason, na.rm = TRUE)
        is.na(sum.rows) <- rowSums(is.na(ySeason)) == ncol(ySeason)
        
        ##number of sites sampled
        out.freqs[i, 1] <- sum(!is.na(sum.rows))
        ##detections
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
        out.det <- list("hist.table.full" = hist.table.full,
                        "hist.table.seasons" = hist.table.seasons,
                        "out.freqs" = out.freqs, "out.props" = out.props,
                        "n.seasons" = n.seasons,
                        "n.visits.season" = n.visits.season,
                        "n.species" = 1, "missing.seasons" = y.seasonsNA)
        class(out.det) <- "detHist"
        return(out.det)
}



#############################
#############################
##TO CHANGE HERE

##for unmarkedFrameOccuMulti
detHist.unmarkedFrameOccuMulti <- function(object, ...) {

    ##extract species detection data
    speciesList <- object@ylist
    speciesNames <- names(object@ylist)
    if(is.null(speciesNames)) {
        speciesNames <- paste("species", 1:nspecies, sep = "")
    }
    nspecies <- length(speciesList)
    n.seasons <- 1
    nsites <- nrow(speciesList[[1]])
    nvisits <- ncol(speciesList[[1]])

    ##visits per season
    n.visits.season <- nvisits/n.seasons

    ##generic name to include in detection history
    genericNames <- letters[1:nspecies]

    ##combine detection histories of each species
    histList <- vector(mode = "list", length = nspecies)
    for(sp in 1:nspecies) {
        detVector <- as.vector(speciesList[[sp]])
        histList[[sp]] <- ifelse(detVector == 1, genericNames[sp], detVector)
    }

    comboDet <- do.call("paste", c(histList, sep = ""))
    ##number of co-occurrences in any given survey across sites
    coOcc <- table(comboDet)
    
    comboMat <- matrix(comboDet, nrow = nsites, ncol = nvisits)
    ##detection histories
    comboHist <- apply(comboMat, MARGIN = 1, FUN = function(i) paste(i, collapse = "'"))
    hist.table.full <- table(comboHist)
    
    ##for each season, determine frequencies
    out.freqs <- matrix(data = NA, ncol = 2, nrow = nspecies)
    colnames(out.freqs) <- c("sampled", "detected")
    rownames(out.freqs) <- speciesNames

    ##create a matrix with proportion of sites with detections
    ##based on raw data
    out.props <- matrix(NA, nrow = nrow(out.freqs), ncol = 1)
    colnames(out.props) <- "naive.occ"
    rownames(out.props) <- rownames(out.freqs)
        
    hist.table.species <- vector(mode = "list", length = nspecies)
    names(hist.table.species) <- speciesNames
    for(i in 1:nspecies) {

        yMat <- speciesList[[i]]
        
        ##summarize detection histories
        hist.full <- apply(X = yMat, MARGIN = 1, FUN = function(i) paste(i, collapse = ""))
        hist.table.species[[i]] <- table(hist.full, deparse.level = 0)
        
        ##determine proportion of sites with at least 1 detection
        det.sum <- apply(X = speciesList[[i]], MARGIN = 1, FUN = function(i) ifelse(sum(i, na.rm = TRUE) > 0, 1, 0))

        ##check sites with observed detections and deal with NA's
        sum.rows <- rowSums(speciesList[[i]], na.rm = TRUE)
        is.na(sum.rows) <- rowSums(is.na(speciesList[[i]])) == ncol(yMat)

        ##number of sites sampled
        out.freqs[i, 1] <- sum(!is.na(sum.rows))
        ##number of sites with at least 1 detection
        out.freqs[i, 2] <- sum(det.sum)

        ##proportion of sites with detections
        out.props[i, 1] <- out.freqs[i, 2]/out.freqs[i, 1]
        
    }

    ##add frequencies of co-occurrences
    hist.table.species$coOcc <- coOcc
    
    out.det <- list("hist.table.full" = hist.table.full,
                    "hist.table.species" = hist.table.species,
                    "out.freqs" = out.freqs, "out.props" = out.props,
                    "n.seasons" = n.seasons,
                    "n.visits.season" = n.visits.season,
                    "n.species" = nspecies, "missing.seasons" = FALSE)
    class(out.det) <- "detHist"
    return(out.det)
}



##for occuMulti
detHist.unmarkedFitOccuMulti <- function(object, ...) {

    ##extract species detection data
    speciesList <- object@data@ylist
    speciesNames <- names(object@data@ylist)
    if(is.null(speciesNames)) {
        speciesNames <- paste("species", 1:nspecies, sep = "")
    }
    nspecies <- length(speciesList)
    n.seasons <- 1
    nsites <- nrow(speciesList[[1]])
    nvisits <- ncol(speciesList[[1]])

    ##visits per season
    n.visits.season <- nvisits/n.seasons

    ##generic name to include in detection history
    genericNames <- letters[1:nspecies]

    ##combine detection histories of each species
    histList <- vector(mode = "list", length = nspecies)
    for(sp in 1:nspecies) {
        detVector <- as.vector(speciesList[[sp]])
        histList[[sp]] <- ifelse(detVector == 1, genericNames[sp], detVector)
    }

    comboDet <- do.call("paste", c(histList, sep = ""))
    ##number of co-occurrences in any given survey across sites
    coOcc <- table(comboDet)
    
    comboMat <- matrix(comboDet, nrow = nsites, ncol = nvisits)
    ##detection histories
    comboHist <- apply(comboMat, MARGIN = 1, FUN = function(i) paste(i, collapse = "'"))
    hist.table.full <- table(comboHist)
    
    ##for each season, determine frequencies
    out.freqs <- matrix(data = NA, ncol = 2, nrow = nspecies)
    colnames(out.freqs) <- c("sampled", "detected")
    rownames(out.freqs) <- speciesNames

    ##create a matrix with proportion of sites with detections
    ##based on raw data
    out.props <- matrix(NA, nrow = nrow(out.freqs), ncol = 1)
    colnames(out.props) <- "naive.occ"
    rownames(out.props) <- rownames(out.freqs)
        
    hist.table.species <- vector(mode = "list", length = nspecies)
    names(hist.table.species) <- speciesNames
    
    for(i in 1:nspecies) {

        yMat <- speciesList[[i]]
        
        ##summarize detection histories
        hist.full <- apply(X = yMat, MARGIN = 1, FUN = function(i) paste(i, collapse = ""))
        hist.table.species[[i]] <- table(hist.full, deparse.level = 0)
        
        ##determine proportion of sites with at least 1 detection
        det.sum <- apply(X = speciesList[[i]], MARGIN = 1, FUN = function(i) ifelse(sum(i, na.rm = TRUE) > 0, 1, 0))

        ##check sites with observed detections and deal with NA's
        sum.rows <- rowSums(speciesList[[i]], na.rm = TRUE)
        is.na(sum.rows) <- rowSums(is.na(speciesList[[i]])) == ncol(yMat)

        ##number of sites sampled
        out.freqs[i, 1] <- sum(!is.na(sum.rows))
        ##number of sites with at least 1 detection
        out.freqs[i, 2] <- sum(det.sum)

        ##proportion of sites with detections
        out.props[i, 1] <- out.freqs[i, 2]/out.freqs[i, 1]
        
    }

    ##add frequencies of co-occurrences
    hist.table.species$coOcc <- coOcc
    
    out.det <- list("hist.table.full" = NULL,
                    "hist.table.species" = hist.table.species,
                    "out.freqs" = out.freqs, "out.props" = out.props,
                    "n.seasons" = n.seasons,
                    "n.visits.season" = n.visits.season,
                    "n.species" = nspecies, "missing.seasons" = FALSE)
    class(out.det) <- "detHist"
    return(out.det)
}



##for occuMS
detHist.unmarkedFrameOccuMS <- function(object, ...) {

    ##extract data
    yMat <- object@y
    nsites <- nrow(yMat)
    n.seasons <- object@numPrimary
    nvisits <- ncol(yMat)
    ##visits per season
    n.visits.season <- nvisits/n.seasons
    seasonNames <- paste("season", 1:n.seasons, sep = "")
    ##no missing season when single season
    y.seasonsNA <- FALSE
        
    ##for each season, determine frequencies
    out.seasons <- vector(mode = "list", length = n.seasons)
    hist.table.seasons <- vector(mode = "list", length = n.seasons)
    names(hist.table.seasons) <- seasonNames

    if(n.seasons == 1) {
        ##summarize detection histories
        hist.full <- apply(X = yMat, MARGIN = 1, FUN = function(i) paste(i, collapse = ""))
        hist.table.full <- table(hist.full, deparse.level = 0)

        out.freqs <- matrix(data = NA, ncol = 2, nrow = 1)
        colnames(out.freqs) <- c("sampled", "detected")
        rownames(out.freqs) <- "Season-1"

        ##sequence of visits
        vis.seq <- seq(from = 1, to = nvisits, by = n.visits.season)
        for(i in 1:n.seasons) {
            col.start <- vis.seq[i]
            col.end <- col.start + (n.visits.season - 1)
            ySeason <- yMat[, col.start:col.end]
            ##summarize detection histories
            det.hist <- apply(X = ySeason, MARGIN = 1, FUN = function(i) paste(i, collapse = ""))
            hist.table.seasons[[i]] <- table(det.hist, deparse.level = 0)

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

        out.props <- matrix(NA, nrow = 1, ncol = 1)
        colnames(out.props) <- "naive.occ"
        rownames(out.props) <- rownames(out.freqs)
        out.props[, 1] <- out.freqs[, 2]/out.freqs[, 1]
    }

    if(n.seasons > 1) {

        ##summarize detection histories
        ##starting and ending columns
        colStarts <- seq(from = 1, to = nvisits, by = n.visits.season)
        colEnds <- colStarts + (n.visits.season - 1)
        yrows <- list( )
        yMat.seasons <- vector(mode = "list", length = n.seasons)
        
        ##subsequent seasons
        for(i in 1:n.seasons) {
            yrows[[i]] <- apply(yMat[, colStarts[i]:colEnds[i]], MARGIN = 1, FUN = function(i) paste(i, collapse = ""))
            
            yMat.seasons[[i]] <- yMat[, colStarts[i]:colEnds[i]]
            
        }

        ##check if any seasons were not sampled
        y.seasonsNA <- sapply(yMat.seasons, FUN = function(i) all(is.na(i)))
                
        ##organize and paste rows
        hist.full <- do.call(what = "paste", args = c(yrows, sep = "'"))
        hist.table.full <- table(hist.full, deparse.level = 0)

        
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
            ##summarize detection histories
            det.hist <- apply(X = ySeason, MARGIN = 1, FUN = function(i) paste(i, collapse = ""))
            hist.table.seasons[[i]] <- table(det.hist, deparse.level = 0)

            ##determine number of sites with at least 1 detection
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
    
    out.det <- list("hist.table.full" = hist.table.full,
                    "hist.table.seasons" = hist.table.seasons,
                    "out.freqs" = out.freqs, "out.props" = out.props,
                    "n.seasons" = n.seasons,
                    "n.visits.season" = n.visits.season,
                    "n.species" = 1,
                    "missing.seasons" = y.seasonsNA)
    class(out.det) <- "detHist"
    return(out.det)
}



##for occuMS
detHist.unmarkedFitOccuMS <- function(object, ...) {


    ##extract data
    yMat <- object@data@y
    nsites <- nrow(yMat)
    n.seasons <- object@data@numPrimary
    nvisits <- ncol(yMat)
    ##visits per season
    n.visits.season <- nvisits/n.seasons
    seasonNames <- paste("season", 1:n.seasons, sep = "")
    ##no missing season when single season
    y.seasonsNA <- FALSE
    
    ##for each season, determine frequencies
    out.seasons <- vector(mode = "list", length = n.seasons)
    hist.table.seasons <- vector(mode = "list", length = n.seasons)
    names(hist.table.seasons) <- seasonNames
    
    if(n.seasons == 1) {
        ##summarize detection histories
        hist.full <- apply(X = yMat, MARGIN = 1, FUN = function(i) paste(i, collapse = ""))
        hist.table.full <- table(hist.full, deparse.level = 0)

        out.freqs <- matrix(data = NA, ncol = 2, nrow = 1)
        colnames(out.freqs) <- c("sampled", "detected")
        rownames(out.freqs) <- "Season-1"

        ##sequence of visits
        vis.seq <- seq(from = 1, to = nvisits, by = n.visits.season)
        for(i in 1:n.seasons) {
            col.start <- vis.seq[i]
            col.end <- col.start + (n.visits.season - 1)
            ySeason <- yMat[, col.start:col.end]
            ##summarize detection histories
            det.hist <- apply(X = ySeason, MARGIN = 1, FUN = function(i) paste(i, collapse = ""))
            hist.table.seasons[[i]] <- table(det.hist, deparse.level = 0)

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

        out.props <- matrix(NA, nrow = 1, ncol = 1)
        colnames(out.props) <- "naive.occ"
        rownames(out.props) <- rownames(out.freqs)
        out.props[, 1] <- out.freqs[, 2]/out.freqs[, 1]
    }

    if(n.seasons > 1) {

        ##summarize detection histories
        ##starting and ending columns
        colStarts <- seq(from = 1, to = nvisits, by = n.visits.season)
        colEnds <- colStarts + (n.visits.season - 1)
        yrows <- list( )
        yMat.seasons <- vector(mode = "list", length = n.seasons)
        
        ##subsequent seasons
        for(i in 1:n.seasons) {
            yrows[[i]] <- apply(yMat[, colStarts[i]:colEnds[i]], MARGIN = 1, FUN = function(i) paste(i, collapse = ""))

            yMat.seasons[[i]] <- yMat[, colStarts[i]:colEnds[i]]
        }

        ##check if any seasons were not sampled
        y.seasonsNA <- sapply(yMat.seasons, FUN = function(i) all(is.na(i)))
        
        ##organize and paste rows
        hist.full <- do.call(what = "paste", args = c(yrows, sep = "'"))
        hist.table.full <- table(hist.full, deparse.level = 0)

        
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
            ##summarize detection histories
            det.hist <- apply(X = ySeason, MARGIN = 1, FUN = function(i) paste(i, collapse = ""))
            hist.table.seasons[[i]] <- table(det.hist, deparse.level = 0)

            ##determine number of sites with at least 1 detection
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
    

    out.det <- list("hist.table.full" = hist.table.full,
                    "hist.table.seasons" = hist.table.seasons,
                    "out.freqs" = out.freqs, "out.props" = out.props,
                    "n.seasons" = n.seasons,
                    "n.visits.season" = n.visits.season,
                    "n.species" = 1,
                    "missing.seasons" = y.seasonsNA)
    class(out.det) <- "detHist"
    return(out.det)
}



##ISSUE WITH DATA TYPE:
##goccu uses unmarkedMultFrame format (same as colext)

##currently returns list of detection histories pooled across devices and list of detection histories specific to each device (scale)
##goccu
detHist.unmarkedFitGOccu <- function(object, ...) {

    ##extract data
    yMat <- object@data@y
    nsites <- nrow(yMat)
    n.seasons <- object@data@numPrimary
    ##number of seasons (corrected)
    n.seasons.adj <- 1
    nvisits <- ncol(yMat)
    ##visits per season
    n.visits.season <- nvisits/n.seasons
    seasonNames <- paste("season", 1:n.seasons, sep = "")

    ##summarize detection histories
    ##starting and ending columns
    colStarts <- seq(from = 1, to = nvisits, by = n.visits.season)
    colEnds <- colStarts + (n.visits.season - 1)
    yrows <- list( )

    ##add check for seasons not sampled
    y.seasons <- list( )
       
    ##subsequent seasons
    for(i in 1:n.seasons) {
        yrows[[i]] <- apply(yMat[, colStarts[i]:colEnds[i]], MARGIN = 1, FUN = function(i) paste(i, collapse = ""))
        y.seasons[[i]] <- yMat[, colStarts[i]:colEnds[i]]
    }

    ##check if any seasons were not sampled
    y.seasonsNA <- sapply(y.seasons, FUN = function(i) all(is.na(i)))


###################
#################
    ##compute detections across secondary periods
    hist.reduced <- matrix(data = NA,
                           nrow = nsites,
                           ncol = n.seasons)
    for(i in 1:n.seasons) {
        hist.reduced[, i] <- apply(y.seasons[[i]], 1, FUN = function(x) ifelse(sum(x) > 0, 1, 0))
    }

    ##summarize detection histories
    hist.full <- apply(X = hist.reduced, MARGIN = 1, FUN = function(i) paste(i, collapse = ""))
    hist.table.full <- table(hist.full, deparse.level = 0)

    ##CODE CHANGED HERE FOR TABLE MODIFIED TO DISPLAY A SINGLE ROW
    ##for each season, determine frequencies
    ##out.seasons <- vector(mode = "list", length = n.seasons)
    hist.table.seasons <- vector(mode = "list", length = n.seasons)
    names(hist.table.seasons) <- seasonNames

        
    ##CODE CHANGED HERE FOR TABLE MODIFIED TO DISPLAY A SINGLE ROW
    out.freqs <- matrix(data = NA, ncol = 2, nrow = n.seasons)
    colnames(out.freqs) <- c("sampled", "detected")
    rownames(out.freqs) <- paste("Season-", 1:n.seasons, sep = "")

    ##sequence of visits
    vis.seq <- seq(from = 1, to = nvisits, by = n.visits.season)
    for(i in 1:n.seasons) {
        col.start <- vis.seq[i]
        col.end <- col.start + (n.visits.season - 1)
        ySeason <- yMat[, col.start:col.end]
        ##summarize detection histories
        det.hist <- apply(X = ySeason, MARGIN = 1, FUN = function(i) paste(i, collapse = ""))
        hist.table.seasons[[i]] <- table(det.hist, deparse.level = 0)

        ##determine proportion of sites with at least 1 detection
        det.sum <- apply(X = ySeason, MARGIN = 1, FUN = function(i) ifelse(sum(i, na.rm = TRUE) > 0, 1, 0))

        ##check sites with observed detections and deal with NA's
        sum.rows <- rowSums(ySeason, na.rm = TRUE)
        is.na(sum.rows) <- rowSums(is.na(ySeason)) == ncol(ySeason)
        
        ##number of sites sampled
        out.freqs[i, 1] <- sum(!is.na(sum.rows))
        ##detections
        out.freqs[i, 2] <- sum(det.sum)
        
     }

    ##create a matrix with proportion of sites with colonizations
    ##and extinctions based on raw data
    out.props <- matrix(NA, nrow = nrow(out.freqs), ncol = 1)
    colnames(out.props) <- "naive.occ"
    rownames(out.props) <- rownames(out.freqs)
    out.props[, 1] <- out.freqs[, 2]/out.freqs[, 1]
    
    out.det <- list("hist.table.full" = hist.table.full,
                    "hist.table.seasons" = hist.table.seasons,
                    "out.freqs" = out.freqs, "out.props" = out.props,
                    "n.seasons" = n.seasons,
                    "n.visits.season" = n.visits.season,
                    "n.species" = 1, "missing.seasons" = y.seasonsNA)
    class(out.det) <- "detHist"
    return(out.det)
}



##for unmarkedFrameOccuComm
detHist.unmarkedFrameOccuComm <- function(object, ...) {

    ##extract species detection data
    speciesList <- object@ylist
    speciesNames <- names(object@ylist)
    if(is.null(speciesNames)) {
        speciesNames <- paste("species", 1:nspecies, sep = "")
    }
    nspecies <- length(speciesList)
    n.seasons <- 1
    nsites <- nrow(speciesList[[1]])
    nvisits <- ncol(speciesList[[1]])

    ##visits per season
    n.visits.season <- nvisits/n.seasons

    ##generic name to include in detection history
    genericNames <- 1:nspecies

    ##for each season, determine frequencies
    out.freqs <- matrix(data = NA, ncol = 2, nrow = nspecies)
    colnames(out.freqs) <- c("sampled", "detected")
    rownames(out.freqs) <- speciesNames

    ##create a matrix with proportion of sites with detections
    ##based on raw data
    out.props <- matrix(NA, nrow = nrow(out.freqs), ncol = 1)
    colnames(out.props) <- "naive.occ"
    rownames(out.props) <- rownames(out.freqs)
        
    hist.table.species <- vector(mode = "list", length = nspecies)
    names(hist.table.species) <- speciesNames
    for(i in 1:nspecies) {

        yMat <- speciesList[[i]]
        
        ##summarize detection histories
        hist.full <- apply(X = yMat, MARGIN = 1, FUN = function(i) paste(i, collapse = ""))
        hist.table.species[[i]] <- table(hist.full, deparse.level = 0)
        
        ##determine proportion of sites with at least 1 detection
        det.sum <- apply(X = speciesList[[i]], MARGIN = 1, FUN = function(i) ifelse(sum(i, na.rm = TRUE) > 0, 1, 0))

        ##check sites with observed detections and deal with NA's
        sum.rows <- rowSums(speciesList[[i]], na.rm = TRUE)
        is.na(sum.rows) <- rowSums(is.na(speciesList[[i]])) == ncol(yMat)

        ##number of sites sampled
        out.freqs[i, 1] <- sum(!is.na(sum.rows))
        ##number of sites with at least 1 detection
        out.freqs[i, 2] <- sum(det.sum)

        ##proportion of sites with detections
        out.props[i, 1] <- out.freqs[i, 2]/out.freqs[i, 1]
        
    }

    ##add frequencies of co-occurrences
    ##hist.table.species$coOcc <- coOcc
    
    out.det <- list("hist.table.full" = NULL,
                    "hist.table.species" = hist.table.species,
                    "out.freqs" = out.freqs, "out.props" = out.props,
                    "n.seasons" = n.seasons,
                    "n.visits.season" = n.visits.season,
                    "n.species" = nspecies, "missing.seasons" = FALSE)
    class(out.det) <- "detHist"
    return(out.det)
}



##for occuComm
detHist.unmarkedFitOccuComm <- function(object, ...) {

    ##extract species detection data
    speciesList <- object@data@ylist
    speciesNames <- names(object@data@ylist)
    if(is.null(speciesNames)) {
        speciesNames <- paste("species", 1:nspecies, sep = "")
    }
    nspecies <- length(speciesList)
    n.seasons <- 1
    nsites <- nrow(speciesList[[1]])
    nvisits <- ncol(speciesList[[1]])

    ##visits per season
    n.visits.season <- nvisits/n.seasons

    ##generic name to include in detection history
    genericNames <- 1:nspecies

    ##for each season, determine frequencies
    out.freqs <- matrix(data = NA, ncol = 2, nrow = nspecies)
    colnames(out.freqs) <- c("sampled", "detected")
    rownames(out.freqs) <- speciesNames

    ##create a matrix with proportion of sites with detections
    ##based on raw data
    out.props <- matrix(NA, nrow = nrow(out.freqs), ncol = 1)
    colnames(out.props) <- "naive.occ"
    rownames(out.props) <- rownames(out.freqs)
        
    hist.table.species <- vector(mode = "list", length = nspecies)
    names(hist.table.species) <- speciesNames
    for(i in 1:nspecies) {

        yMat <- speciesList[[i]]
        
        ##summarize detection histories
        hist.full <- apply(X = yMat, MARGIN = 1, FUN = function(i) paste(i, collapse = ""))
        hist.table.species[[i]] <- table(hist.full, deparse.level = 0)
        
        ##determine proportion of sites with at least 1 detection
        det.sum <- apply(X = speciesList[[i]], MARGIN = 1, FUN = function(i) ifelse(sum(i, na.rm = TRUE) > 0, 1, 0))

        ##check sites with observed detections and deal with NA's
        sum.rows <- rowSums(speciesList[[i]], na.rm = TRUE)
        is.na(sum.rows) <- rowSums(is.na(speciesList[[i]])) == ncol(yMat)

        ##number of sites sampled
        out.freqs[i, 1] <- sum(!is.na(sum.rows))
        ##number of sites with at least 1 detection
        out.freqs[i, 2] <- sum(det.sum)

        ##proportion of sites with detections
        out.props[i, 1] <- out.freqs[i, 2]/out.freqs[i, 1]
        
    }

    ##add frequencies of co-occurrences
    ##hist.table.species$coOcc <- coOcc
    
    out.det <- list("hist.table.full" = NULL,
                    "hist.table.species" = hist.table.species,
                    "out.freqs" = out.freqs, "out.props" = out.props,
                    "n.seasons" = n.seasons,
                    "n.visits.season" = n.visits.season,
                    "n.species" = nspecies, "missing.seasons" = FALSE)
    class(out.det) <- "detHist"
    return(out.det)
}



##print method
print.detHist <- function(x, digits = 2, ...) {
    ##convert NA to . for nicer printing
    hist.names <- names(x$hist.table.full)
    if(!is.null(x$hist.table.full)) {
        names(x$hist.table.full) <- gsub(pattern = "NA",
                                         replacement = ".",
                                         x = hist.names)
    }
    if(identical(x$n.seasons, 1)) {
        if(x$n.species > 1) {
            nspecies <- x$n.species
            speciesNames <- rownames(x$out.freqs)

            ##convert NA to . for nicer printing
            for(d in 1:nspecies) {
                hist.names <- names(x$hist.table.species[[d]])
                names(x$hist.table.species[[d]]) <- gsub(pattern = "NA",
                                                         replacement = ".",
                                                         x = hist.names)
            }
            
            if(!is.null(x$hist.table.full)) {
                ##species code in detection histories
                speciesCode <- character( )
                for(j in 1:nspecies) {
                    speciesCode[j] <- paste(speciesNames[j], " (", letters[j], ")", sep = "")
                }
                cat("\nSummary of detection histories: \n")
                num.chars <- nchar(paste(names(x$hist.table.full), collapse = ""))
                if(num.chars >= 80) {
                    cat("\nNote:  Detection histories exceed 80 characters and are not displayed\n")
                } else {
                    cat("(")
                    cat(speciesCode, sep = ", ")
                    cat(")\n")
                
                    out.mat <- matrix(x$hist.table.full, nrow = 1)
                    colnames(out.mat) <- names(x$hist.table.full)
                    rownames(out.mat) <- "Frequency"
                    print(out.mat)
                }
            }

            cat("\nSpecies-specific detection histories: \n")
            cat("\n")
            for(i in 1:nspecies) {
                cat(speciesNames[i], "\n")
                temp.tab <- x$hist.table.species[[i]]
                out.mat <- matrix(temp.tab, nrow = 1)
                colnames(out.mat) <- names(temp.tab)
                rownames(out.mat) <- "Frequency"
                print(out.mat)
                cat("--------\n\n")
            }
            
            if(!is.null(x$hist.table.full)) {
                cat("Frequency of co-occurrence among sites: \n")
                cat("(")
                cat(speciesCode, sep = ", ")
                cat(")\n")
            
                occ.tab <- x$hist.table.species$coOcc
                occ.mat <- matrix(occ.tab, nrow = 1)
                colnames(occ.mat) <- names(occ.tab)
                rownames(occ.mat) <- "Frequency"
                print(occ.mat)
            }
            
            cat("\nProportion of sites with at least one detection:\n")
            print(x$out.props[, "naive.occ", drop = FALSE], digits)
            cat("\n")

            cat("Frequencies of sites with detections:\n")
            ##add matrix of frequencies
            print(x$out.freqs)
            
        } else {
        
            cat("\nSummary of detection histories: \n")
            out.mat <- matrix(x$hist.table.full, nrow = 1)
            colnames(out.mat) <- names(x$hist.table.full)
            rownames(out.mat) <- "Frequency"
            print(out.mat)
            cat("\nProportion of sites with at least one detection:\n", round(x$out.props[, "naive.occ"], digits), "\n\n")

            cat("Frequencies of sites with detections:\n")
            ##add matrix of frequencies
            print(x$out.freqs)

        }
    }

    ##multiseason model
    if(x$n.seasons > 1 && ncol(x$out.freqs) > 2) {
        ##convert NA to . for nicer printing
        for(d in 1:x$n.seasons) {
            hist.names <- names(x$hist.table.seasons[[d]])
            names(x$hist.table.seasons[[d]]) <- gsub(pattern = "NA",
                                                                replacement = ".",
                                                                x = hist.names)
        }

        cat("\nSummary of detection histories (", x$n.seasons, " seasons combined): \n", sep ="")
        ##determine number of characters
        num.chars <- nchar(paste(names(x$hist.table.full), collapse = ""))
        if(num.chars >= 80) {
            cat("\nNote:  Detection histories exceed 80 characters and are not displayed\n")
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
                
            cat("\nSeason-specific detection histories: \n")
            cat("\n")
            for(i in 1:x$n.seasons) {
                if(!x$missing.seasons[i]) {
                    cat("Season", i, "\n")
                } else {
                    cat("Season", i, "(no sites sampled)", "\n")
                }
                
                temp.tab <- x$hist.table.seasons[[i]]
                out.mat <- matrix(temp.tab, nrow = 1)
                colnames(out.mat) <- names(temp.tab)
                rownames(out.mat) <- "Frequency"
                print(out.mat)
                cat("--------\n\n")
            }
        }

        ##cat("\n")
        cat("Frequencies of sites with detections, extinctions, and colonizations:\n")
        ##add matrix of frequencies
        print(x$out.freqs)
    }


    ##goccu multiscale model
    if(x$n.seasons > 1 && ncol(x$out.freqs) == 2) {

        ##modify names
        orig.names <- names(x$hist.table.seasons)
        ##modify to "device"
        new.names <- gsub(pattern = "season", replacement = "device",
                          x = orig.names)
        ##convert NA to . for nicer printing
        for(d in 1:x$n.seasons) {
            hist.names <- names(x$hist.table.seasons[[d]])
            names(x$hist.table.seasons[[d]]) <- gsub(pattern = "NA",
                                                     replacement = ".",
                                                     x = hist.names)
        }

        cat("\nSummary of detection histories across scales (", x$n.seasons, " devices combined): \n", sep ="")
        ##determine number of characters
        num.chars <- nchar(paste(names(x$hist.table.full), collapse = ""))
        if(num.chars >= 80) {
            cat("\nNote:  Detection histories exceed 80 characters and are not displayed\n")
        } else {
            out.mat <- matrix(x$hist.table.full, nrow = 1)
            colnames(out.mat) <- names(x$hist.table.full)
            rownames(out.mat) <- "Frequency"
            print(out.mat)
        }

        ##if some seasons have not been sampled
        if(any(x$missing.seasons)) {
            if(sum(x$missing.seasons) == 1) {
                cat("\nNote: device", which(x$missing.seasons), "did not sample\n")
            } else {
                cat("\nNote: devices",
                    paste(which(x$missing.seasons), sep = ", "),
                    "did not sample\n")
            }
        }
        
                
        cat("\nDevice-specific detection histories: \n")
        cat("\n")
        for(i in 1:x$n.seasons) {
            if(!x$missing.seasons[i]) {
                cat("Device", i, "\n")
            } else {
                cat("Device", i, "(no sites sampled)", "\n")
            }
                
            temp.tab <- x$hist.table.seasons[[i]]
            out.mat <- matrix(temp.tab, nrow = 1)
            colnames(out.mat) <- names(temp.tab)
            rownames(out.mat) <- "Frequency"
            print(out.mat)
            cat("--------\n\n")
        }
    
        ##cat("\n")
        cat("Frequencies of sites with detections for each device:\n")
        ##add matrix of frequencies
        print(x$out.freqs)
    }
}
