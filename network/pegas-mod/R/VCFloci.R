## VCFloci.R (2019-11-15)

##   Handling VCF Files

## Copyright 2015-2019 Emmanuel Paradis

## This file is part of the R-package `pegas'.
## See the file ../DESCRIPTION for licensing issues.

.VCFconnection <- function(file)
{
    file <- path.expand(file) # fix by Frederic Michaud
    remote <- if (length(grep("^(ht|f)tp(s|):", file))) TRUE else FALSE
    GZ <- if (length(grep("\\.gz$", file))) TRUE else FALSE
    if (GZ) {
        file <- if (remote) url(file) else gzfile(file)
        file <- gzcon(file)
    }
    x <- readChar(file, 16L, TRUE)
    if (!identical(x, "##fileformat=VCF"))
        stop("file apparently not in VCF format")
    file
}

.getMETAvcf <- function(file)
{
    f <- .VCFconnection(file)
    HEADER <- character(0)
    skip <- 0L
    repeat {
        x <- scan(file, "", sep = "\n", n = 1L, skip = skip, quiet = TRUE)
        if (substr(x, 1, 6) == "#CHROM") break
        skip <- skip + 1L
        HEADER <- c(HEADER, x)
    }
    HEADER <- paste(paste0(HEADER, "\n"), collapse = "")
    j <- nchar(HEADER, "bytes") + 1L + nchar(x, "bytes")
    list(HEADER = HEADER, LABELS = x, position = j)
}

VCFheader <- function(file) .getMETAvcf(file)$HEADER

VCFlabels <- function(file) strsplit(.getMETAvcf(file)$LABELS, "\t")[[1]][-(1:9)]

VCFloci <- function(file, what = "all", chunck.size = 1e9, quiet = FALSE)
{
    meta <- .getMETAvcf(file)
    f <- .VCFconnection(file)
    GZ <- if (inherits(f, "connection")) TRUE else FALSE

    FIELDS <- c("CHROM", "POS", "ID", "REF", "ALT",
                "QUAL", "FILTER", "INFO", "FORMAT")

    what <-
        if (identical(what, "all")) 1:9 else match(what, FIELDS)

    obj <- vector("list", 9L)

    if (GZ) open(f) else {
        sz <- file.info(file)$size
        if (is.na(sz))
            stop(paste("cannot find information on file", sQuote(file)))
        left.to.scan <- sz
    }

    if (!quiet) cat("Scanning file", file, "\n")
    scanned <- 0 # NOT integer
    ncycle <- 0L
    FROM <- integer()
    TO <- integer()
    CHUNCK.SIZES <- numeric()

    repeat {
        if (!GZ) {
            if (left.to.scan > chunck.size) {
                left.to.scan <- left.to.scan - chunck.size
            } else {
                chunck.size <- left.to.scan
                left.to.scan <- 0L
            }
            Y <- .Call(read_bin_pegas, file, chunck.size, scanned)
        } else {
            Y <- readBin(f, "raw", chunck.size)
            if (!length(Y)) break
        }

        nY <- length(Y)
        scanned <- scanned + nY

        if (!quiet) {
            if (GZ) cat("\r", scanned/1e6, "Mb")
            else cat("\r", scanned/1e6, "/", sz/1e6, "Mb")
        }

        ncycle <- ncycle + 1L

        if (ncycle == 1) {
            skip <- meta$position - 3L
            nCol <- length(gregexpr("\t", meta$LABELS)[[1]]) + 1L
        } else skip <- 0L

        hop <- 2L * nCol - 1L
        EOL <- .Call(findEOL_C, Y, skip, hop) # can multiply 'hop' by 2 if diploid

        ck <- nY
        if (exists("trail", inherits = FALSE)) {
            x <- c(trail, Y[1:EOL[1L]])
            for (i in what) {
                tmp <-
                    if (i %in% c(2, 6)) .Call(extract_POS, x, c(0L, length(trail)), i - 1L)
                    else .Call(extract_REF, x, c(0L, length(trail)), i - 1L)
                obj[[i]] <- c(obj[[i]], tmp)
            }
            ck <- ck + length(trail)
            rm(trail)
            extra.locus <- 1L
        } else extra.locus <- 0L

        nEOL <- length(EOL)
        if (EOL[nEOL] != nY) {
            trail <- Y[(EOL[nEOL] + 1L):nY]
            ck <- ck - length(trail)
        }

        from <- if (ncycle == 1) 1L else TO[ncycle - 1L] + 1L
        to <- from + nEOL - 2L + extra.locus
        FROM <- c(FROM, from)
        TO <- c(TO, to)
        CHUNCK.SIZES <- c(CHUNCK.SIZES, ck)

        ## we assume there are no extra blank line!

        for (i in what) {
            tmp <-
                if (i %in% c(2, 6)) .Call(extract_POS, Y, EOL, i - 1L)
                else .Call(extract_REF, Y, EOL, i - 1L)
            obj[[i]] <- c(obj[[i]], tmp)
        }

        if (!GZ && !left.to.scan) break
    }

    if (GZ) close(f)
    else if (!quiet) cat("\r", scanned/1e6, "/", sz/1e6, "Mb")
    if (!quiet) cat("\nDone.\n")

    assign(file, data.frame(FROM = FROM, TO = TO, CHUNCK.SIZES = CHUNCK.SIZES),
           envir = .cacheVCF)

    names(obj) <- FIELDS
    obj <- obj[!sapply(obj, is.null)]
    obj <- as.data.frame(obj, stringsAsFactors = FALSE)
    class(obj) <- c("VCFinfo", "data.frame")
    obj
}

print.VCFinfo <- function(x, ...)
{
    n <- length(x[[1]])
    if (n < 10) print(as.data.frame(x)) else {
        x <- x[c(1:5, (n - 4):n), , drop = FALSE]
        x <- apply(x, 2, as.character)
        x <- as.data.frame(x, stringsAsFactors = FALSE)
        x[5:6, ] <- ""
        row.names(x) <- c(1:4, "....", ".....", (n - 3):n)
        print(x)
    }
}

is.snp.VCFinfo <- function(x)
{
    REF <- x$REF
    if (is.null(REF)) stop("no REF allele(s)")
    ALT <- x$ALT
    if (is.null(ALT)) stop("no ALT allele(s)")
    nchar(REF) == 1 & nchar(ALT) == 1
}

rangePOS <- function(x, from, to)
{
    POS <- x$POS
    if (is.null(POS)) stop("no POS(isition)")
    which(from <= POS & POS <= to)
}

selectQUAL <- function(x, threshold = 20)
{
    QUAL <- x$QUAL
    if (is.null(QUAL)) stop("no QUAL(ility)")
    which(QUAL >= threshold)
}

getINFO <- function(x, what = "DP", as.is = FALSE)
{
    INFO <- x$INFO
    if (is.null(INFO)) stop("no INFO")
    regexp <- paste0("^.*", what, "=")
    tmp <- gsub(regexp, "", INFO)
    tmp <- gsub(";.+$", "", tmp)
    if (!as.is) {
        op <- options(warn = -1)
        on.exit(options(op))
        if (!all(is.na(as.numeric(tmp[1:10]))))
            tmp <- as.numeric(tmp)
    }
    tmp
}
