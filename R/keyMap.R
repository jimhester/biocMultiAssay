setGeneric("mapIDToRanges", signature="db",
    function(db, ...) standardGeneric("mapIDToRanges")
)

#setMethod("mapIDToRanges", "Multiassay",
    #function(db, ...)
    #{
        ## Extract x, keys, keytype etc. from the Multiassay class
        ## Dispatch to method below based on 'x':
        #mapIDToRanges(db, keytype, columns, level, ...)
    #}
#)

setMethod("mapIDToRanges", "TxDb",
          function(db,
                   vals,
                   type = c("cds", "exon", "tx", "gene"),
                   columns = NULL,
                   ...)
{
    assert(is.list(vals) && is.named(vals),
        "'vals' must be a named list")

    assert(is.null(columns) || is.character(columns),
           "'column' must be 'NULL' or a character vector")

    type <- match.arg(type)

    fun <- switch(type,
                  cds = GenomicFeatures::cds,
                  exon = GenomicFeatures::exons,
                  tx = GenomicFeatures::transcripts,
                  gene = GenomicFeatures::genes)

    res <- fun(db, vals, columns = unique(c(names(vals), columns)), ...)
    matches <- BiocGenerics::match(mcols(res)[[names(vals)]], vals[[1]])
    ranges <- rep(res, vapply(matches, length, integer(1)))
    mcols(ranges)[[names(vals)]] <- NULL
    res2 <- GRangesList(split(ranges, vals[[1]][BiocGenerics::unlist(matches)]))
    res2[match(vals[[1]], names(res2))]
})

setGeneric("mapRangesToID", signature="db",
    function(db, ...) standardGeneric("mapRangesToID")
)

setMethod("mapRangesToID", "TxDb",
          function(db,
                   ranges,
                   type = c("cds", "exon", "tx", "gene"),
                   id = NULL,
                   columns = NULL,
                   ...)
{
    type <- match.arg(type)
    assert(is.null(columns) || is.character(columns),
           "'column' must be 'NULL' or a character vector")

    fun <- switch(type,
                  cds = GenomicFeatures::cds,
                  exon = GenomicFeatures::exons,
                  tx = GenomicFeatures::transcripts,
                  gene = GenomicFeatures::genes)

    all <-
        if (is.null(columns)) {
            fun(db)
        } else {
            fun(db, columns = columns)
        }

    hits <- findOverlaps(ranges, all, ...)

    res <- mcols(all[subjectHits(hits)])

    # Add the input ids and put them as the first column
    res$id <- names(ranges)[queryHits(hits)]
    res[c(NCOL(res), seq_along(res)[-NCOL(res)])]
})

#setMethod("mapIDToRanges", "OrganismDb",
    #function(db, keys, keytype = NULL, columns = NULL, level =
#c("exon", "cds", "gene", etc.), ...)
    #{
        #stopifnot(keytype)
        ## Check that keytype is in columns(x)
        ## Check that keytype is a valid map between OrgDb and TxDb objects
            #eg, ENSEMBLTRANS -> ENTREZID is not
        ## Decide what columns are needed to create GRanges based on level
        ## Dispatch to low level mapper
    #}
#)

#setMethod("mapIDToRanges", "biomaRt",
    #function(x, keys, keytype=someDefault, columns=someDefault, level =
#c("exon", "cds", "gene", etc.), ...)
    #{
     #...
    #}
#)
#' Map a set of keys and return a corresponding range
#' @param db Database to use
#' @param keytype keytype passed to \code{select()}, if \code{NULL} uses the first keytype.
#' @param keys keys passed to \code{select()}, if \code{NULL}, uses the first 3 keys.
#' @param columns type of column, one of '
#' @return \code{GRangesList} object with \code{names} corresponding to the
#' \code{keys} used.
#' @examples
#'
#' mapKey(Homo.sapiens, "ENSEMBLTRANS")
#' mapKey(Homo.sapiens)
#' # also works with orgDb objects, although the seqinfo is not populated
#' mapKey(org.Hs.eg.db)
#mapKey <- function((db), key, ...)

#TxDb <- function(db, keys = "ENL", keytype = c("cds", "exon", "tx", "gene"), by = keytype) { #}c("cds", "exon", "tx", "gene")) {
    #if (keytype == by) {
        #switch(keytype,
               #tx = transcript,
               #exon = exon,
               #)
    #} else {
    #switch(keytype,
           #tx = transcriptsBy(db, keys, by = by, ...)
           #exon = exonsBy(db, keys, by = by, ...)
         #)
    #}

mapKey <- function(db, vals, valstype = c("cds", "exon", "transcript", "gene"), type = c("cds", "exon", "transcript", "gene"), columns = NULL) {

    if (is.null(columns)) {
        type <- match.arg(type)

        prefix <- switch(type,
            cds = "CDS",
            exon = "EXON",
            transcript = "TX"
            )

        columns <- paste0(prefix, c("CHROM", "START", "END", "STRAND"))
    }

    if (is.null(keytype)) {
        keytype <- keytypes(db)[1L]
    }

    if (is.null(keys)) {
        keys <- head(keys(db, keytype), n = 3)
    }
    seqInfo <- tryCatch(seqinfo(db), error = function(e) NULL)

    res <- suppressMessages(
        AnnotationDbi::select(db,
            keytype = keytype,
            columns = columns,
            keys = keys))

    split(
        makeGRangesFromDataFrame(res[, -1, drop = FALSE],
            seqinfo = seqInfo,
            keep.extra.columns = TRUE,
            seqnames.field = columns[1],
            start.field = columns[2],
            end.field = columns[3],
            strand.field = columns[4]),
        res[[keytype]])
}

assert <- function(x, message) {
    if(!x) {
        stop(message, call. = FALSE)
    }
}

is.named <- function(x) {
    nm <- names(x)
    !is.null(nm) && all(!is.na(nm) & nm != "")
}

