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
mapKey <- function(db, keytype = NULL, keys = NULL, type = c("cds", "exon", "transcript"), columns = NULL) {

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
