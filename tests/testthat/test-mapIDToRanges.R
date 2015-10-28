context("mapIDToRanges")
txdb <- makeTxDbFromGRanges(readRDS("gr.rds"))

test_that("Improper inputs throw errors", {
    expect_error(mapIDToRanges(txdb, vals = ""),
        "must be a named list")

    expect_error(mapIDToRanges(txdb, vals = "ENST000000271582"),
        "must be a named list")

    expect_error(mapIDToRanges(txdb, vals = list("ENST000000271582")),
        "must be a named list")

    expect_error(mapIDToRanges(txdb, vals = list(tx_name = "ENST000000271582"), column = 1),
        "'column' must be 'NULL' or a character vector")
})

test_that("it returns results in same order as input", {
    vals = list(tx_name = c("ENST00000371582", "ENST00000371588", "ENST00000494752", "ENST00000614008", "ENST00000496771"))
    res <- mapIDToRanges(txdb, vals = vals, type = "tx")
    expect_equal(names(res), vals[[1]])

    # shuffle the order and make sure it remains equivalent
    for (i in seq_len(10)) {
        vals <- sample(vals)
        res <- mapIDToRanges(txdb, vals = vals, type = "tx")
        expect_equal(names(res), vals[[1]])
    }
})

test_that("it returns NA if no results are found", {
    vals <- list(tx_name = c("ENST00000371582", "NOT_FOUND", "ENST00000494752"))
    res <- mapIDToRanges(txdb, vals = vals, type = "tx")
    expect_equal(names(res)[2], NA_character_)

    # shuffle the order and make sure it remains equivalent
    for (i in seq_len(10)) {
        vals <- sample(vals)
        res <- mapIDToRanges(txdb, vals = vals, type = "tx")
        expect_equal(names(res)[grep("NOT_FOUND", vals[[1]])], NA_character_)
    }
})

test_that("it returns duplicate ranges if needed", {
    # both of these transcripts are from the same gene
    vals <- list(tx_name = c("ENST00000371582", "ENST00000494752"))
    res <- mapIDToRanges(txdb, vals = vals, type = "gene")

    #names match input
    expect_equal(names(res), vals[[1]])

    # but values are the same
    expect_equivalent(res[1], res[2])
})
