test_that("extracted values are cleaned consistently", {
    clean_values <- getFromNamespace(".clean_extracted_values", "GPlinksR")

    result <- clean_values(c("TP53", NA, "", "TP53"), "failed")

    expect_identical(result, "TP53")
    expect_error(clean_values(c(NA, ""), "failed"), "failed", fixed = TRUE)
})


test_that("rowData columns are extracted consistently", {
    row_data_info <- getFromNamespace(".row_data_info", "GPlinksR")
    extract_column <- getFromNamespace(".extract_row_data_column", "GPlinksR")
    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(seq_len(4), nrow = 2)),
        rowData = S4Vectors::DataFrame(symbol = c("TP53", "GATA1"))
    )

    info <- row_data_info(se)

    expect_identical(extract_column(info, "symbol", "gene_col"),
                     c("TP53", "GATA1"))
    expect_error(
        extract_column(info, "missing", "gene_col"),
        "'gene_col' was not found in rowData",
        fixed = TRUE
    )
})


test_that("peak format validation preserves caller-specific errors", {
    validate_peak_format <- getFromNamespace(
        ".validate_peak_format",
        "GPlinksR"
    )

    expect_no_error(validate_peak_format("chr1:100-200", "Bad peak: %s"))
    expect_error(
        validate_peak_format("not-a-peak", "Bad peak: %s"),
        "Bad peak: not-a-peak",
        fixed = TRUE
    )
})


test_that("PEREGRINE resources reuse the shared cache helper", {
    get_cached_file <- getFromNamespace(
        ".get_cached_peregrine_file",
        "GPlinksR"
    )
    cache_dir <- tempfile("gplinksr-cache-")
    dir.create(cache_dir)
    source_file <- tempfile("gplinksr-source-")
    writeLines("test resource", source_file)
    bfc <- BiocFileCache::BiocFileCache(cache = cache_dir, ask = FALSE)

    first <- NULL
    expect_message(
        first <- get_cached_file(
            source_file,
            "test-resource",
            "Downloading test resource",
            "Using cached test resource",
            bfc = bfc
        ),
        "Downloading test resource"
    )
    second <- NULL
    expect_message(
        second <- get_cached_file(
            source_file,
            "test-resource",
            "Downloading test resource",
            "Using cached test resource",
            bfc = bfc
        ),
        "Using cached test resource"
    )

    expect_identical(first, second)
    expect_true(file.exists(first))
})
