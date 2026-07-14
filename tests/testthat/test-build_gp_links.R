# tests/testthat/test-build_gp_links.R

.local_peregrine_files <- function() {
    enh_file <- tempfile(fileext = ".tsv")
    coords_file <- tempfile()

    write.table(
        data.frame(
            enhancer = "EH1",
            gene = "gene;HGNC=11998",
            stringsAsFactors = FALSE
        ),
        file = enh_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )

    write.table(
        data.frame(
            chr = "chr1",
            start = 1050L,
            end = 1150L,
            enhancer = "EH1"
        ),
        file = coords_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE
    )

    list(enh_file = enh_file, coords_file = coords_file)
}


# 1. Valid deterministic example (success case)

test_that(
    "build_gp_links returns valid outputs for local example peaks and genes",
    {
        local_files <- .local_peregrine_files()
        test_mart <- structure(list(host = "mirror"), class = "test_mart")
        testthat::local_mocked_bindings(
            .map_genes_to_hgnc = function(gn, gene_map, mart) {
                expect_null(gene_map)
                expect_identical(mart, test_mart)
                data.frame(
                    hgnc_symbol = "TP53",
                    hgnc_id = "HGNC:11998",
                    stringsAsFactors = FALSE
                )
            },
            .build_promoter_links = function(pk_gr, gn, edb) {
                data.frame(
                    Peak = character(0),
                    Gene = character(0),
                    Src = character(0),
                    stringsAsFactors = FALSE
                )
            },
            .build_closest_links = function(pk_gr, gn, edb) {
                data.frame(
                    Peak = character(0),
                    Gene = character(0),
                    Src = character(0),
                    stringsAsFactors = FALSE
                )
            },
            .package = "GPlinksR"
        )

        test_peaks <- c(
            "chr1:1000-1100",
            "chr1:5000-5100"
        )

        test_genes <- "TP53"

        result <- build_gp_links(
            pk = test_peaks,
            gn = test_genes,
            enh_file = local_files$enh_file,
            enh_coords_file = local_files$coords_file,
            mart = test_mart
        )

        expect_s3_class(result, "data.frame")
        expect_true(all(c("Peak", "Gene", "Src") %in% colnames(result)))
        expect_gt(nrow(result), 0)
        expect_false(any(is.na(result$Gene)))
        expect_true(all(unique(result$Src) %in% c("enh", "prom", "clo")))
        expect_true(all(result$Peak %in% test_peaks))
        expect_true(any(result$Src == "enh" & result$Gene == "TP53"))
    }
)


test_that("build_gp_links accepts local gene and enhancer resources", {
    local_files <- .local_peregrine_files()
    gene_map <- data.frame(
        hgnc_symbol = "TP53",
        hgnc_id = "HGNC:11998"
    )
    testthat::local_mocked_bindings(
        .build_promoter_links = function(pk_gr, gn, edb) {
            data.frame(
                Peak = character(),
                Gene = character(),
                Src = character()
            )
        },
        .build_closest_links = function(pk_gr, gn, edb) {
            data.frame(
                Peak = character(),
                Gene = character(),
                Src = character()
            )
        },
        .package = "GPlinksR"
    )

    result <- build_gp_links(
        pk = "chr1:1000-1100",
        gn = "TP53",
        enh_file = local_files$enh_file,
        enh_coords_file = local_files$coords_file,
        gene_map = gene_map
    )

    expect_equal(result$Peak, "chr1:1000-1100")
    expect_equal(result$Gene, "TP53")
    expect_equal(result$Src, "enh")
})


# 2. Missing gene input (gn)

test_that("build_gp_links errors when gn is missing", {
    expect_error(
        build_gp_links(pk = "chr1:1000-2000"),
        "'gn' must be provided"
    )
})


# 3. gn is not a character vector

test_that("build_gp_links errors when gn is not a character vector", {
    expect_error(
        build_gp_links(
            pk = "chr1:1000-2000",
            gn = 123
        ),
        "'gn' must be a character vector"
    )
})


# 4. gn contains NA or empty strings

test_that("build_gp_links errors when gn contains NA or empty strings", {
    expect_error(
        build_gp_links(
            pk = "chr1:1000-2000",
            gn = c("TP53", NA)
        ),
        "contains missing or empty entries"
    )

    expect_error(
        build_gp_links(
            pk = "chr1:1000-2000",
            gn = c("TP53", "")
        ),
        "contains missing or empty entries"
    )
})


# 5. Missing peak input (pk)

test_that("build_gp_links errors when pk is missing", {
    expect_error(
        build_gp_links(gn = "TP53"),
        "'pk' must be provided"
    )
})


# 6. pk is wrong type

test_that("build_gp_links errors when pk is not character or data.frame", {
    expect_error(
        build_gp_links(
            pk = 123,
            gn = "TP53"
        ),
        "'pk' must be a character vector or a data.frame"
    )
})


# 7. Invalid peak coordinate format

test_that("build_gp_links errors for malformed peak strings", {
    expect_error(
        build_gp_links(
            pk = "1:1000-2000",
            gn = "TP53"
        ),
        "Invalid peak format"
    )
})


# 8. pk data.frame missing required columns

test_that("build_gp_links errors when pk data.frame lacks required columns", {
    bad_pk <- data.frame(a = 1, b = 2)

    expect_error(
        build_gp_links(
            pk = bad_pk,
            gn = "TP53"
        ),
        "Peak data.frame must contain"
    )
})


# 9. Invalid PEREGRINE enhancer version

test_that("build_gp_links errors for invalid PEREGRINE version", {
    expect_error(
        build_gp_links(
            pk = "chr1:1000-2000",
            gn = "TP53",
            ver = 16
        ),
        "'ver' must be one of 17, 18, or 19"
    )
})


# 10. Nonexistent enhancer file

test_that("build_gp_links errors when enh_file does not exist", {
    expect_error(
        build_gp_links(
            pk = "chr1:1000-2000",
            gn = "TP53",
            enh_file = "not_a_real_file.tsv"
        ),
        "does not exist"
    )
})


test_that("build_gp_links errors when enhancer coordinates do not exist", {
    expect_error(
        build_gp_links(
            pk = "chr1:1000-2000",
            gn = "TP53",
            enh_coords_file = "not_a_real_file.tsv"
        ),
        "enh_coords_file.*does not exist"
    )
})


test_that("build_gp_links validates a supplied gene map", {
    expect_error(
        build_gp_links(
            pk = "chr1:1000-2000",
            gn = "TP53",
            gene_map = data.frame(symbol = "TP53", id = "HGNC:11998")
        ),
        "must contain columns"
    )
})


# 11. Gene symbols cannot be mapped to HGNC

test_that("build_gp_links errors when gene symbols cannot be mapped to HGNC", {
    local_files <- .local_peregrine_files()

    expect_error(
        build_gp_links(
            pk = "chr1:1000-2000",
            gn = "THISISNOTAREALGENE",
            enh_file = local_files$enh_file,
            enh_coords_file = local_files$coords_file,
            gene_map = data.frame(
                hgnc_symbol = "TP53",
                hgnc_id = "HGNC:11998"
            )
        ),
        "No HGNC IDs found"
    )
})


test_that("peak preparation keeps standard human chromosomes", {
    prepare_peak_granges <- getFromNamespace(
        ".prepare_peak_granges",
        "GPlinksR"
    )
    peaks <- data.frame(
        chr = c("chr1", "chr23", "chrUn_KI270442v1"),
        start = c(100L, 200L, 300L),
        end = c(150L, 250L, 350L)
    )

    result <- prepare_peak_granges(peaks)

    expect_equal(as.character(GenomicRanges::seqnames(result)), "chr1")
})


test_that("seqlevel styles are aligned without external lookups", {
    set_seqlevels_style <- getFromNamespace(
        ".set_seqlevels_style",
        "GPlinksR"
    )
    subject <- GenomicRanges::GRanges(
        seqnames = c("1", "MT"),
        ranges = IRanges::IRanges(start = c(1L, 1L), width = 10L)
    )
    query <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chrM"),
        ranges = IRanges::IRanges(start = c(1L, 1L), width = 10L)
    )

    result <- set_seqlevels_style(subject, query)

    expect_equal(GenomeInfoDb::seqlevels(result), c("chr1", "chrM"))
})


test_that("closest links assign every peak among the input genes", {
    closest_links_from_tss <- getFromNamespace(
        ".closest_links_from_tss",
        "GPlinksR"
    )
    peaks <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1"),
        ranges = IRanges::IRanges(start = c(110L, 10010L), width = 10L),
        id = c("near_unselected", "near_selected")
    )
    tss <- GenomicRanges::GRanges(
        seqnames = c("1", "1"),
        ranges = IRanges::IRanges(start = c(100L, 10000L), width = 1L),
        strand = "+",
        symbol = c("UNSELECTED", "SELECTED")
    )

    result <- closest_links_from_tss(peaks, tss, gn = "SELECTED")

    expect_equal(
        result,
        data.frame(
            Peak = c("near_unselected", "near_selected"),
            Gene = c("SELECTED", "SELECTED"),
            Src = c("clo", "clo"),
            stringsAsFactors = FALSE
        )
    )
    expect_equal(nrow(result), length(peaks))
})


test_that("peak centers follow the paper midpoint definition", {
    peak_center_ranges <- getFromNamespace(
        ".peak_center_ranges",
        "GPlinksR"
    )
    peaks <- GenomicRanges::GRanges(
        seqnames = c("chr1", "chr1"),
        ranges = IRanges::IRanges(
            start = c(100L, 100L),
            end = c(199L, 200L)
        ),
        id = c("even_width", "odd_width")
    )

    result <- peak_center_ranges(peaks)

    expect_equal(IRanges::start(result), c(150L, 150L))
    expect_equal(IRanges::width(result), c(1L, 1L))
    expect_equal(result$id, peaks$id)
})


test_that("combined links retain every evidence source", {
    combine_link_tables <- getFromNamespace(
        ".combine_link_tables",
        "GPlinksR"
    )
    link <- data.frame(Peak = "chr1:100-200", Gene = "GENE1")
    enhancer <- transform(link, Src = "enh")
    promoter <- transform(link, Src = "prom")
    closest <- transform(link, Src = "clo")

    result <- combine_link_tables(enhancer, promoter, closest)

    expect_equal(nrow(result), 3L)
    expect_equal(sort(result$Src), c("clo", "enh", "prom"))
})


test_that("packaged example links match the executable offline example", {
    data("gp_example_inputs", package = "GPlinksR")
    data("gp_example_links", package = "GPlinksR")
    enh_file <- system.file(
        "extdata", "gp_example_enhancer_links.tsv",
        package = "GPlinksR", mustWork = TRUE
    )
    coords_file <- system.file(
        "extdata", "gp_example_enhancer_coordinates.tsv",
        package = "GPlinksR", mustWork = TRUE
    )
    gene_map <- data.frame(
        hgnc_symbol = c("TTLL10", "SAMD11"),
        hgnc_id = c("HGNC:26693", "HGNC:28706")
    )

    result <- build_gp_links(
        pk = gp_example_inputs$pk[seq_len(3)],
        gn = gp_example_inputs$gn[seq_len(4)],
        enh_file = enh_file,
        enh_coords_file = coords_file,
        gene_map = gene_map
    )

    expect_equal(result, gp_example_links)
})
