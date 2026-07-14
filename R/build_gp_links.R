utils::globalVariables(c("Gene"))

#' Build Gene-Peak Network (Enhancer, Promoter, Closest)
#'
#' @param pk Either:
#'   (1) A character vector of genomic coordinates ("chr1:1000-1200"), or
#'   (2) A data.frame containing peak coordinates with columns
#'       "chr","start","end" or a "PeakRegion" column
#'       ("chr1:1000-1200" style).
#' @param gn Character vector of gene symbols
#'   (HGNC official symbols, e.g. "TP53", "FMNL2"). Each input peak is linked
#'   to the nearest transcript TSS among these genes.
#' @param ver Integer (17, 18, or 19) - PANTHER enhancer-gene link version.
#' @param enh_file Optional local path to enhancer-gene link file (.tsv).
#' @param enh_coords_file Optional local path to the corresponding PEREGRINE
#'   enhancer-coordinate file. The file must contain chromosome, start, end,
#'   and enhancer identifier columns without a header. If `NULL`, the hg38
#'   PEREGRINE coordinate file is downloaded or retrieved from cache.
#' @param gene_map Optional data.frame with columns `hgnc_symbol` and `hgnc_id`.
#'   When supplied, this local mapping is used instead of querying BioMart.
#' @param mart Optional `biomaRt::Mart` object used for gene identifier mapping.
#'   If `NULL`, the default Ensembl human genes mart is used. Supply an object
#'   created by [biomaRt::useEnsembl()] to select a mirror or custom host.
#'
#' @return A data.frame with columns: Peak, Gene, Src. Closest-gene mappings
#'   assign each peak center to its single nearest transcript TSS across the
#'   genes in `gn`. When a peak-gene pair has multiple evidence sources, each
#'   source is returned as a separate row.
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @export
#'
#'
#' @examples
#' data("gp_example_inputs", package = "GPlinksR")
#' pk <- gp_example_inputs$pk[seq_len(3)]
#' gn <- gp_example_inputs$gn[seq_len(4)]
#'
#' enh_file <- system.file(
#'     "extdata", "gp_example_enhancer_links.tsv",
#'     package = "GPlinksR", mustWork = TRUE
#' )
#' enh_coords_file <- system.file(
#'     "extdata", "gp_example_enhancer_coordinates.tsv",
#'     package = "GPlinksR", mustWork = TRUE
#' )
#' gene_map <- data.frame(
#'     hgnc_symbol = c("TTLL10", "SAMD11"),
#'     hgnc_id = c("HGNC:26693", "HGNC:28706")
#' )
#'
#' gp <- build_gp_links(
#'     pk, gn,
#'     enh_file = enh_file,
#'     enh_coords_file = enh_coords_file,
#'     gene_map = gene_map
#' )
#' gp
#'
#' if (interactive()) {
#'     mirror_mart <- biomaRt::useEnsembl(
#'         biomart = "genes",
#'         dataset = "hsapiens_gene_ensembl",
#'         mirror = "useast"
#'     )
#'     gp <- build_gp_links(pk, gn, mart = mirror_mart)
#'     head(gp)
#' }
build_gp_links <- function(
    pk,
    gn,
    ver = 19,
    enh_file = NULL,
    enh_coords_file = NULL,
    gene_map = NULL,
    mart = NULL
) {
    message("Checking inputs...")
    .validate_build_gp_inputs(
        pk = pk,
        gn = gn,
        ver = ver,
        enh_file = enh_file,
        enh_coords_file = enh_coords_file,
        gene_map = gene_map
    )
    msg <- paste0(
        "Input validation passed: ", length(pk),
        " peaks and ", length(gn), " genes provided."
    )
    message(msg)
    message("Preparing peaks...")
    pk_gr <- .prepare_peak_granges(pk)
    message("Building enhancer-based links...")
    if (is.null(enh_file)) {
        enh_file <- get_peregrine_file(ver)
    }
    df_enh <- .build_enhancer_links(
        pk_gr = pk_gr,
        gn = gn,
        enh_file = enh_file,
        enh_coords_file = enh_coords_file,
        gene_map = gene_map,
        mart = mart
    )
    .report_link_count("Enhancer-gene links found: ", df_enh)
    message("Building promoter-based links...")
    edb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
    df_prom <- .build_promoter_links(pk_gr = pk_gr, gn = gn, edb = edb)
    .report_link_count("Promoter-gene links found: ", df_prom)
    message("Building closest-gene links...")
    df_clo <- .build_closest_links(pk_gr = pk_gr, gn = gn, edb = edb)
    .report_link_count("Closest-gene links found: ", df_clo)
    df_all <- .combine_link_tables(df_enh, df_prom, df_clo)
    .report_link_count("Total combined links: ", df_all)
    df_all
}


.report_link_count <- function(label, links) {
    message(label, nrow(links))
}


.validate_build_gp_inputs <- function(
    pk,
    gn,
    ver,
    enh_file,
    enh_coords_file,
    gene_map
) {
    .validate_gene_input(gn)
    .validate_peak_input(pk)
    if (!ver %in% c(17, 18, 19)) {
        stop("'ver' must be one of 17, 18, or 19 (PEREGRINE version numbers).")
    }

    .validate_optional_file(enh_file, "enh_file")
    .validate_optional_file(enh_coords_file, "enh_coords_file")

    if (!is.null(gene_map)) {
        .prepare_gene_map(gene_map, gn)
    }
}


.validate_optional_file <- function(path, argument) {
    if (is.null(path) || file.exists(path)) {
        return(invisible(NULL))
    }

    msg <- paste0("The provided '", argument, "' path does not exist: ", path)
    stop(msg)
}


.validate_gene_input <- function(gn) {
    if (missing(gn) || is.null(gn) || length(gn) == 0) {
        msg <- paste0(
            "'gn' must be provided - a non-empty character vector ",
            "of gene symbols."
        )
        stop(msg)
    }

    if (!is.character(gn)) {
        stop("'gn' must be a character vector (e.g., c('TP53', 'FMNL2')).")
    }

    if (any(is.na(gn) | gn == "")) {
        msg <- paste0(
            "'gn' contains missing or empty entries. ",
            "Please remove or replace them."
        )
        stop(msg)
    }
}


.validate_peak_input <- function(pk) {
    if (missing(pk) || is.null(pk) || length(pk) == 0) {
        msg <- paste0(
            "'pk' must be provided - either a character vector ",
            "or a data.frame."
        )
        stop(msg)
    }

    if (is.character(pk)) {
        .validate_peak_strings(pk)
        return(invisible(NULL))
    }

    if (is.data.frame(pk)) {
        .validate_peak_dataframe(pk)
        return(invisible(NULL))
    }

    pk_class <- paste(class(pk), collapse = ", ")
    msg <- paste0(
        "'pk' must be a character vector or a data.frame, not ", pk_class
    )
    stop(msg)
}


.validate_peak_strings <- function(pk) {
    message_template <- paste0(
        "Invalid peak format: '%s'. Expected 'chr#:start-end' ",
        "(e.g., 'chr1:1000-1200')."
    )
    .validate_peak_format(pk, message_template)
}


.validate_peak_format <- function(pk, message_template) {
    is_valid <- grepl("^chr[0-9XYMT]+:[0-9]+-[0-9]+$", pk)
    if (all(is_valid)) {
        return(invisible(NULL))
    }

    stop(sprintf(message_template, pk[!is_valid][1]))
}


.validate_peak_dataframe <- function(pk) {
    cols <- colnames(pk)
    has_required <- "PeakRegion" %in% cols ||
        all(c("chr", "start", "end") %in% cols)

    if (has_required) {
        return(invisible(NULL))
    }

    msg <- paste0(
        "Peak data.frame must contain either 'PeakRegion' or columns ",
        "'chr','start','end'."
    )
    stop(msg)
}


.prepare_peak_granges <- function(pk) {
    if (is.data.frame(pk)) {
        pk_dt <- .peak_dataframe_to_dt(pk)
    } else if (is.character(pk)) {
        pk_dt <- .peak_vector_to_dt(pk)
    } else {
        msg <- paste0(
            "Unsupported peak input type. Provide a character vector ",
            "or data.frame."
        )
        stop(msg)
    }

    pk_dt$start <- as.integer(pk_dt$start)
    pk_dt$end <- as.integer(pk_dt$end)
    pk_dt$id <- paste0(pk_dt$chr, ":", pk_dt$start, "-", pk_dt$end)
    pk_dt$chr <- ifelse(
        !grepl("^chr", pk_dt$chr),
        paste0("chr", pk_dt$chr),
        pk_dt$chr
    )

    pk_gr <- makeGRangesFromDataFrame(
        pk_dt,
        seqnames.field = "chr",
        start.field = "start",
        end.field = "end",
        keep.extra.columns = TRUE
    )

    valid_chr <- GenomeInfoDb::standardChromosomes(
        pk_gr,
        species = "Homo sapiens"
    )
    GenomeInfoDb::keepSeqlevels(
        pk_gr,
        valid_chr,
        pruning.mode = "coarse"
    )
}


.peak_dataframe_to_dt <- function(pk) {
    if ("PeakRegion" %in% colnames(pk)) {
        pk_split <- data.table::tstrsplit(pk$PeakRegion, "[:-]")
        return(data.table::data.table(
            chr = pk_split[[1]],
            start = pk_split[[2]],
            end = pk_split[[3]]
        ))
    }

    if (all(c("chr", "start", "end") %in% colnames(pk))) {
        return(data.table::data.table(
            chr = pk$chr,
            start = pk$start,
            end = pk$end
        ))
    }

    msg <- paste0(
        "Peak input must contain either 'PeakRegion' or columns ",
        "'chr','start','end'."
    )
    stop(msg)
}


.peak_vector_to_dt <- function(pk) {
    pk_split <- data.table::tstrsplit(pk, "[:-]")
    data.table::data.table(
        chr = pk_split[[1]],
        start = pk_split[[2]],
        end = pk_split[[3]]
    )
}


.build_enhancer_links <- function(
    pk_gr,
    gn,
    enh_file,
    enh_coords_file = NULL,
    gene_map = NULL,
    mart = NULL
) {
    enh_raw <- .read_enhancer_link_table(enh_file)

    gn_map <- .map_genes_to_hgnc(
        gn,
        gene_map = gene_map,
        mart = mart
    )
    gn_ids <- sub("HGNC:", "", gn_map$hgnc_id)
    enh_sub <- enh_raw[enh_raw$HGNC %in% gn_ids, ]

    enh_reg <- .read_enhancer_coordinates(enh_coords_file)
    enh_full <- merge(enh_sub, enh_reg, by = "enhancer", all.x = TRUE)
    enh_full$sym <- gn_map$hgnc_symbol[match(enh_full$HGNC, gn_ids)]

    enh_gr <- makeGRangesFromDataFrame(
        enh_full,
        seqnames.field = "chr",
        start.field = "start",
        end.field = "end",
        keep.extra.columns = TRUE
    )

    enh_gr <- .align_seqlevels(enh_gr, pk_gr)

    .overlap_df(
        query_gr = pk_gr,
        subject_gr = enh_gr,
        subject_gene = S4Vectors::mcols(enh_gr)$sym,
        src = "enh"
    )
}


.read_enhancer_link_table <- function(enh_file) {
    enh_raw <- data.table::fread(
        enh_file,
        sep = "\t",
        header = TRUE,
        fill = TRUE,
        quote = ""
    )
    enh_raw$HGNC <- sub(".*HGNC=([0-9]+).*", "\\1", enh_raw$gene)
    enh_raw$enhancer <- as.character(enh_raw$enhancer)
    enh_raw
}


.read_enhancer_coordinates <- function(enh_coords_file) {
    message("Loading enhancer coordinates...")
    if (is.null(enh_coords_file)) {
        enh_coords_file <- .get_peregrine_enhancer_coords_file()
    }
    data.table::fread(
        enh_coords_file,
        sep = "\t",
        header = FALSE,
        col.names = c("chr", "start", "end", "enhancer")
    )
}


.map_genes_to_hgnc <- function(gn, gene_map = NULL, mart = NULL) {
    if (!is.null(gene_map)) {
        return(.prepare_gene_map(gene_map, gn))
    }

    if (is.null(mart)) {
        mart <- biomaRt::useMart(
            "ensembl",
            dataset = "hsapiens_gene_ensembl"
        )
    }
    gn_map <- biomaRt::getBM(
        attributes = c("hgnc_symbol", "hgnc_id"),
        filters = "hgnc_symbol",
        values = gn,
        mart = mart
    )

    if (nrow(gn_map) == 0) {
        stop("No HGNC IDs found for the provided gene symbols.")
    }

    gn_map
}


.prepare_gene_map <- function(gene_map, gn) {
    if (!is.data.frame(gene_map)) {
        stop("'gene_map' must be a data.frame.")
    }

    required <- c("hgnc_symbol", "hgnc_id")
    missing_columns <- setdiff(required, colnames(gene_map))
    if (length(missing_columns) > 0) {
        msg <- paste0(
            "'gene_map' must contain columns: ",
            paste(required, collapse = ", "), "."
        )
        stop(msg)
    }

    gene_map <- data.frame(
        hgnc_symbol = as.character(gene_map$hgnc_symbol),
        hgnc_id = as.character(gene_map$hgnc_id),
        stringsAsFactors = FALSE
    )
    complete <- !is.na(gene_map$hgnc_symbol) &
        nzchar(gene_map$hgnc_symbol) &
        !is.na(gene_map$hgnc_id) &
        nzchar(gene_map$hgnc_id)
    gene_map <- unique(gene_map[complete, , drop = FALSE])
    gene_map <- gene_map[gene_map$hgnc_symbol %in% gn, , drop = FALSE]

    if (nrow(gene_map) == 0) {
        stop("No HGNC IDs found for the provided gene symbols.")
    }

    gene_map
}


.build_promoter_links <- function(pk_gr, gn, edb) {
    tx_gr <- ensembldb::transcripts(
        edb,
        columns = c("tx_id", "gene_id", "gene_name", "seq_name")
    )

    prom_gr <- .promoter_ranges(tx_gr)
    prom_gr$symbol <- tx_gr$gene_name
    prom_gr <- prom_gr[prom_gr$symbol %in% gn]
    prom_gr <- .align_seqlevels(prom_gr, pk_gr)

    .overlap_df(
        query_gr = pk_gr,
        subject_gr = prom_gr,
        subject_gene = prom_gr$symbol,
        src = "prom",
        ignore_strand = TRUE
    )
}


.promoter_ranges <- function(tx_gr) {
    prom_gr <- withCallingHandlers(
        GenomicRanges::promoters(tx_gr, upstream = 2000, downstream = 200),
        warning = function(w) {
            if (grepl("out-of-bound", conditionMessage(w))) {
                invokeRestart("muffleWarning")
            }
        }
    )
    GenomicRanges::trim(prom_gr)
}


.build_closest_links <- function(pk_gr, gn, edb) {
    tx_gr <- ensembldb::transcripts(
        edb,
        columns = c("tx_id", "gene_id", "gene_name", "seq_name")
    )
    tss_gr <- .transcript_tss(tx_gr)
    .closest_links_from_tss(pk_gr, tss_gr, gn)
}


.transcript_tss <- function(tx_gr) {
    tss_gr <- GenomicRanges::promoters(tx_gr, upstream = 0, downstream = 1)
    tss_gr$symbol <- tx_gr$gene_name

    valid <- !is.na(tss_gr$symbol) & nzchar(tss_gr$symbol)
    tss_gr <- tss_gr[valid]
    key <- paste(
        GenomeInfoDb::seqnames(tss_gr),
        IRanges::start(tss_gr),
        IRanges::end(tss_gr),
        GenomicRanges::strand(tss_gr),
        tss_gr$symbol,
        sep = ":"
    )
    tss_gr[!duplicated(key)]
}


.closest_links_from_tss <- function(pk_gr, tss_gr, gn) {
    tss_gr <- tss_gr[tss_gr$symbol %in% gn]
    tss_gr <- .align_seqlevels(tss_gr, pk_gr)
    pk_gr4 <- .subset_common_seqlevels(pk_gr, tss_gr)

    pk_gr4 <- GenomicRanges::trim(pk_gr4)
    tss_gr <- GenomicRanges::trim(tss_gr)

    if (length(tss_gr) == 0 || length(pk_gr4) == 0) {
        return(.empty_link_df())
    }

    pk_center_gr <- .peak_center_ranges(pk_gr4)
    nearest_hits <- GenomicRanges::nearest(
        pk_center_gr,
        tss_gr,
        ignore.strand = TRUE
    )

    data.frame(
        Peak = pk_center_gr$id,
        Gene = tss_gr$symbol[nearest_hits],
        Src = "clo",
        stringsAsFactors = FALSE
    ) |>
        dplyr::filter(!is.na(Gene)) |>
        dplyr::distinct()
}


.peak_center_ranges <- function(pk_gr) {
    center <- IRanges::start(pk_gr) + floor(IRanges::width(pk_gr) / 2)
    GenomicRanges::GRanges(
        seqnames = GenomeInfoDb::seqnames(pk_gr),
        ranges = IRanges::IRanges(start = center, width = 1),
        strand = "*",
        id = pk_gr$id
    )
}


.align_seqlevels <- function(subject_gr, query_gr) {
    subject_gr <- .set_seqlevels_style(subject_gr, query_gr)

    common <- intersect(
        GenomeInfoDb::seqlevels(query_gr),
        GenomeInfoDb::seqlevels(subject_gr)
    )
    GenomeInfoDb::keepSeqlevels(subject_gr, common, pruning.mode = "coarse")
}


.set_seqlevels_style <- function(subject_gr, query_gr) {
    subject_levels <- GenomeInfoDb::seqlevels(subject_gr)
    query_levels <- GenomeInfoDb::seqlevels(query_gr)
    subject_has_chr <- any(startsWith(subject_levels, "chr"))
    query_has_chr <- any(startsWith(query_levels, "chr"))

    if (query_has_chr && !subject_has_chr) {
        renamed_levels <- paste0("chr", subject_levels)
        renamed_levels[subject_levels == "MT"] <- "chrM"
        GenomeInfoDb::seqlevels(subject_gr) <- renamed_levels
    } else if (!query_has_chr && subject_has_chr) {
        renamed_levels <- sub("^chr", "", subject_levels)
        renamed_levels[renamed_levels == "M"] <- "MT"
        GenomeInfoDb::seqlevels(subject_gr) <- renamed_levels
    }

    subject_gr
}


.subset_common_seqlevels <- function(query_gr, subject_gr) {
    common <- intersect(
        GenomeInfoDb::seqlevels(query_gr),
        GenomeInfoDb::seqlevels(subject_gr)
    )
    GenomeInfoDb::keepSeqlevels(query_gr, common, pruning.mode = "coarse")
}


.overlap_df <- function(query_gr,
                        subject_gr,
                        subject_gene,
                        src,
                        ignore_strand = FALSE) {
    query_gr2 <- .subset_common_seqlevels(query_gr, subject_gr)
    hits <- GenomicRanges::findOverlaps(
        query_gr2,
        subject_gr,
        ignore.strand = ignore_strand
    )

    if (length(hits) == 0) {
        return(.empty_link_df())
    }

    data.frame(
        Peak = query_gr2$id[S4Vectors::queryHits(hits)],
        Gene = subject_gene[S4Vectors::subjectHits(hits)],
        Src = src,
        stringsAsFactors = FALSE
    ) |>
        dplyr::distinct()
}


.empty_link_df <- function() {
    data.frame(
        Peak = character(0),
        Gene = character(0),
        Src = character(0),
        stringsAsFactors = FALSE
    )
}


.combine_link_tables <- function(df_enh, df_prom, df_clo) {
    dplyr::bind_rows(df_enh, df_prom, df_clo) |>
        dplyr::distinct()
}
