# Helpers for creating SRA submissions

instruments <- list(
    ILLUMINA = c("HiSeq X Five",
                 "HiSeq X Ten",
                 "Illumina Genome Analyzer",
                 "Illumina Genome Analyzer II",
                 "Illumina Genome Analyzer IIx",
                 "Illumina HiScanSQ",
                 "Illumina HiSeq 1000",
                 "Illumina HiSeq 1500",
                 "Illumina HiSeq 2000",
                 "Illumina HiSeq 2500",
                 "Illumina HiSeq 3000",
                 "Illumina HiSeq 4000",
                 "Illumina iSeq 100",
                 "Illumina NovaSeq 6000",
                 "Illumina MiniSeq",
                 "Illumina MiSeq",
                 "NextSeq 500",
                 "NextSeq 550"),
    HELICOS = c("Helicos HeliScope"),
    ABI_SOLID = c("AB 5500 Genetic Analyzer",
                  "AB 5500xl Genetic Analyzer",
                  "AB 5500x-Wl Genetic Analyzer",
                  "AB SOLiD 3 Plus System",
                  "AB SOLiD 4 System",
                  "AB SOLiD 4hq System",
                  "AB SOLiD PI System",
                  "AB SOLiD System",
                  "AB SOLiD System 2.0",
                  "AB SOLiD System 3.0"),
    COMPLETE_GENOMICS = c("Complete Genomics"),
    PACBIO_SMRT = c("PacBio RS",
                    "PacBio RS II",
                    "PacBio Sequel"),
    ION_TORRENT = c("Ion Torrent PGM",
                    "Ion Torrent Proton",
                    "Ion Torrent S5 XL",
                    "Ion Torrent S5"),
    CAPILLARY = c("AB 310 Genetic Analyzer",
                  "AB 3130 Genetic Analyzer",
                  "AB 3130xL Genetic Analyzer",
                  "AB 3500 Genetic Analyzer",
                  "AB 3500xL Genetic Analyzer",
                  "AB 3730 Genetic Analyzer",
                  "AB 3730xL Genetic Analyzer"),
    OXFORD_NANOPORE = c("GridION", "MinION", "PromethION"),
    BGISEQ = c("BGISEQ-500")
)

gut_envo <- c(
    env_broad_scale = paste0("environment associated with an animal part ",
                             "or small animal [ENVO:01001055]"),
    env_local_scale = "intestine environment [ENVO:2100002]",
    env_medium = "fecal material [ENVO:00002003]"
)

presets = list(
    `human gut 16S` = c(
        gut_envo,
        organism = "human gut metagenome",
        host = "Homo sapiens",
        library_strategy = "AMPLICON",
        library_source = "GENOMIC",
        library_selection = "PCR",
        filetype = "fastq",
        usage = paste0("You are now ready for submission. ",
                       "Go to https://submit.ncbi.nlm.nih.gov/subs/sra/, ",
                       "log in and click on `New submission`. ",
                       "Fill in the general data for your project in steps ",
                       "1 through 3. In step 4 ",
                       "choose `Genome, metagenome or marker sequences ",
                       "(MIxS compliant)` and `Survey-related Marker ",
                       "Sequences MIMARKS`. In step 5 and 6 ",
                       "upload the respective files in %s. In step 7 you ",
                       "can directly upload the `*.tar.gz` submission ",
                       "package. Just click on `continue` another time to ",
                       "have the archive unpacked as indicated.")
    ),
    `mouse gut 16S` = c(
        gut_envo,
        organism = "mouse gut metagenome",
        host = "Mus musculus",
        library_strategy = "AMPLICON",
        library_source = "GENOMIC",
        library_selection = "PCR",
        filetype = "fastq",
        usage = paste0("You are now ready for submission. ",
                       "Go to https://submit.ncbi.nlm.nih.gov/subs/sra/, ",
                       "log in and click on `New submission`. ",
                       "Fill in the general data for your project in steps ",
                       "1 through 3. In step 4 ",
                       "choose `Genome, metagenome or marker sequences ",
                       "(MIxS compliant)` and `Survey-related Marker ",
                       "Sequences MIMARKS`. In step 5 and 6 ",
                       "upload the respective files in %s. In step 7 you ",
                       "can directly upload the `*.tar.gz` submission ",
                       "package. Just click on `continue` another time to ",
                       "have the archive unpacked as indicated.")
    ),
    `human gut metagenome` = c(
        gut_envo,
        organism = "human gut metagenome",
        host = "Homo sapiens",
        library_strategy = "WGS",
        library_source = "METAGENOMIC",
        library_selection = "RANDOM",
        filetype = "fastq",
        usage = paste0("You are now ready for submission. ",
                       "**Please make sure you only upload sequence files ",
                       "where human sequences have been removed and (you can ",
                       "use the `filter_reference` workflow for that)!**",
                       "Go to https://submit.ncbi.nlm.nih.gov/subs/sra/, ",
                       "log in and click on `New submission`. ",
                       "Fill in the general data for your project in steps ",
                       "1 through 3. In step 4 ",
                       "choose `Genome, metagenome or marker sequences ",
                       "(MIxS compliant)` and `Environmental/Metagenome ",
                       "Genomic Sequences MIMS`. In step 5 and 6 ",
                       "upload the respective files in %s. In step 7 you ",
                       "can directly upload the `*.tar.gz` submission ",
                       "package. Just click on `continue` another time to ",
                       "have the archive unpacked as indicated.")
    ),
    `mouse gut metagenome` = c(
        gut_envo,
        organism = "mouse gut metagenome",
        host = "Mus musculus",
        library_strategy = "WGS",
        library_source = "METAGENOMIC",
        library_selection = "RANDOM",
        filetype = "fastq",
        usage = paste0("You are now ready for submission. ",
                       "Go to https://submit.ncbi.nlm.nih.gov/subs/sra/, ",
                       "log in and click on `New submission`. ",
                       "Fill in the general data for your project in steps ",
                       "1 through 3. In step 4 ",
                       "choose `Genome, metagenome or marker sequences ",
                       "(MIxS compliant)` and `Environmental/Metagenome ",
                       "Genomic Sequences MIMS`. In step 5 and 6 ",
                       "upload the respective files in %s. In step 7 you ",
                       "can directly upload the `*.tar.gz` submission ",
                       "package. Just click on `continue` another time to ",
                       "have the archive unpacked as indicated.")
    ),
    `human gut RNA-Seq` = c(
        gut_envo,
        organism = "human gut metagenome",
        host = "Homo sapiens",
        library_strategy = "RNA-Seq",
        library_source = "METATRANSCRIPTOMIC",
        library_selection = "RANDOM",
        filetype = "fastq",
        usage = paste0("You are now ready for submission. ",
                       "**Please make sure you only upload sequence files ",
                       "where human sequences have been removed and (you can ",
                       "use the `filter_reference` workflow for that)!**",
                       "Go to https://submit.ncbi.nlm.nih.gov/subs/sra/, ",
                       "log in and click on `New submission`. ",
                       "Fill in the general data for your project in steps ",
                       "1 through 3. In step 4 ",
                       "choose `Genome, metagenome or marker sequences ",
                       "(MIxS compliant)` and `Environmental/Metagenome ",
                       "Genomic Sequences MIMS`. In step 5 and 6 ",
                       "upload the respective files in %s. In step 7 you ",
                       "can directly upload the `*.tar.gz` submission ",
                       "package. Just click on `continue` another time to ",
                       "have the archive unpacked as indicated.")
    ),
    `mouse gut RNA-Seq` = c(
        gut_envo,
        organism = "mouse gut metagenome",
        host = "Mus musculus",
        library_strategy = "RNA-Seq",
        library_source = "METATRANSCRIPTOMIC",
        library_selection = "RANDOM",
        filetype = "fastq",
        usage = paste0("You are now ready for submission. ",
                       "Go to https://submit.ncbi.nlm.nih.gov/subs/sra/, ",
                       "log in and click on `New submission`. ",
                       "Fill in the general data for your project in steps ",
                       "1 through 3. In step 4 ",
                       "choose `Genome, metagenome or marker sequences ",
                       "(MIxS compliant)` and `Environmental/Metagenome ",
                       "Genomic Sequences MIMS`. In step 5 and 6 ",
                       "upload the respective files in `%s`. In step 7 you ",
                       "can directly upload the `*.tar.gz` submission ",
                       "package. Just click on `continue` another time to ",
                       "have the archive unpacked as indicated.")
    )
)

validate <- function(config) {
    if (is.null(config$metadata)) {
        stop("Need to specify a metadata files in config :(")
    }
    if (is.null(config$title)) {
        stop(paste0("Please specify a title in the form ",
                    "`{methodology} of {organism}: {sample info}`"))
    }
    if (!config$preset %in% names(presets)) {
        stop(sprintf("Not a supported presets. Allowed are: %s",
                     paste0(names(presets), collapse = ",")))
    }
    if (!toupper(config$platform) %in% names(instruments)) {
        stop(sprintf("Not a recognized platform, please pick one of: %s",
             paste0(names(instruments), collapse = ",")))
    }
    if (!config$instrument_model %in%
        instruments[[toupper(config$platform)]]) {
        stop(sprintf("Not a recognized model, please pick one of: %s",
             paste0(instruments[[toupper(config$platform)]], collapse = ",")))
    }
}

#' Build a configuration for the SRA submission workflow.
#'
#' This can be saved and passed on to others to ensure reproducibility.
#'
#' @param ... Any arguments are used to update the default configuration. See
#'  the example below. Optional.
#' @return A list with the parameters used in the transcript counting
#'  workflow.
#' @export
#' @examples
#'  config <- config_sra(date_col = "collection_date")
config_sra <- function(...) {
    config <- list(
        metadata = NULL,
        id_col = "id",
        date_col = "date",
        country = "USA",
        preset = "human gut 16S",
        out_dir = "sra",
        title = NULL,
        platform = "ILLUMINA",
        instrument_model = "Illumina MiSeq",
        bioproject = NULL,
        country = "USA",
        latitude = 42.36,
        longitude = -71.0941,
        make_package = TRUE
    )
    args <- list(...)
    for (arg in names(args)) {
        config[[arg]] <- args[[arg]]
    }
    return(config)
}


#' Prepare submission files for the NCBI sequence read archive (SRA).
#'
#' @param object An experiment data table as returned by any alignment method
#'  like \code{\link{align_short_reads}} or \code{\link{align_long_reads}} .
#' @param config A configuration as generated by \code{\link{config_count}}.
#' @return A list containing the used alignments and the transcript counts in
#'  `counts`.
#'
#' @export
#' @importFrom utils tar
#' @importFrom stringr str_replace_all
sra_submission <- function(object, config) {
    files <- get_files(object)
    validate(config)
    if (!dir.exists(config$out_dir)) {
        dir.create(config$out_dir, recursive = TRUE)
    }
    meta <- as.data.table(config$metadata)
    preset <- presets[[config$preset]]
    files <- copy(files)
    setkey(files, id)
    if (!all(meta[[config$id_col]] %in% files$id)) {
        stop("Some ids in the metadata have no matching files :(")
    }
    if (!all(files$id %in% meta[[config$id_col]])) {
        file.info(paste0("Some file ids are not mentioned in the metadata. ",
                         "Will assume that is intended and omit those."))
    }
    files <- files[config$metadata[[config$id_col]]]
    date <- as.Date(meta[[config$date_col]])
    date <- format(date, "%Y-%m-%d")
    sample_data <- data.table(
        sample_name = files$id,
        organism = preset["organism"],
        collection_date = date,
        env_broad_scale = preset["env_broad_scale"],
        env_local_scale = preset["env_local_scale"],
        env_medium = preset["env_medium"],
        host = preset["host"],
        geo_loc_name = config$country,
        lat_lon = paste(round(abs(config$latitude), 4),
                        ifelse(config$latitude >= 0, "N", "S"),
                        round(abs(config$longitude), 4),
                        ifelse(config$longitude >= 0, "E", "W")),
        source_material_id = files$id
    )
    additional <- names(meta)[!names(meta) %in%
                              c(config$id_col, config$date_col)]
    additional <- meta[, additional, with = FALSE]
    names(additional) <- tolower(str_replace_all(names(additional),
                                 "[^A-Za-z0-9_]", "_"))
    sample_data <- cbind(sample_data, additional)

    sra_metadata <- data.table(
        sample_name = files$id,
        library_ID = files$id,
        filename = basename(files$forward),
        title = config$title,
        library_strategy = preset["library_strategy"],
        library_source = preset["library_source"],
        library_selection = preset["library_selection"],
        filetype = preset["filetype"],
        library_layout = ifelse("reverse" %in% names(files),
                                "paired", "single"),
        platform = toupper(config$platform),
        instrument_model = config$instrument_model,
        design_description = "see manuscript for methods"
    )
    if (!is.null(config$bioproject)) {
        sample_data[, "bioproject_accession" := config$bioproject]
        sra_metadata[, "bioproject_accession" := config$bioproject]
    }
    upload = files$forward
    if ("reverse" %in% names(files)) {
        sra_metadata[, "filename2" := basename(files$reverse)]
        upload <- c(upload, files$reverse)
    }

    if (config$make_package) {
        if (Sys.getenv("tar") == "") {
            tar <- Sys.getenv("TAR")
        } else {
            tar <- Sys.getenv("tar")
        }
        flog.info("Packing submission files to %s.",
                file.path(config$out_dir, "sra_files.tar.gz"))
        ret <- tar(file.path(config$out_dir, "sra_files.tar.gz"),
            files = upload,
            compression = "gzip", tar = tar)
        if (ret != 0) {
            stop("tar command failed")
        }
    }
    flog.info("Writing biosample attributes to %s.",
              file.path(config$out_dir, "05_biosample_attributes.tsv"))
    fwrite(sample_data, sep="\t", file.path(config$out_dir,
                                            "05_biosample_attributes.tsv"))
    flog.info("Writing file metadata to %s.",
              file.path(config$out_dir, "06_sra_metadata.tsv"))
    fwrite(sra_metadata, sep="\t",
           file.path(config$out_dir, "06_sra_metadata.tsv"))
    flog.info(sprintf(preset["usage"], config$out_dir))
    artifact <- list(
        files = files,
        biosample_attributes = sample_data,
        sra_metadata = sra_metadata,
        upload = ifelse(config$make_package,
                        file.path(config$out_dir, "sra_files.tar.gz"),
                        upload)
    )
}
