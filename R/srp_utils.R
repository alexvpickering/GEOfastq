#' Get GSE text from GEO
#'
#' @param gse_name  GEO study name to get metadata for
#'
#' @return Character vector of lines on GSE record.
#' @export
#'
#' @examples
#' gse_text <- crawl_gse('GSE111459')
#'
crawl_gse <- function(gse_name) {

  # get html text for GSE page
  gse_url  <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", gse_name, "&targ=self&form=text&view=full")

  gse_text <- NULL
  attempt <- 1
  while(is.null(gse_text) && attempt <= 3) {
    con <- url(gse_url)
    try(gse_text <- readLines(con))
    if(is.null(gse_text)) Sys.sleep(15)
    close(con)
  }
  return(gse_text)
}

#' Extract GSMs needed to download RNA-seq data for a series
#'
#'
#' @param gse_text GSE text returned from \code{\link{crawl_gse}}
#'
#' @return Character vector of sample GSMs for the series \code{gse_name}
#' @export
#'
#' @examples
#'
#' gse_text <- crawl_gse('GSE111459')
#' gsm_names <- extract_gsms(gse_text)
#'
extract_gsms <- function(gse_text) {

  # GSM names
  gsm_lines <- grep('^!Series_sample_id', gse_text)
  gsm_names <- gsub('^!Series_sample_id = (GSM\\d+)$', '\\1', gse_text[gsm_lines])

  return(gsm_names)
}

#' Crawls SRX pages for each GSM to get metadata.
#'
#' Goes to each GSM page to get SRX then to each SRX page to get some more metadata.
#'
#'
#' @param samples Character vector of GSMs.
#' @param max.workers Maximum number of parallel workers to split task betweem
#' @importFrom foreach %dopar%
#'
#' @return data.frame
#' @export
#'
#' @examples
#' srp_meta <- crawl_gsms("GSM3031462")
#'
#' # returns NULL because records on dbGAP for privacy reasons
#' srp_meta <- crawl_gsms("GSM2439650")
#'
crawl_gsms <- function(gsm_names, max.workers = 50) {

  nsamp <- length(gsm_names)
  cat(nsamp, 'GSMs to process\n')

  cl <- parallel::makeCluster(min(max.workers, nsamp))
  doParallel::registerDoParallel(cl)

  srp_meta <- foreach::foreach(j=1:nsamp, .combine = plyr::rbind.fill) %dopar% {

    # save in srp_meta
    srp_meta <- data.frame(stringsAsFactors = FALSE)

    gsm_name <- gsm_names[j]
    # get html text
    gsm_url  <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", gsm_name, '&targ=self&form=text&view=full')
    gsm_text <- NULL
    attempt <- 1
    while(is.null(gsm_text) && attempt <= 3) {
      con <- url(gsm_url)
      try(gsm_text <- readLines(con))
      if(is.null(gsm_text)) Sys.sleep(5)
      close(con)
    }

    # get SRA number for this GSM
    # some won't (e.g. submitted to dbGAP for privacy reasons)
    has.srx <- grep('!Sample_relation = SRA: .+?SRX\\d+', gsm_text)
    if (!length(has.srx)) return(NULL)

    experiment <- gsub('^.+?(SRX\\d+)$', '\\1', gsm_text[has.srx])

    info <- gsub('^!Sample_', '', gsm_text[-1])
    cols <- gsub('^(.+?) = .+?$', '\\1', info)
    cols <- make.unique(cols)
    vals <- gsub('^.+? = (.+?$)', '\\1', info)
    names(vals) <- cols


    if (!is.na(experiment)) {

      # extract SRR runs info -----
      srx_url  <- paste0("https://www.ncbi.nlm.nih.gov/sra/", experiment, '[accn]?report=FullXml')

      srx_html <- NULL
      attempt <- 1
      while(is.null(srx_html) && attempt <= 3) {
        try(srx_html <- xml2::read_html(srx_url))
        if(is.null(srx_html)) Sys.sleep(5)
      }
      srx_text <- rvest::html_text(srx_html)

      runs <- stringr::str_extract_all(srx_text, '<PRIMARY_ID>SRR\\d+</PRIMARY_ID>')[[1]]
      runs <- gsub('<PRIMARY_ID>(SRR\\d+)</PRIMARY_ID>', '\\1', unique(runs))
      taxon_id <- stringr::str_extract(srx_text, '<TAXON_ID>\\d+</TAXON_ID>')
      library_source <- stringr::str_extract(srx_text, '<LIBRARY_SOURCE>.+?</LIBRARY_SOURCE>')
      library_layout <- stringr::str_extract(srx_text, '<LIBRARY_LAYOUT>.+?</LIBRARY_LAYOUT>')
      library_layout <- stringr::str_extract(library_layout, 'SINGLE|PAIRED')

      if (length(runs)) {

        # add info to srp_meta
        srp_meta[runs, 'run'] <- runs
        srp_meta[runs, 'experiment'] <- experiment
        srp_meta[runs, 'gsm_name'] <- gsm_name
        for (col in cols) srp_meta[runs, col] <- vals[col]

        srp_meta[runs, 'library_source'] <- gsub('<LIBRARY_SOURCE>(.+?)</LIBRARY_SOURCE>', '\\1', library_source)
        srp_meta[runs, 'library_layout'] <- library_layout
        srp_meta[runs, 'taxon_id'] <- gsub('<TAXON_ID>(\\d+)</TAXON_ID>', '\\1', taxon_id)
        srp_meta[runs, 'ebi_dir']  <- sapply(runs, get_dldir, 'ebi')
        srp_meta[runs, 'ncbi_dir'] <- sapply(runs, get_dldir, 'ncbi')
      }
    }
    return(srp_meta)
  }
  parallel::stopCluster(cl)
  row.names(srp_meta) <- srp_meta$run
  return(srp_meta)
}


#' Gets part of path to download bulk RNAseq sample from EBI or NCBI
#'
#' @param srr SRR/ERR run name
#' @param type Either \code{'ebi'} or \code{'ncbi'}
#'
#' @return String path used by \code{\link{get_fastqs}}.
#' @export
get_dldir <- function(srr, type = c('ebi', 'ncbi')) {


  dir1 <- substr(srr, 1, 6)

  if (type[1] == 'ebi') {
    digits  <- gsub('^SRR|^ERR', '', srr)
    ndigits <- nchar(digits)

    if (ndigits == 7) {
      dir2 <- paste0('00', substr(digits, 7, 7))

    } else if (ndigits == 8) {
      dir2 <- paste0('0', substr(digits, 7, 8))

    } else if (ndigits == 9) {
      dir2 <- substr(digits, 7, 9)

    } else if (ndigits == 6) {
      return(file.path(dir1, srr))
    }

    return(file.path(dir1, dir2, srr))

  } else if (type[1] == 'ncbi') {
    return(file.path(dir1, srr))
  }

}


#' Download and RNA-seq fastq data from EBI
#'
#' First tries to get RNA-Seq fastq files from EBI.
#'
#' @param srp_meta \code{data.frame} with SRP meta info. Returned from \code{\link{crawl_gsms}}.
#' @param data_dir Path to folder that fastq files will be downloaded to. Will be created if doesn't exist.
#' @param method One of \code{'aspera'} or \code{'ftp'}. \code{'aspera'} is generally faster but requires the
#'  ascp command line utility to be on your path and in the authors experience frequently stalls.
#' @param max_rate Used when \code{method = 'aspera'} only. Sets the target transfer rate. If downloads are stalling,
#'  reducing the default \code{'1g'} down to e.g. \code{'300m'} may help.
#'
#' @export
get_fastqs <- function(srp_meta, data_dir, method = c('ftp', 'aspera'), max_rate = '1g') {

  # setup fastq directory
  dir.create(data_dir)

  # seperate runs based on GSM (can be multiple per GSM)
  srr_names <- srp_meta$run
  gsm_names <- unique(srp_meta$gsm_name)
  srr_names_list <- lapply(gsm_names, function(gsm_name) srr_names[srp_meta$gsm_name %in% gsm_name])
  names(srr_names_list) <- gsm_names

  method <- method[1]
  ngsm <- length(gsm_names)

  res <- c()
  for (i in 1:ngsm) {
    # download everything
    srr_names <- srr_names_list[[i]]
    resi <- c()

    for (srr_name in srr_names) {
      # try to get fastq from ebi
      resi <- c(
        resi,
        tryCatch(get_ebi_fastqs(srp_meta, srr_name, data_dir, method = method, max_rate = max_rate),
                 error = function(e) return(1)))
    }
    names(resi) <- srr_names
    res <- c(res, resi)
  }
  return(res)
}

#' Download fastqs from EBI
#'
#' Much faster to use aspera than ftp
#'
#' @param srr_name Run accession as string.
#' @inheritParams get_fastqs

#' @export
#'
get_ebi_fastqs <- function(srp_meta, srr_name, data_dir, method = c('ftp', 'aspera'), max_rate = '1g') {
  url <- paste0('ftp://ftp.sra.ebi.ac.uk/vol1/fastq/', srp_meta[srr_name, 'ebi_dir'], '/')
  resp <- RCurl::getURL(url)
  resp <- strsplit(resp, '\n')[[1]]
  resp <- strsplit(resp, ' +')
  fnames <- sapply(resp, `[`, 9)
  fsizes <- sapply(resp, `[`, 5)

  if (method[1] == 'aspera') {
    ascp_path <- system('which ascp', intern = TRUE)
    ascp_pubkey <- gsub('bin/ascp$', 'etc/asperaweb_id_dsa.openssh', ascp_path)

    # only overwrite if different from source
    ascpCMD <- paste('ascp --overwrite=diff -k1 -QT -l', max_rate, '-P33001 -i', ascp_pubkey)
    files <- paste0('era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/', srp_meta[srr_name, 'ebi_dir'], '/', fnames)
    for (i in seq_along(files)) {
      destfile <- file.path(data_dir, fnames[i])
      if (file.exists(destfile) && file.size(destfile) == fsizes[i]) res <- 0
      else res <- ascpR(ascpCMD, files[i], data_dir)
    }

  } else if (method[1] == 'ftp') {
    files <- paste0('ftp://ftp.sra.ebi.ac.uk/vol1/fastq/', srp_meta[srr_name, 'ebi_dir'], '/', fnames)
    for (i in seq_along(files)) {
      # check for existing file with same size
      destfile <- file.path(data_dir, fnames[i])
      if (file.exists(destfile) && file.size(destfile) == fsizes[i]) res <- 0
      else res <- download.file(files[i], destfile)
    }
  }
  return(res)
}

#' Utility function to run aspera
#'
#' @param ascpCMD aspera command as string.
#' @param file Url to aspera file to download.
#' @param destDir Path to directory to download \code{files} into.
#'
#' @export
ascpR <- function (ascpCMD, file, destDir = getwd()) {
  ascp_cmd <- paste(ascpCMD, src, destDir, sep = " ")
  ret <- system(ascp_cmd)
  return(ret)
}
