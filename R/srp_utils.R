
#' Get GSMs needed to download RNA-seq data for a series
#'
#' Goes to a series GSE page to get sample GSMs.
#'
#' @param gse_name GEO study name to get metadata for
#'
#' @return Character vector of sample GSMs for the series \code{gse_name}
#' @export
#'
#' @examples
#'
#' gsm_names <- get_gsms('GSE111459')
#'
get_gsms <- function(gse_name, data_dir = NULL) {

  if (!is.null(data_dir)) {
    gse_dir <- file.path(data_dir, gse_name)
    if (!dir.exists(gse_dir)) dir.create(gse_dir)

    # load srp meta from file if exists
    srp_meta_path <- file.path(gse_dir, 'srp_meta.rds')
    if (file.exists(srp_meta_path)) return(readRDS(srp_meta_path))
  }

  # get GSM names ----

  # get html text for GSE page
  gse_url  <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", gse_name)

  gse_html <- NULL
  attempt <- 1
  while(is.null(gse_html) && attempt <= 3) {
    try(gse_html <- xml2::read_html(gse_url))
    if(is.null(gse_html)) Sys.sleep(15)
  }
  gse_text <- rvest::html_text(gse_html)

  # GSM names
  samples <- stringr::str_extract(gse_text, stringr::regex('\nSamples.+?\nRelations', dotall = TRUE))
  samples <- stringr::str_extract_all(samples, 'GSM\\d+')[[1]]

  return(samples)
}

#' Crawls SRX pages for each GSM to get metadata.
#'
#' Goes to each GSM page to get SRX then to each SRX page to get some more metadata.
#'
#'
#' @param samples Character vector of GSMs.
#' @importFrom foreach %dopar%
#'
#' @return data.frame
#' @export
#'
#' @examples
#' gsm_names <- get_gsms('GSE111459')
#' srp_meta <- crawl_gsms(gsm_names)
#'
crawl_gsms <- function(gsm_names) {

  nsamp <- length(gsm_names)
  cat(nsamp, 'GSMs to process\n')

  cl <- parallel::makeCluster(nsamp)
  doParallel::registerDoParallel(cl)

  srp_meta <- foreach::foreach(j=1:nsamp, .combine = rbind) %dopar% {

    # save in srp_meta
    srp_meta <- data.frame(stringsAsFactors = FALSE)

    gsm_name <- gsm_names[j]
    # get html text
    gsm_url  <- paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", gsm_name, '&targ=self&form=text&view=full')
    gsm_html <- NULL
    attempt <- 1
    while(is.null(gsm_html) && attempt <= 3) {
      try(gsm_html <- xml2::read_html(gsm_url))
      if(is.null(gsm_html)) Sys.sleep(5)
    }
    gsm_text <- rvest::html_text(gsm_html)

    # get SRA number for this GSM
    experiment <- stringr::str_extract(gsm_text, 'SRX\\d+\r')
    experiment <- gsub('(SRX\\d+)\r', '\\1', experiment)

    info <- strsplit(gsm_text, '!Sample_')[[1]]
    info <- gsub('\r\n', '', info)

    cols <- gsub('^(.+?) = .+?$', '\\1', info)[-1]
    cols <- make.unique(cols)
    vals <- gsub('^.+? = (.+?$)', '\\1', info)[-1]
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
  return(srp_meta)
}


#' Gets part of path to download bulk RNAseq sample from EBI or NCBI
#'
#' @param srr SRR run name
#' @param type Either \code{'ebi'} or \code{'ncbi'}
#'
#' @return String path used by \code{\link{get_fastqs}}.
#' @export
get_dldir <- function(srr, type = c('ebi', 'ncbi')) {


  dir1 <- substr(srr, 1, 6)

  if (type[1] == 'ebi') {
    digits  <- gsub('^SRR', '', srr)
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


#' Download and fasterq-dump RNA-seq data from GEO
#'
#' First tries to get RNA-Seq fastq files from EBI. If not successful, gets SRA from GEO and converts to fastq.gz.
#'
#' @param gse_name GSE name. Will create folder with this name in \code{data_dir} and download data there.
#' @param srp_meta \code{data.frame} with SRP meta info. Returned from \code{\link{get_srp_meta}}.
#' @param data_dir Path to folder that \code{gse_name} folder will be created in.
#' @param method One of \code{'aspera'} or \code{'ftp'}. \code{'aspera'} is generally faster but requires the
#'  ascp command line utility to be on your path. \code{'ftp'} is recommended if multiple fastqs are required as
#'  a parallel for loop can give speeds comparable to \code{'aspera'}.
#'
#' @export
get_fastqs <- function(gse_name, srp_meta, data_dir = getwd(), method = c('aspera', 'ftp')) {

  # setup gse directory
  gse_dir  <- file.path(data_dir, gse_name)
  dir.create(gse_dir)


  # seperate runs based on GSM (can be multiple per GSM)
  srr_names <- srp_meta$run
  gsm_names <- unique(srp_meta$gsm_name)
  srr_names_list <- lapply(gsm_names, function(gsm_name) srr_names[srp_meta$gsm_name %in% gsm_name])
  names(srr_names_list) <- gsm_names

  # download everything
  for (i in seq_along(srr_names_list)) {
    srr_names <- srr_names_list[[i]]

    for (srr_name in srr_names) {
      # try to get fastq from ebi
      get_ebi_fastqs(srp_meta, srr_name, gse_dir, method = method[1])
    }
  }
  return(NULL)
}

#' Download fastqs from EBI
#'
#' Much faster to use aspera than ftp
#'
#' @param srp_meta Result from \code{get_srp_meta}.
#' @param srr_name Run accession as string.
#' @param gse_dir Folder to save fastq.gz files in
#' @param method One of either \code{'aspera'} (Default) or \code{'ftp'}.
#'
#' @export
#'
get_ebi_fastqs <- function(srp_meta, srr_name, gse_dir, method = c('aspera', 'ftp')) {
  url <- paste0('ftp://ftp.sra.ebi.ac.uk/vol1/fastq/', srp_meta[srr_name, 'ebi_dir'], '/.')
  fnames <- unlist(strsplit(RCurl::getURL(url, dirlistonly = TRUE), '\n'))

  if (method[1] == 'aspera') {
    ascp_path <- system('which ascp', intern = TRUE)
    ascp_pubkey <- gsub('bin/ascp$', 'etc/asperaweb_id_dsa.openssh', ascp_path)

    ascpCMD <- paste('ascp --overwrite=diff -k1 -QT -l 1g -P33001 -i', ascp_pubkey)
    files <- paste0('era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/', srp_meta[srr_name, 'ebi_dir'], '/', fnames)
    ascpR(ascpCMD, files, gse_dir)

  } else if (method[1] == 'ftp') {
    files <- paste0('ftp://ftp.sra.ebi.ac.uk/vol1/fastq/', srp_meta[srr_name, 'ebi_dir'], '/', fnames)
    for (i in seq_along(files)) {
      download.file(files[i], file.path(gse_dir, fnames[i]))
    }
  }
}

#' Utility function to run aspera
#'
#' @param ascpCMD aspera command as string.
#' @param files Urls to aspera files to download.
#' @param destDir Path to directory to download \code{files} into.
#'
#' @export
ascpR <- function (ascpCMD, files, destDir = getwd()) {
  ret <- c()
  for (src in files) {
    ascp_cmd <- paste(ascpCMD, src, destDir, sep = " ")
    ret <- c(ret, system(ascp_cmd))
  }
  return(ret)
}
