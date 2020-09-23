#' Fix overlapping characteristics columns
#'
#' @param pdata Data.frame returned by \code{\link{crawl_gsms}}
#'
#' @return \code{pdata} with characteristics columns aligned
#' @export
#'
#' @examples
#'
#' gse_name <- 'GSE69529'
#' gse_text <- crawl_gse(gse_name)
#' gsm_names <- extract_gsms(gse_text)
#'
#' # note characteristics columns misaligned
#' srp_meta <- crawl_gsms(gsm_names)
#'
#' # this fixes them
#' srp_meta <- fix_characteristics(srp_meta)
#
fix_characteristics <- function(pdata) {

  # characteristics column names
  char_cols  <- grep('characteristics', colnames(pdata), value = TRUE)
  if (!length(char_cols)) return(pdata)

  # ensure space after colons on characteristics columns
  pdata[, char_cols] <- lapply(pdata[, char_cols, drop=FALSE], function(col) gsub(':([^ ])', ': \\1', col))

  char_names <- lapply(pdata[, char_cols], function(col) stringr::str_extract_all(col, '^[^:]+: ', simplify = TRUE))
  char_names <- do.call(cbind, char_names)
  char_names <- unique(char_names[char_names != ""])
  char_names <- char_names[!is.na(char_names)]

  # merge characteristics columns with characteristics vals
  char_vals <- apply(pdata[, char_cols, drop=FALSE], 1, function(row) paste(row[grepl('^[^:]+: ', row)], collapse = " "))

  if (!length(char_names)) {
    pdata[, char_cols] <- NULL
    return(pdata)
  }

  # char_names in each char_val
  in_vals  <- lapply(char_vals, function(char_val) sapply(char_names, function(char_name) grepl(paste0('\\b', esc(char_name)), char_val)))
  in_names <- lapply(in_vals, function(x) names(x)[x])


  # split based on characteristic names
  char_vals <- mapply(function(char_val, in_name) {

    split1 <- sapply(in_name, function(char_name) strsplit(char_val, esc(char_name))[[1]][2])
    split2 <- sapply(split1, function(x) strsplit(x, paste(esc(char_names), collapse='|'))[[1]][1])
    gsub('^[ ,]+|[ ,]+$', '', split2)

  }, char_vals, in_names, SIMPLIFY = FALSE)

  # merge with pdata
  new_cols <- paste0('characteristics_', seq_along(char_names))
  pdata[, new_cols] <- NA

  for (i in seq_along(in_vals)) {
    pdata[i, new_cols[in_vals[[i]]]] <- paste0(char_names[in_vals[[i]]], char_vals[[i]])
  }

  pdata[, char_cols] <- NULL

  return(pdata)
}

#' Preserve escape characters in grep patterns
#'
#' @param string Character vector with pattern
#'
#' @return \code{string} with escape characters preserved
#' @export
#'
esc <- function(string) {
  gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", string)
}
