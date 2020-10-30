test_that("crawl_gse can get GSE metadata", {
  expect_error(crawl_gse('GSE111459'), NA)
})

test_that("extract_gsms parses GSMs", {
    expect_equal(extract_gsms("!Series_sample_id = GSM3926903"),
                 'GSM3926903')
})

test_that("crawl_gsms gets GSM metadata", {
    srp_meta <- crawl_gsms('GSM3926903')
    expect_equal(nrow(srp_meta), 1)

    expect_equal(srp_meta$ebi_dir, 'SRR964/004/SRR9640294')
})

test_that("get_fastqs can download a fastq files", {
    skip_on_bioc()

    # fake srp_meta
    run = 'SRR014242'
    srp_meta <- data.frame(
        run  = run,
        gsm_name = run,
        ebi_dir = get_dldir(run),
        row.names = run,
        stringsAsFactors = FALSE)

    data_dir <- tempdir()
    res <- get_fastqs(srp_meta, data_dir)
    expect_equal(res, c(SRR014242 = 0))
})
