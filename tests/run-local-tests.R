pkgload::load_all(".", quiet = TRUE)
testthat::test_dir("tests/testthat")
