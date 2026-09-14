#!/usr/bin/env Rscript
# Build locally using an explicitly selected, already installed candidate.
# No install, deployment, remote task, or production analysis is performed.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript --vanilla tools/build_documentation.R PACKAGE_DIR NEW_OUTPUT_DIR LIBRARY_DIR")
}
root <- normalizePath(args[1], mustWork = TRUE)
out <- path.expand(args[2])
lib <- normalizePath(args[3], mustWork = TRUE)
if (file.exists(out)) stop("Output exists; choose a new versioned directory: ", out)
if (!dir.exists(dirname(out))) stop("Output parent directory does not exist.")
out <- file.path(normalizePath(dirname(out)), basename(out))
if (startsWith(out, paste0(root, .Platform$file.sep))) {
  stop("Write the site outside the frozen source candidate.")
}
# Confine build caches to this validation run, never the user's global cache.
Sys.setenv(R_USER_CACHE_DIR = file.path(dirname(out), "documentation_cache"))
required <- c("pkgdown", "rmarkdown", "knitr", "yaml")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing maintainer dependencies: ", paste(missing, collapse = ", "))
if (!rmarkdown::pandoc_available()) stop("Pandoc is unavailable.")
.libPaths(c(lib, .libPaths()))
installed <- find.package("LHWL", lib.loc = lib, quiet = TRUE)
if (!length(installed)) stop("Install the validated candidate into LIBRARY_DIR first.")
meta <- read.dcf(file.path(root, "DESCRIPTION"))
installed_meta <- read.dcf(file.path(installed, "DESCRIPTION"))
if (meta[1, "Version"] != installed_meta[1, "Version"]) {
  stop("Installed LHWL version differs from source.")
}
receipt <- file.path(lib, "LHWL_SOURCE_SIGNATURE.rds")
# Paths are listed independently of timestamps; compare exact executable source.
source_paths <- sort(c("DESCRIPTION", "NAMESPACE",
  paste0("R/", list.files(file.path(root, "R"), pattern = "[.]R$")),
  paste0("src/", list.files(file.path(root, "src"),
                           pattern = "[.](cpp|h)$|^Makevars([.]win)?$"))
))
signature <- unname(tools::md5sum(file.path(root, source_paths)))
names(signature) <- source_paths
if (!file.exists(receipt) || !identical(readRDS(receipt), signature)) {
  stop("Missing or stale installation receipt LHWL_SOURCE_SIGNATURE.rds; ",
       "create it only after a successful install of these exact sources.")
}
status <- system2(file.path(R.home("bin"), "Rscript"), c(
  "--vanilla", shQuote(file.path(root, "tools", "update_documentation.R")),
  "--check", shQuote(root)
))
if (status != 0L) stop("Citation consistency check failed.")
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           VECLIB_MAXIMUM_THREADS = "1", R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
pkg <- pkgdown::as_pkgdown(root, override = list(destination = out))
pkgdown::build_site(
  pkg, examples = FALSE, seed = 20260910L, lazy = FALSE,
  preview = FALSE, devel = FALSE, new_process = FALSE, install = FALSE
)
cat("Local documentation written to:", out, "\n")
