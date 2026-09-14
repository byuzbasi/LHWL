#!/usr/bin/env Rscript
# Local maintainer utility, excluded from the CRAN source tarball.
# inst/CITATION is the bibliographic source of truth; DESCRIPTION owns version.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L || !args[1] %in% c("--check", "--write")) {
  stop("Usage: Rscript --vanilla tools/update_documentation.R --check|--write PACKAGE_DIR")
}
mode <- args[1]
root <- normalizePath(args[2], mustWork = TRUE)
if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Maintainer dependency 'yaml' is required; this script never installs it.")
}
meta <- as.list(read.dcf(file.path(root, "DESCRIPTION"))[1, ])
stopifnot(identical(meta$Package, "LHWL"), identical(meta$Encoding, "UTF-8"))
refs <- utils::readCitationFile(file.path(root, "inst", "CITATION"), meta = meta)
entries <- unclass(refs)
keys <- vapply(entries, function(x) attr(x, "key"), character(1))
types <- vapply(entries, function(x) attr(x, "bibtype"), character(1))
if (anyDuplicated(keys) || any(!grepl("^[A-Za-z][A-Za-z0-9_-]*$", keys))) {
  stop("Every citation needs a unique stable alphanumeric BibTeX key.")
}
if (sum(keys == "LHWL") != 1L || types[keys == "LHWL"] != "Manual") {
  stop("Keep the software Manual record with key LHWL.")
}
if (any(!types %in% c("Manual", "Article", "Misc", "Book"))) {
  stop("Unsupported reference type; extend the explicit CFF mapping before use.")
}
for (entry in entries) {
  if (is.null(entry$title) || is.null(entry$year) || is.null(entry$author)) {
    stop("Every record needs a title, author list and year.")
  }
  if (!is.null(entry$url) && !grepl("^https://", entry$url)) {
    stop("Reference URLs must use HTTPS.")
  }
  if (!is.null(entry$doi) && !grepl("^10\\.[0-9]{4,9}/\\S+$", entry$doi)) {
    stop("DOI must be a bare 10.xxxx/... identifier.")
  }
}
authors_cff <- function(x) {
  lapply(seq_along(x), function(i) {
    a <- x[[i]]
    out <- list("family-names" = a$family,
                "given-names" = paste(a$given, collapse = " "))
    if (is.null(a$family)) out <- list(name = paste(a$given, collapse = " "))
    out
  })
}
software <- entries[[which(keys == "LHWL")]]
cff <- list(
  "cff-version" = "1.2.0",
  message = paste(
    "Please cite LHWL and, when available, the primary LHWR methods paper.",
    "For an extension, additionally cite the corresponding methods paper."
  ),
  type = "software", title = software$title, version = meta$Version,
  authors = authors_cff(software$author),
  "repository-code" = software$url, url = software$url,
  license = "GPL-3.0-or-later",
  abstract = paste(
    "Local regression and penalized local regression over representations on",
    "the Lorentz model of two-dimensional hyperbolic space, with bandwidth and",
    "regularization selection, diagnostics, and visualization."
  ),
  keywords = c("hyperbolic geometry", "local regression", "Lorentz model",
               "penalized regression", "variable selection")
)
papers <- entries[keys != "LHWL"]
if (length(papers)) {
  cff$references <- lapply(papers, function(e) {
    type <- switch(attr(e, "bibtype"),
                   Article = "article", Book = "book",
                   Misc = "generic", Manual = "manual")
    out <- list(type = type, title = e$title, authors = authors_cff(e$author),
                year = as.integer(e$year))
    for (field in c("doi", "url", "journal", "volume", "issue")) {
      source <- if (field == "issue") "number" else field
      if (!is.null(e[[source]])) out[[field]] <- e[[source]]
    }
    if (!is.null(e$publisher)) out$publisher <- list(name = e$publisher)
    if (!is.null(e$pages)) {
      pages <- strsplit(e$pages, "--|-", perl = TRUE)[[1]]
      out$start <- pages[1]
      if (length(pages) == 2L) out$end <- pages[2]
    }
    if (!is.null(e$note)) out$notes <- e$note
    out
  })
}

readme_path <- file.path(root, "README.md")
readme <- readLines(readme_path, warn = FALSE, encoding = "UTF-8")
start <- which(readme == "<!-- LHWL-REFERENCES:START -->")
end <- which(readme == "<!-- LHWL-REFERENCES:END -->")
if (length(start) != 1L || length(end) != 1L || start >= end) {
  stop("README must contain exactly one ordered citation marker pair.")
}
formatted <- vapply(seq_along(refs), function(i) {
  paste(format(refs[[i]], style = "text"), collapse = " ")
}, character(1))
reference_block <- c(
  "<!-- Generated from inst/CITATION by tools/update_documentation.R. -->",
  "", "### Citation records", "",
  paste0("- ", formatted), ""
)
generated_readme <- c(readme[seq_len(start)], reference_block,
                      readme[seq.int(end, length(readme))])
outputs <- list(
  "vignettes/references.bib" = paste(
    c("% Generated from inst/CITATION; do not edit independently.",
      as.character(utils::toBibtex(refs))), collapse = "\n"
  ),
  "CITATION.cff" = sub("\n$", "", yaml::as.yaml(cff, indent.mapping.sequence = TRUE)),
  "README.md" = paste(generated_readme, collapse = "\n")
)
changed <- character()
for (name in names(outputs)) {
  target <- file.path(root, name)
  wanted <- strsplit(enc2utf8(outputs[[name]]), "\n", fixed = TRUE)[[1]]
  current <- if (file.exists(target)) {
    readLines(target, warn = FALSE, encoding = "UTF-8")
  } else character()
  if (!identical(current, wanted)) {
    changed <- c(changed, name)
    if (mode == "--write") writeLines(wanted, target, useBytes = TRUE)
  }
}
if (mode == "--check" && length(changed)) {
  stop("Generated documentation is stale: ", paste(changed, collapse = ", "),
       ". Review the source records, then run --write.")
}
cat(if (mode == "--write") "Updated" else "Verified",
    length(outputs), "generated documents from", length(refs), "citation record(s).\n")
if (length(changed)) cat("Changed:", paste(changed, collapse = ", "), "\n")
