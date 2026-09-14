# Maintaining LHWL documentation and citations

This is a local, manual release workflow. It does not publish to GitHub,
deploy a website, submit to CRAN, install dependencies, or start a cloud task.
Keep historical manuscript replication bundles and validated release
snapshots unchanged. Start each revision in a new versioned source directory.

The canonical package repository is https://github.com/byuzbasi/LHWL. Use
that address for package URLs, issue reports, software citations and GitHub
installation examples. LHWR remains the method name. The repository is public;
CRAN submission remains a separate, explicitly authorized step.
Keep private anonymous-review OSF links out of
the public package, website and arXiv source files.

## Documentation sources

- `vignettes/lorentz-hgwr-introduction.Rmd`: the main English guide.
- `vignettes/lhwr-diagnostics-and-tuning.Rmd`: LHWR kernel and bandwidth
  contracts, selection criteria, F1--F3, local inference and Moran/LISA.
- `vignettes/lhwr-graphics.Rmd`: manuscript-style LHWR figure recipes and
  non-overwriting exports; helpers in this guide are not new package APIs.
- `vignettes/lhwpr-lasso.Rmd`: the separate Lasso extension tutorial.
- `inst/CITATION`: the authoritative bibliographic records.
- `DESCRIPTION`: the authoritative package version and package metadata.
- `_pkgdown.yml` and `pkgdown/extra.css`: local navigation and style.
- `tools/update_documentation.R`: checks or regenerates the shared bibliography,
  GitHub citation metadata and the marked README citation section.
- `tools/build_documentation.R`: builds local HTML using an explicitly selected
  installation with a matching source receipt; refuses existing output paths.

The normal R package build creates the installed HTML vignettes. pkgdown uses
the same R Markdown sources for its articles. Do not hand-edit generated HTML,
`vignettes/references.bib`, `CITATION.cff` or the marked README citation block.
The site and maintainer tools are excluded from the CRAN source archive;
the Rmd files, bibliography and standard vignette outputs are included.
Documentation tools are maintainer dependencies, not runtime Imports.

## Adding a paper

1. Verify the actual public record: title, full author list, year, publication
   status, arXiv identifier/version and DOI when available. A manuscript being
   prepared or reviewed is not a published article. Do not invent identifiers.
2. Add a `bibentry()` in `inst/CITATION`, with a stable unique key. Keep the
   software `Manual` record with key `LHWL`. Supported record types are
   `Manual`, `Article`, `Misc` (including preprints), and `Book`. Every record
   requires `title`, `author`, `year` and `key`; use HTTPS URLs and bare DOIs.
   Public manuscript records should include a resolvable URL or DOI.
3. Treat LHWR as the primary methodological foundation. Cite the Lasso
   paper additionally for Lasso functionality. Add later penalty or application
   papers only to their relevant sections; do not imply that one paper proves
   all implemented extensions. A new paper does not automatically make all
   its authors package authors.
4. Add the corresponding `@key` reference and a concise explanation in the
   relevant Rmd section. Replace any now-obsolete pending-publication wording.
   Once the primary LHWR record is public, describe it explicitly as the
   primary methods reference in the main guide and citation instructions.
5. Run the citation generator with `--write`, inspect its diff, then run
   `--check`. It renders all registered citations into the bibliography,
   README citation block and CFF references without changing package authors.
   A preprint-to-journal update should normally retain its citation key.
6. Update NEWS. For an already released CRAN package, increment DESCRIPTION's
   version before preparing a new submission. Unpublished candidate revisions
   may retain their intended version, distinguished by directory and hashes.
7. Rebuild, inspect and check the exact candidate. GitHub documentation can be
   updated separately from a CRAN release; an edit on GitHub does not update
   CRAN's installed vignettes or citation record. Batch routine documentation
   changes into appropriate CRAN releases rather than making automatic
   submissions for every small edit.

The generator verifies structure and synchronization, not the truth of
publication metadata or the availability of a URL. That verification remains
an explicit editorial gate.

## Copyable local workflow for this candidate

These commands use a reviewed local source checkout and the already installed
maintainer dependencies (`pkgdown`, `knitr`, `rmarkdown`, `yaml`) plus Pandoc.
Do not install or update dependencies implicitly.

```sh
cd /absolute/path/to/LHWL_source
Rscript --vanilla tools/update_documentation.R --write .
Rscript --vanilla tools/update_documentation.R --check .
```

Run `--write` only in a new editable candidate after reviewing the citation
source. To build the site using an existing verified installation, choose a
new output directory each time:

```sh
cd /absolute/path/to/LHWL_source
Rscript --vanilla tools/build_documentation.R \
  . \
  /absolute/path/to/LHWL_site_new \
  /absolute/path/to/LHWL_doc_library_new
```

The site-build command requires a verified installation of the same candidate.
It expects
`LHWL_doc_library_new/LHWL_SOURCE_SIGNATURE.rds`, written by local validation
only after successful installation of that candidate. It contains named
source-content checksums for DESCRIPTION, NAMESPACE, R files and native source
and Makevars files. Version equality alone is insufficient because two
unpublished candidates may share a version. The release validation also
records SHA-256 manifests of source and outputs. Reinstall and create a new
receipt after changing executable source or DESCRIPTION; do not bypass a
receipt mismatch. On a different machine use its own working R toolchain.

For a new reviewed candidate, installation and receipt creation must refer to
that same source directory. The following R code records the receipt only
after the explicit local source installation succeeds. Use a new library
directory, with a working compiler configuration appropriate to your host:

```r
pkg_src <- normalizePath(".")
doc_lib <- file.path(dirname(pkg_src), "LHWL_doc_library_new")
stopifnot(!dir.exists(doc_lib))
dir.create(doc_lib)
status <- system2(file.path(R.home("bin"), "R"), c(
  "CMD", "INSTALL", "--preclean",
  paste0("--library=", shQuote(doc_lib)), shQuote(pkg_src)
))
stopifnot(status == 0L)
paths <- sort(c("DESCRIPTION", "NAMESPACE",
  paste0("R/", list.files(file.path(pkg_src, "R"), pattern = "[.]R$")),
  paste0("src/", list.files(file.path(pkg_src, "src"),
    pattern = "[.](cpp|h)$|^Makevars([.]win)?$"))
))
receipt <- unname(tools::md5sum(file.path(pkg_src, paths)))
names(receipt) <- paths
saveRDS(receipt, file.path(doc_lib, "LHWL_SOURCE_SIGNATURE.rds"))
```

Run these R commands from the intended new source candidate, never from a
frozen historical bundle. Source installation can leave compiled build
products in that editable copy; exclude them from distributed source and
preserve the final cleaned-source manifest.

The build does not install or deploy the package and disables reference-page
example evaluation; the guide examples execute normally. Standard R package
checks test the reference examples separately. HTML math may use MathJax's
standard rendering assets; document examples themselves do not fetch data.
Local source operation is not a promise of fully offline browser math assets.
The pkgdown builder may also request package metadata from CRAN or
Bioconductor and standard static web dependencies. Build caches are scoped
to `documentation_cache` next to the output directory. A restricted-network
failure is not an example failure; preserve its log and retry the documented
build in a new output directory with authorized network access.

## Before publication

- Confirm package/version/citation metadata agree across all outputs.
- Check every executable guide chunk; label non-executed installation and
  alternative-protocol snippets honestly. Do not hide errors or warnings.
- Review table widths, figure captions, equations, links and narrow-screen
  layout. Use `knitr::kable()` for vignette tables.
- Run a source build and `R CMD check --as-cran`; record host-specific notes,
  environment workarounds and unperformed checks.
- Verify hashes and preserve the final source archive and local site.
- Confirm the public repository and documentation URL before advertising them.
- Obtain separate authorization for GitHub publication, hosting or CRAN
  submission. Multi-platform and internet-enabled incoming checks remain
  release gates; a successful local documentation build does not replace them.

There is no recurring monitor or automatic publication job. New papers enter
this workflow when their verified metadata are provided.
