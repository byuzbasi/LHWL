## September 14, 2026: LHWL 0.3.1 GitHub candidate

The canonical public repository is https://github.com/byuzbasi/LHWL.
This revision includes the approved Gaussian/exponential Moran and LISA
diagnostic cache correction and updated installation/diagnostic documentation.
Only the private diagnostic cache helper changes executable behavior relative
to the preceding candidate; estimators, distances, data and APIs are unchanged.
This file records local validation, not a completed CRAN submission.

## Test environment

- Local macOS Sequoia 15.7.3, Apple silicon, R 4.6.0.
- Serial native build using a validation-only local Fortran-path correction.
- No system/user configuration or dependency installations were changed.

## R CMD check results

Full source build and `R CMD check --as-cran`: 0 errors | 0 warnings | 3 notes.

The three notes concern the local validation environment:

1. The system clock could not be verified remotely.
2. HTML validation was skipped because the installed HTML Tidy is too old.
3. Apple's build tools created an `xcrun_db` entry in the temporary directory.

Remote incoming checks were disabled explicitly; repository-index access
attempts also failed in the restricted environment. These local results do
not establish successful internet-enabled incoming checks.

The host R installation's default Fortran library directory does not exist.
Validation therefore reused a known working compiler/library configuration
outside the package. Distributed Makevars contain no host-specific paths.
This is not an unqualified clean-host or cross-platform portability result.

## Tests and documentation

- 753 installed test assertions passed, with zero failures and warnings.
- Five development-only R-reference tests were skipped because those reference
  files are deliberately outside the distributed package. They were not rerun
  for this diagnostic-only fix; their prior validation is historical evidence.
- The 101 new targeted assertions passed after reproducing 56 failures against
  the preceding implementation. They cover compact/noncompact kernels,
  cached/full/explicit weights, adaptive/fixed bandwidths, row standardization,
  and global/local permutation results on small deterministic examples.
- All four vignettes built and rebuilt successfully.
- Help examples, documentation consistency, and the PDF manual build passed.
- The shared README/CFF/BibTeX citation generator passed its synchronization check.

## Gates still pending

- Verified public arXiv identifiers for the primary LHWR and additional Lasso
  papers, to be added only after the corresponding records exist.
- Windows/Linux and current R-devel checks of the exact final source archive.
- Internet-enabled CRAN incoming checks, including current package-name checks.
- Explicit CRAN-submission authorization. GitHub publication approval does not
  authorize CRAN submission, arXiv deposit, or journal submission.
