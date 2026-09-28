## Submission

This is the first CRAN submission of metaweave, version 0.4.1.

The package reconstructs spatial ecological networks from species distributions
and supplied interaction models. It includes species-level probability and
stochastic-block backends, directed constrained network ensembles, summaries,
and spatial mapping. A small self-contained vignette uses synthetic inputs.
Publication case studies and their empirical inputs are not included.

## Test environments and results

* Windows 11 x64, R 4.6.1 (2026-06-24 ucrt):
  R CMD check --as-cran: 0 errors, 0 warnings, 1 note (New submission).
  PDF and HTML manual checks, examples, tests and vignette rebuild passed.
* Windows 11 x64, R-devel (2026-09-25 r90590 ucrt):
  R CMD check --as-cran with _R_CHECK_CRAN_INCOMING_=false: Status OK.
  Incoming checks were disabled only for this devel run after repeated SSL
  failures retrieving CRAN's archive index. Incoming checks passed on R 4.6.1.
  PDF/HTML manuals, examples, tests and vignette rebuild passed.
  This run used the same recorded R-4.6 binary dependency library as the release run.
* Windows 11 x64, R 4.5.1: package build, examples, tests and vignette passed.
  An additional full check passed the PDF manual; HTML validation was skipped
  in that earlier run because Tidy had not yet been installed.

The release and devel test runs each passed 131 expectations with no failures,
warnings or skips. No Linux or macOS results are claimed; these remain pending.

## Notes

The single note from the full R-release check reports a new submission.
The maintainer is Gabriel Munoz <gmunozacevedo@proton.me>, sole author and
copyright holder, affiliated with Concordia University. The license is MIT.

Official current CRAN, CRAN archive, and Bioconductor release/devel software
indices were checked on 2026-09-27; no case-insensitive metaweave match was found.
