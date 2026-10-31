# Working on `celestial`

`celestial` is an R package of astronomical utilities. Most functionality is in
`R/`, package metadata is in `DESCRIPTION` and `NAMESPACE`, and user-facing
function documentation is in `man/`. The `src/` directory contains the
Rcpp-backed implementation and generated Rcpp bindings.

## Tests and checks

- Add package tests under `tests/testthat/`, with descriptive
  `test-<area>.R` filenames. Use `test_that()` and `expect_*()` from testthat.
- Run the focused suite with `Rscript -e 'devtools::test()'` and the package
  check with `R CMD check .` when the required R dependencies are available.
- Keep tests deterministic; avoid network access and uncontrolled randomness.
- Keep `Suggests` and CI test dependencies in sync when adding test tools.

## Changes

- Preserve existing public function names, arguments, units, and return
  structures unless a change explicitly calls for an API change.
- Update the relevant Rd documentation when changing a public function.
- Do not hand-edit `R/RcppExports.R` or `src/RcppExports.cpp`; regenerate these
  bindings with Rcpp tooling when the C++ interface changes.
- Keep changes focused, and add tests for changed behavior and boundary cases.
