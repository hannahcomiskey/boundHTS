Version 1.0.2 -------------------------------------------------------------------

Updates following CRAN submission 10-03-2026

1. Description file updated to included a more useful description of the package.
2. References for the methods will be available soon, and the description file 
will update accordingly.
3. Examples added to exported functions.
4. message() used instead of cat() in `rBeta_4p`.
4. The package uses functions from the ExtDist package via Imports and explicit 
namespace calls (ExtDist::).
The function `dBeta_4p` has been removed. Therefore, no code has been copied or 
modified from ExtDist, and its authors are not listed in Authors@R.
5. Wrapper functions (`Poisson_convolution`, `Beta_convolution`, 
`ZIB_convolution`, and `ZOIB_convolution`) added to streamline user experience 
for generating convolutions.
6. Vignettes are updated to use these new wrapper functions. 

> devtools::check(cran=T)
── R CMD check results ──────────────────────────────────────────────────────────────────────── boundHTS 1.0.2 ────
Duration: 1m 1.6s

Version 1.0.1 -------------------------------------------------------------------
This is a new submission.

All checks were performed using R CMD check --as-cran.
No ERRORs or WARNINGs were found.

The package has been tested on Windows, macOS, and Linux via GitHub Actions.

Thank you for your time.

> devtools::check(cran = T)
── R CMD check results ─────────────────────────────────────────────────────────────────────────────────────── boundHTS 1.0.0 ────
Duration: 23.5s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔
> devtools::test()
ℹ Testing boundHTS
✔ | F W  S  OK | Context
✔ |          7 | Beta_convolution_density
✔ |          9 | dBeta_4p
✔ |          9 | dZOIB_4p
✔ |          7 | moment_condition_tilting
✔ |          7 | rBeta_4p
Error: Dimension mismatch across inputs. Check data, phi and weights.
✔ |         12 | rZOIB_4p
✔ |          9 | tilted_density_cont
✔ |          9 | tilted_density_discrete
✔ |          8 | ZOIB_convolution_density

══ Results ═══════════════════════════════════════════════════════════════════════════════════════════════════════════════════════
[ FAIL 0 | WARN 0 | SKIP 0 | PASS 77 ]


