Version 1.0.0 -------------------------------------------------------------------
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


