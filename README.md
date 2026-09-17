# LRC-BART historical borrowing

Revision of the MAP-BART historical-borrowing work, renamed LRC-BART. This
branch supersedes the `mapbart-*` layout on `main`.

## Layout

| Path | Contents |
|---|---|
| `lrcBART/` | LRC-BART C++ engine (`clrcbart.cpp`, `cess.cpp`, `include.lrc/`) |
| `wBART/`, `aBART/` | Gaussian and AFT comparator engines |
| `bartModelMatrix.R` | Shared design-matrix helper |
| `lrcbart-sim-gaussian{,-single-arm}/` | Gaussian simulations |
| `lrcbart-sim-survival{,-single-arm}/` | Survival simulations |
| `lrcbart-case-study-mm/` | Myeloma application (EloKRd vs UCMM) |
| `manuscript/` | Paper, appendix, and build scripts |
| `*_MIGRATION_PLAN.md` | Plans governing the migration |

## Running

Each project is driven by its own `run_all.R`:

```bash
cd lrcbart-sim-gaussian
Rscript run_all.R
```

Scripts locate the repository root by walking up from the script's own
location (and, failing that, the working directory) until they find the
`.lrcbart-root` marker. You can therefore invoke them from any directory:

```bash
Rscript /path/to/lrcbart-historical-borrowing/lrcbart-sim-gaussian/run_all.R
```

No paths need editing.

Requires R with `Rcpp` and `RcppEigen`; the C++ engines compile on first use
via `Rcpp::sourceCpp`.

## Not in this repository

`res/`, `inserts/`, and the simulation `data/` directories hold regenerable
output and ship as empty placeholders. Rebuild simulated data with each
project's `data_gen_p10.R`, then run `run_all.R`.

`lrcbart-sim-gaussian/PSCL.R` needs **psrwe** (Wang et al.), a third-party
package that is not vendored here. Point `PSRWE_DIR` at a local checkout:

```bash
PSRWE_DIR=/path/to/psrwe Rscript run_all.R
```

Otherwise it is looked for beside this repository. `PSCL.R` sources that
checkout's `R/` directory directly, so it runs your local psrwe source rather
than an installed build.

`MAP.R` and `ess_local/ess_cal_map.R` use the **installed** `RBesT` package
(`library(RBesT)`), not a local copy. Results depend on the installed version,
so record it when reporting: `packageVersion("RBesT")`.

Every other comparator is included in this repository.

Patient-level inputs are excluded; `lrcbart-case-study-mm/data_cleaned/`
carries only the de-identified merged analysis sets.
