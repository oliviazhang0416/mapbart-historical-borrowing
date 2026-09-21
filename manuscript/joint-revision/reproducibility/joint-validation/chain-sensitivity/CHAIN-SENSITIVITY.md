# Saved-chain constant-effect sensitivity

This is a descriptive chain-index sensitivity check from the reviewed saved-trace summary. Each chain-specific RMSE uses all 120 constant-effect method--dataset groups: 40 Compatible, 40 One-variable and 40 Two-variable. No ATE convergence flags or other diagnostics were used to exclude rows.

Input: `05-writing/journal-revision/reproducibility/joint-validation/chain-sensitivity/INPUT-constant_chain_rmse.csv` (SHA-256 `59d8a7bccd019b8b0ece7b968e6791237f549c109d6383c64624a33304dd5c17`). The upstream diagnostic follow-up that produced this input read the saved 1,800 raw RDS results; this derivation performs no MCMC and no refit.

The percentage is $100(1-\mathrm{RMSE}_{\mathrm{LRC}}/\mathrm{RMSE}_{\mathrm{comparator}})$. Chain 0 in the upstream file is the existing pooled three-chain reference; the table below reports only individual saved chains 1--3.

| Chain | Groups | Joint LRC RMSE | Joint source RMSE | Separated source RMSE | LRC vs separated | LRC vs joint source |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 120 | 0.199402950331 | 0.214420818863 | 0.209068948516 | 4.623354% | 7.003923% |
| 2 | 120 | 0.199120605479 | 0.214919910289 | 0.209618300883 | 5.008005% | 7.351252% |
| 3 | 120 | 0.201897059958 | 0.215474179897 | 0.210701071381 | 4.178437% | 6.301043% |

The chain-index reductions range from 4.178% to 5.008% against Separated source-BART and from 6.301% to 7.351% against Joint source-BART. This narrows the descriptive sensitivity of the pooled RMSE direction across saved chain indices; it is not a convergence proof, a replacement for the ATE diagnostic screen, or evidence that the chains have mixed adequately.

The upstream frozen raw files, hashes and validation outputs remain unchanged.
