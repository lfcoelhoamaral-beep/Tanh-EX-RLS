# Reviewer 1 — reproducible FIR parameter validation

This is NEW code prepared for the revision. It has not been executed in MATLAB in this environment. It does not change or replace the submitted results. The three supplied helper functions are included unchanged. Requires MATLAB and Signal Processing Toolbox (fir1).

## Run in MATLAB
Extract this package into a new folder, open that folder as MATLAB's Current Folder and run:

```matlab
reviewer1_parameter_validation('smoke')
```
Smoke uses three configurations per method and three realizations in two scenarios. It retains the full 2000-sample trajectory and the original measurement windows. Its selections are not suitable for the manuscript.

If smoke finishes without an execution error, run:

```matlab
diary('parameter_validation_log.txt')
version
ver
reviewer1_parameter_validation('validate')
diary off
```

The full validation evaluates 24 configurations per method, 30 realizations per scenario and nine unique scenarios (union of the alpha and sigma sweeps). It may take substantial time; no measured duration is available. Progress is printed and a checkpoint is saved after each candidate. Automatic checkpoint resumption is not implemented. Re-running overwrites the output files in the relevant results directory.

After inspecting the selected parameters and failure counts, run the independent test:

```matlab
diary('independent_test_log.txt')
reviewer1_parameter_validation('test')
diary off
```

The test reads frozen selections from results_full/selection.mat and runs 200 realizations per scenario on independently seeded data. It never reselects parameters from test performance. Test mode cannot load smoke selections.

## Selection protocol
- Same generated inputs, true systems and disturbances for every method/configuration within each stage.
- Validation seeds: 310000+10000*scenario+run. Independent test: 910000+10000*scenario+run. Both differ from the submitted seed 2026.
- Candidate seed: 73017; original caller RNG restored when function exits.
- Scenario rows: (2,.1), (1.8,.1), (1.5,.1), (1.2,.1), (1,.1), (1.5,.02), (1.5,.05), (1.5,.2), (1.5,.3).
- Same search budget of 24 configurations per method, including the submitted configuration as candidate 1. Budget equality means number of configurations/realizations; it does not mean equal wall-clock time or identical search dimensionality.
- EX-RLS and Tanh-EX-RLS share the same randomly sampled q, r0 and initial uncertainty values for corresponding candidates; Tanh additionally samples gamma.
- Scores are the equal average over nine scenarios and three window statistics (901:1000, 1001:1200, 1901:2000). Each window statistic averages the dB curve AFTER the median over realizations. Lower is better.
- Any execution exception or nonfinite trajectory disqualifies that candidate (score Inf). Failures are reported, never dropped from the median. These checks do not prove positive definiteness or mathematical stability.
- Ties select the first candidate. One configuration per method is frozen across all alpha and sigma scenarios.
- The score and search ranges are specified before running the new validation. Do not adjust them in response to independent test results.

## Search ranges

| Method | Tuned parameters |
|---|---|
| EX-RLS | q: 1e-8 to 1e-4; r0: 1e-3 to 1e-1; initialP: 10 to 1e4 |
| Tanh-EX-RLS | Same q/r0/initialP ranges; gamma: 0.25 to 16 |
| EX-RNQ | beta: 0.95 to 0.9995; q: 1e-8 to 1e-4; j: 2 or 3; initialP: 10 to 1e4 |
| RLM | lambda/lambdaSigma: 0.95 to 0.9995; Nw: 8,14,20,32; initialP: 10 to 1e4 |

Positive scale parameters use log-uniform random sampling; forgetting factors use uniform sampling. These ranges are declared revision design choices, not empirically demonstrated universal optima. State-transition matrix A=I, true systems, observation scenarios and RLM Hampel threshold multipliers are fixed. EX-RLS has no forgetting factor in the supplied state-space recursion. Gamma is unused by EX-RLS; r0 is unused by the supplied EX-RNQ implementation.

All methods initialize w=0. Submitted initial matrices (1000I for the EX methods, 20I for RLM) are included in candidate 1; random candidates tune the common initial matrix scale.

## Outputs to return for review
After smoke: console output and any error message.
After full validation: results_full/candidate_parameters.csv, validation_scores.csv and selection.mat.
After independent test: results_full/independent_test_summary.csv and both execution logs.

Validation checkpoint includes scores, failures, candidate definitions and the accumulated scenario rows. Selected and tested parameters are retained in MAT files.

## Limits and manuscript use
The search identifies the best tested candidate under the declared criterion. It does not establish the best achievable performance of any method, exhaust the parameter domain, or provide convergence/stability analysis. Thirty validation realizations provide a finite and potentially noisy selection estimate. Test results remain necessary.

This driver covers FIR identification only. The nonlinear AR(1) channel experiment has a different metric (output MedSE), true model and state settings and must be tuned/tested separately before making an article-wide claim of comparable tuning. The preliminary gamma=4 procedure remains undocumented by the supplied scripts; this NEW validation will provide a reproducible procedure going forward.

Do not state in the reviewer response that the search has been executed until logs and outputs are available. Preserve the submitted tables; any changed tables must be explicitly identified as revised experiments.
