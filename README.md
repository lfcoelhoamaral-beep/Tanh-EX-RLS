# Tanh-EX-RLS

MATLAB implementations and reproducible experiments for **A Robust Tanh-Based Extended Recursive Least-Squares Algorithm Under alpha-Stable Noise**, by Luis Fernando Coelho Amaral (UFMA).

This repository contains the revised FIR and nonlinear channel experiments comparing EX-RLS, Tanh-EX-RLS, EX-RNQ and RLM. The manuscript is under revision; no acceptance or publication is claimed.

## Requirements

The author ran the experiments in MATLAB R2026a Update 2 on Windows 11. The FIR experiment requires Signal Processing Toolbox (`fir1`). The channel driver uses MATLAB functions and the supplied helpers.

## Reproduce the FIR experiments

Set MATLAB Current Folder to `FIR_validation/code` and execute in the Command Window:

```matlab
reviewer1_parameter_validation('smoke')
reviewer1_parameter_validation('validate')
reviewer1_parameter_validation('test')
```

The smoke stage is a small execution check and must not be used as paper evidence. Full validation evaluates 24 configurations per method using 30 runs in each of nine scenarios. It creates `results_full/selection.mat`. The test stage requires this file, freezes the selected settings and uses 200 independent runs per scenario.

## Reproduce the channel experiments

Set MATLAB Current Folder to `channel_validation/code` and execute:

```matlab
reviewer1_channel_validation('smoke')
reviewer1_channel_validation('validate')
reviewer1_channel_validation('test')
```

Full validation evaluates 24 configurations per method using 30 runs in each of four scenarios. The output folder is `results_channel_full`; `selection.mat` is required by the test stage, which uses 200 independent runs per scenario.

## Included evidence and outputs

Each experiment directory includes `candidate_parameters.csv`, `validation_scores.csv`, and `independent_test_summary.csv` from the completed author runs. Additional driver-specific documentation is in each `code/README.md`. Running the drivers creates output files inside the corresponding code directory; the included evidence is stored separately.

FIR performance measures coefficient error (MedSD); channel performance measures a priori output prediction error (MedSE). A median is taken across realizations, converted to dB, then averaged within each temporal window. These different metrics must not be pooled into a single ranking. Search ranges, seed rules and failure handling are declared in the drivers and local README files. Full precision is retained in the CSV candidate settings.

## Scope

Equal candidate counts are not exhaustive optimization or equal runtime budgets. No numerical failure in a finite test is a proof of stability. The theoretical guarantees use exact arithmetic and explicit restrictions; the noisy changing FIR case has no general infinite-horizon convergence guarantee. The channel is a real-output nonlinear prediction benchmark built from an AR(1) complex channel; it does not reproduce the full Jakes autocorrelation or certify physical channel recovery. Zero channel recovery values do not imply instantaneous recovery.

## Attribution

EX-RNQ: Amaral, Lopes and Barros, IEEE Signal Processing Letters (2018), doi:10.1109/LSP.2018.2864609. RLM: Zou, Chan and Ng, IEEE Signal Processing Letters (2000), doi:10.1109/97.873571. These are MATLAB implementations of the cited algorithms, not a claim of original authorship of the baseline methods.

## License

See the MIT LICENSE supplied at the repository root.
