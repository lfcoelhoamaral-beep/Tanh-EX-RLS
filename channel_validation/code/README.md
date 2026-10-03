# Reviewer 1 — separate channel parameter validation

Run from this extracted folder in MATLAB:

```matlab
reviewer1_channel_validation('smoke')
reviewer1_channel_validation('validate')
reviewer1_channel_validation('test')
```

Run smoke first and inspect its output. Full validation: 24 candidates/method, 30 realizations/scenario. Independent test: 200 realizations/scenario, selected configurations frozen. Four scenarios: alpha=2,1.5,1.2,1; sigma=0.1. The search budget matches the FIR protocol in configurations and realizations, not elapsed time or parameter dimension.

Outputs are in results_channel_smoke or results_channel_full. Do not confuse these with the FIR results_full directory. Calling without arguments defaults to smoke. Re-running a stage overwrites its corresponding files; no automatic checkpoint resumption is implemented. Test requires results_channel_full/selection.mat from full validation.

## Model and metric
The model reproduces the supplied rayleigh_tracking.m construction: five real-valued coefficients obtained from the real projection of a complex Gaussian AR(1) process; rho=J0(2*pi*fD*Ts), Ts=0.8 microseconds, fD=100 to 50 Hz at sample 1001; N=2000. This matches one-step correlation only, not the complete Jakes model. A=0.9999999368*I is fixed.

Desired output: d=tanh(u'*h+v). Noise is applied BEFORE saturation. Performance is (d-w_previous'*u)^2, measured BEFORE the update. Median across realizations is converted to dB, then averaged over 901:1000, 1001:1200 and 1901:2000. This evaluates saturated noisy output prediction, not recovery of the physical complex channel. The channel metric contains observation noise and nonlinear modeling error.

Selection score gives equal weights to all four scenarios and all three windows: 12 window averages. The Doppler change is not an abrupt coefficient jump; the transient window does not by itself measure recovery from a jump. RecoverySamples is provided mechanically using the same +3dB/20-sample definition and must be interpreted cautiously for this continuously varying channel.

## Data separation
Validation seed base 410000; test base 1010000. Per-scenario/realization seed: base+10000*scenario+run. Candidate seed 83017. All candidates and methods share the same generated data within a stage. Original caller RNG is restored. Candidate generation is separate from data generation. Test never selects parameters.

## Candidate ranges and baseline
Original channel settings are included as candidate 1: q=3.26e-7, r0=.01, initialP=1000 for EX methods; gamma=4 for Tanh; beta=.995, j=2 for EX-RNQ; RLM lambda=lambdaSigma=.99, Nw=14, initialP=20.

Remaining candidates use log-uniform q in [1e-8,1e-4], r0 in [.001,.1], initialP in [10,1e4], gamma in [.25,16], where applicable. EX-RLS and Tanh share nuisance-parameter triples. EX-RNQ additionally samples beta uniformly in [.95,.9995], j uniformly from {2,3}. RLM samples lambda/lambdaSigma uniformly in [.95,.9995] and Nw from {8,14,20,32}. Fixed RLM threshold multipliers retain the supplied implementation. Unused parameters are NaN in candidate_parameters.csv. This finite search is not a guarantee of global optimality.

Any exception/nonfinite trajectory disqualifies a candidate; failures are reported and not silently excluded from medians. This is a numerical screen, not a stability proof. Ties select the first candidate. Checkpoints are saved after each validation candidate. Independent test prints progress per method.

## Files to return
After validation, send results_channel_full/candidate_parameters.csv and validation_scores.csv. After independent testing, send results_channel_full/independent_test_summary.csv. Retain selection.mat and execution logs on your computer for reproduction.

The new driver has not been executed in MATLAB here. It preserves the supplied three helper files unchanged and retains the original channel equations and a priori metric. It does not modify earlier manuscript tables until author-executed results are reviewed.
