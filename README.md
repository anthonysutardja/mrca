#CS176 HMM Project

## Background
This is a parallized implementation of parameter estimation of a hidden markov model using the Baum-Welch algorithm. There is a special harness `mrca.py` to utilize the parameter estimation for estimating the TMRCA problem.

## Running MRCA
You can run EM on a new set of initial parameters and a FASTA file with the following:

    python mrca.py initial_params.txt fasta_sequence.txt output_estimated_params output_initial_decodings output_estimated_decodings output_likelihood

This will run 15 iterations of EM and output the corresponding files. If you wish to run more iterations or specify a threshold, see the call to `tmrca.estimate_params(max_iter=15)` in `mrca.py`. You can specify either a threshold for the change in log-likelihood or a maximum iteration (or both!).

    Example:
    tmrca.estimate_params(thresh=1e-3, max_iter=100)



## Code Overview
### `class TSequence`
`TSequence` is a wrapper abstraction for performing initial decodings, parameter estimation, and estimated decodings for the TMRCA problem. Every `TSequence` is associated with an initial parameter and observed difference sequence. The class will automatically compute the observed differences.

##### `TSequence.estimate_params`
This function must be called before calling `decode_estimate` and/or accessing the the following class attributes: `estimate`, `likelihood`, and `initial_likelihood`.

### `class EM`
Generically performs parameter estimation for a given initial parameter and observed sequence. The code has been parallelized to calculate transition and emission probabilities. 

### `class Decoding`
Wrapper for calculating all decodings: Viterbi, posterior decoding, and mean expectation.

### `plot` module
Utilizes Matplotlib to plot the estimated decodings for the `mu`, `2mu`, `4mu`, and `5mu` samples. Call `plot_initial` or  `plot_estimate` with `mu`, `2mu`, `4mu`, or `5mu` to see the plots.

    Example
    $ python -i plot.py
    >>> plot_estimate('mu')

**NOTE:** You must have `matplotlib` installed in order to plot.

## Solutions to PS5 (estimates)
Our estimations and plots are stored in the `/estimates/` folder.

##Authors
- Anthony Sutardja
- Kevin Tee
