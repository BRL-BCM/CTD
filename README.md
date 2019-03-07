# Connect The Dots (CTD):: A Generative Model for Modular Perturbation Detection
An R package for probabilistic estimation of multivariate feature sets against a partial correlation network.

## Installation
In R, install the devtools package, and install CTD by install_github(“BRL-BCM/CTD”).


## Look at the package vignette.
It will take you across all the stages in the analysis pipeline, including:

1. Background knowledge graph generation.
2. The encoding algorithm: including generating node permutations using a network walker, converting node permutations into bitstrings, and calculating the minimum encoding length between k codewords.
3. Calculate the probability of a node subset based on the encoding length.
4. Calculate similarity between two node subsets, using a metric based on mutual information.
