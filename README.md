# Connect The Dots (CTD):: a method which interprets multivariate perturbations identified in molecular profiles by identifying highly connected nodes in disease-specific co-perturbation networks
An R package for probabilistic estimation of multivariate feature sets against a partial correlation network.

## Check out CTD via Google Colab (UNDER CONSTRUCTION)
Check it out at https://tinyurl.com/BRL-BCM-CTD.
Here you can reproduce analysis results from Thistlethwaite et al. (2019).

## Using CTD in R.
### Installation
In R, install the devtools package, and install CTD by install_github(“BRL-BCM/CTD”).

### Look at the package Rmd vignette.
Located in /vignette/CTD_Lab-Exercise.html. It will take you across all the stages in the analysis pipeline, including:

1. Background knowledge graph generation.
2. The encoding algorithm: including generating node permutations using a network walker, converting node permutations into bitstrings, and calculating the minimum encoding length between k codewords.
3. Calculate the probability of a node subset based on the encoding length.
4. Calculate similarity between two node subsets, using a metric based on mutual information.

## References
Thistlethwaite, L.R. et al. (2019). CTD: a method which interprets multivariate perturbations identified in molecular profiles by identifiying highly connected nodes in disease specific co-perturbation networks. In preparation.
