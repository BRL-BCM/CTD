# CTD: an information-theoretic method to interpret multivariate perturbations in the context of graphical models with applications in metabolomics and transcriptomics
Our novel network-based approach, CTD, “connects the dots” between metabolite perturbations observed in individual metabolomics profiles and a given disease state by calculating how connected those metabolites are in the context of a disease-specific network.

## Using CTD in R.
### Installation
In R, install the devtools package, and install CTD by install_github(“BRL-BCM/CTD”). 

### Look at the package Rmd vignette.
Located in /vignette/CTD_Lab-Exercise.Rmd. It will take you across all the stages in the analysis pipeline, including:
1. Background knowledge graph generation.
2. The encoding algorithm: including generating node permutations using a network walker, converting node permutations into bitstrings, and calculating the minimum encoding length between k codewords.
3. Calculate the probability of a node subset based on the encoding length.
4. Calculate similarity between two node subsets, using a metric based on mutual information.

## References
Thistlethwaite L.R., Petrosyan V., Li X., Miller M.J., Elsea S.H., Milosavljevic A. (2020). CTD: an information-theoretic method to interpret multivariate perturbations in the context of graphical models with applications in metabolomics and transcriptomics. Manuscript in review.
