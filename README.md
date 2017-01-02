# L2R2
Programs for L2R2 model in the manuscript "Bayesian longitudinal low-rank regression models for imaging genetic data from longitudinal studies"

This L2R2 package is developed by Zhao-Hua Lu, Zakaria Khondker, and Hongtu Zhu from [BIG-S2 lab](http://odin.mdacc.tmc.edu/bigs2/). The package is also available on [BIG-S2 GitHub site](https://github.com/BIG-S2).



To perform a joint analysis of multivariate  neuroimaging phenotypes and candidate genetic markers obtained from longitudinal  studies, we develop a Bayesian longitudinal low-rank regression (L2R2) model. The L2R2 model integrates three key methodologies: a low-rank matrix for approximating the high-dimensional regression coefficient matrices corresponding to the genetic main effects and their interactions with time, penalized splines for characterizing the overall time effect, and a sparse factor analysis model coupled with random effects for capturing within-subject spatio-temporal correlations of longitudinal phenotypes. Posterior computation proceeds via an efficient Markov chain Monte Carlo algorithm.

L2R2 is a function for estimating L2R2 model with MCMC algorithm.

sim1.m is an example that uses the L2R2 function.

1. Specify the input of responses, genetic predictors, prognostic covariates, and covariates of random effects, and other constants.
2. Specify the Hyperparameters for parameters in the mean structure.
3. Specify the Hyperparameters for parameters in the loading structure.
4. Specify the Intitial values of mean parameters.
5. Specify the  Intitial values of Factor Model parameters.
6. Run the sim1.m matlab script. 

The MCMC output are store in the struct "Output", which contains the MCMC samples, posterior mean and SE of the parameters.
The workspace is also stored in a mat file.
