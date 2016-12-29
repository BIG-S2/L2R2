Programs for L2R2 model in the manuscript "Bayesian longitudinal low-rank regression models for imaging genetic data from longitudinal studies"

L2R2 is a function for estimating L2R2 model with MCMC algorithm
sim1.m is an example that use the L2R2 function.

1. Specify the input of responses, genetic predictors, prognostic covariates, and covariates of random effects, and other constants
2. Specify the Hyperparameters for parameters in the mean structure
3. Specify the Hyperparameters for parameters in the loading structure
4. Specify the Intitial values of mean parameters
5. Specify the  Intitial values of Factor Model parameters
6. Run the sim1.m matlab script. 

The MCMC output are store in the struct "Output", which contains the MCMC samples, posterior mean and SE of the parameters.
The workspace is also stored in a mat file.

