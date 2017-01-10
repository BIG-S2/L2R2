path(path,'./MlabFunctions/');

%% Input 

Y0 = dlmread(['./TestData/Y_1.txt'],'\t');Y0(:,size(Y0,2))=[];% N*d matrix of response variable 
X = dlmread(['./TestData/X_1.txt'],'\t');X(:,size(X,2))=[]; % N*p matrix of genetic predictors
W = dlmread(['./TestData/W_1.txt'],'\t');W(:,size(W,2))=[]; % N*q matrix of prognostic covariates, NaN is no prognostic factor is used.
Z = dlmread(['./TestData/Z_1.txt'],'\t');Z(:,size(Z,2))=[]; % N*(n*p_b) matrix of covariates of random effects, NaN is no prognostic factor is used.

TypeCov='Factor'; % Don't change
r = 3; % Rank of the low rank representation of the coefficient matrix
MCMCpara.FacModR = 25; % factor model number of factors
OutputFileName = 'Result';  %filename of output mat file

 
MCMCpara.nBurnin = 4e0; % number of burn-in samples
MCMCpara.nCollect = 4e0; % number of samples after burn-in

%% Hyperparameters for parameters in the mean structure
MeanStructHyperpar.a0 = 1e-6; % 1st gamma hyperparameter of tau_delta
MeanStructHyperpar.b0 = 1e-6; % 2nd gamma hyperparameter of tau_delta
MeanStructHyperpar.c0 = 1e-6; % 1st gamma hyperparameter of tau_gamma
MeanStructHyperpar.d0 = 1e-6; % 2nd gamma hyperparameter of tau_gamma
MeanStructHyperpar.e0 = 1e-6; % 1st gamma hyperparameter of tau_b
MeanStructHyperpar.f0 = 1e-6; % 2nd gamma hyperparameter of tau_b
MeanStructHyperpar.g0 = 1e-6; % Not used now
MeanStructHyperpar.h0 = 1e-6; % Not used now

%% Hyperparameters for parameters in the loading structure
FactorModelHyperpar.v=3;    % v/2 is the shape and rate hyperparameter of the gamma prior for local shrinkage parameter
FactorModelHyperpar.a1=1;   % the shape hyperparameter of the gamma prior for global shrinkage parameter
FactorModelHyperpar.a2=3;   % the rate hyperparameter of the gamma prior for global shrinkage parameter
FactorModelHyperpar.as=3;   % the shape hyperparameter of the gamma prior for the precision of error
FactorModelHyperpar.bs=1;   % the rate hyperparameter of the gamma prior for the precision of error

%% Intitial values of mean parameters
[N,p] = size(X); % Get dimensions
d = size(Y0,2); % dimension of reponses
q = size(W,2); % dimension of genetic predictor, e.g., SNPs and SNP-age interactions

MeanStructPar.U = randn(p,r)*0.01; % low-rank component
MeanStructPar.V = randn(d,r)*0.01; % low-rank component
MeanStructPar.Delta = randn(r,1)*0.01; % low-rank component
MeanStructPar.G = randn(q,d)*0.01; % polynomial and spline coefficients
MeanStructPar.b = randn(size(Z,2), d)*0.01; % random effect coefficients
MeanStructPar.s = 3; % spline degree
MeanStructPar.tau_del = 1;
MeanStructPar.tau_e = 1; % needed for iid case instead of O
MeanStructPar.tau_G = 1;
MeanStructPar.tau_b = 1;

%% Intitial values of Factor Model parameters
k = MCMCpara.FacModR;
FactorModelPar.sig = gamrnd(FactorModelHyperpar.as,1/FactorModelHyperpar.bs,d,1); % diagonals of sigmainv
FactorModelPar.Lambda = zeros(d,k);                       % loadings matrix
FactorModelPar.eta =  mvnrnd(zeros(1,k),eye(k),N);        % latent factors
FactorModelPar.t = gamrnd(3/2,2/3,[d,k]);               % local shrinkage coefficients
FactorModelPar.delta = [gamrnd(FactorModelHyperpar.a1,1);gamrnd(FactorModelHyperpar.a2,1,[k-1,1])]; % gobal shrinkage coefficients multilpliers
FactorModelPar.tau = cumprod(FactorModelPar.delta);         % global shrinkage coefficients

% Flags
Flags.OrthogonalZ = false; % Are the columes of Z orthogonal, Z is the covariate matrix of random effects
Flags.RecordMCSampleb = false; % Record the MCMC samples of random effects !! Large memory required
Flags.RecordMCSampleO = false; % Record the MCMC samples of precision matrix of error !! Large memory required


%% Automatical initialization

% Precision error matrix
Lmsg = bsxfun(@times,FactorModelPar.Lambda,FactorModelPar.sig);
FactorModelPar.O = diag(FactorModelPar.sig) - Lmsg*((eye(k)+FactorModelPar.Lambda'*Lmsg)\(Lmsg'));

% Put all in MCMCpara
MCMCpara.MeanStructHyperpar = MeanStructHyperpar;
MCMCpara.FactorModelHyperpar = FactorModelHyperpar;

MCMCpara.MeanStructPar = MeanStructPar;
MCMCpara.FactorModelPar = FactorModelPar; 

MCMCpara.Flags = Flags;

if isnan(W)
   W = zeros(size(Y0,1),0); 
end
if isnan(Z)
   Z = zeros(size(Y0,1),0); 
end


%% RUN L2R2 %%

t0 = cputime;

tsta1 = cputime;
rng(1);
Output = L2R2(X,W,Z,Y0,r,MCMCpara,TypeCov); %Wish, Diag, Factor, iid(default)
tend1 = cputime;
disp(tend1 - tsta1);


time = (cputime - t0)/60;
disp(['Time = ', num2str(time), ' minutes']);


%% %%%%%%% Saving %%%%%%%%%

dlmwrite(['./B0_1.txt'],Output.PostB,'\t')
dlmwrite(['./VB0_1.txt'],Output.Bsd,'\t')
save(['./',OutputFileName,'.mat'])
