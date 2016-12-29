function Output = L2R2(X,W,Z,Y0,r,MCMCpara,TypeCov) %#codegen
%Generalized Reduced Rank Regression Via Decomposition
%
% Model: Y0 = U*(Delta.*Z)*V + EB + Zb + E
%-----------------------------------------------------------------
%
%USAGE: Output = Bayesian Longitudinal Low Rank Regression results
%INPUT :
%   Y0: N x d, input data matrix.  every column is a rsponse from the same subject.
%   X: N x p, covariate matrix for low rank component.
%   W: N x q, prognostic factor and spline matrix for traditiona component.
%   Z: N x n, covariate matrix for random time coefficient.
%   r: scalar, the rank of the low rank decomposition of the coefficient
%   matrix
%   MCMCpara, input structure of initial values and hyperparameters, an
%   example script for specifying MCMCpara is sim1.R
%   TypeCov, A string, default to be "Factor", to be used later
%
%
%
%
%OUTPUT:
%   Output: a struct data.

% MCMC samples
% Output.U -- Coefficients of linear combinations of covariates
% Output.V -- Coefficients of linear combinations of responses
% Output.Delta -- Coefficients of diagonal elements
% Output.G -- Coefficients of prognostic factors
% Output.b -- random effects
% Output.O -- Precision matrix of error
% 
% 
% Estimated posterior mean
% Output.PostU 
% Output.PostV 
% Output.PostDelta
% Output.PostB -- coefficients between SNPs and ROIs
% Output.PostG 
% Output.Postb 
% Output.PostO 
% 
% Estimated SE
% Output.Usd 
% Output.Vsd 
% Output.Deltasd 
% Output.Osd 
% Output.Bsd 
% Output.Gsd 
% Output.bsd 

% ---------------------
% check input arguments
% ---------------------

%% get dimensions
[N,p] = size(X);
d = size(Y0,2);
q = size(W,2);
n = size(Z,2); % total dimension of the random effects

nBurnin = MCMCpara.nBurnin;
nCollect = MCMCpara.nCollect;
FacModR = MCMCpara.FacModR;


%% Hyperparameters for parameters in the mean structure
MeanStructHyperpar = MCMCpara.MeanStructHyperpar;
a0 = MeanStructHyperpar.a0;
b0 = MeanStructHyperpar.b0;
c0 = MeanStructHyperpar.c0;
d0 = MeanStructHyperpar.d0;
e0 = MeanStructHyperpar.e0;
f0 = MeanStructHyperpar.f0;
g0 = MeanStructHyperpar.g0;
h0 = MeanStructHyperpar.h0;

%% Hyperparameters for parameters in the loading structure
FactorModelHyperpar = MCMCpara.FactorModelHyperpar;

%% Intitial values of mean parameters
MeanStructPar = MCMCpara.MeanStructPar;
U = MeanStructPar.U; % low-rank component
V = MeanStructPar.V; % low-rank component
Delta = MeanStructPar.Delta; % low-rank component
G = MeanStructPar.G; % polynomial and spline coefficients
b = MeanStructPar.b; % random effect coefficients
s = MeanStructPar.s; % for cubic spline
tau_del = MeanStructPar.tau_del;
tau_e = MeanStructPar.tau_e; % needed for iid case instead of O
tau_G = MeanStructPar.tau_G;
tau_b = MeanStructPar.tau_b;

%% Factor Model parameters
FactorModelPar = MCMCpara.FactorModelPar;
O = FactorModelPar.O;

%% get flags
OrthogonalZ = MCMCpara.Flags.OrthogonalZ;
RecordMCSampleb = MCMCpara.Flags.RecordMCSampleb;
RecordMCSampleO = MCMCpara.Flags.RecordMCSampleO;



%% pre-computing some matrix products
XX = X'*X;
WW = W'*W;
ZZ = Z'*Z;
Id = eye(d);
Ip = eye(p);
Iq = eye(q);
In = eye(n);

ZZ2 = sum(Z.^2,1);
[svdU0, svdD0, ~] = svd(X');

svdD = svdD0;
D2 = diag(svdD).^2;
svdU = svdU0;
svdUt = svdU';

%% Allocate output arrays;
SumU = zeros(p,r);
SumSU = zeros(p,r);
SumV = zeros(d,r);
SumSV = zeros(d,r);
SumDelta = zeros(r,1);
SumSDelta = zeros(r,1);
SumB = zeros(p,d);
SumSB = zeros(p,d);
SumG = zeros(q,d);
SumSG = zeros(q,d);
Sumb = zeros(n,d);
SumSb = zeros(n,d);
SumO = zeros(d,d);
SumSO = zeros(d,d);


Output.U = zeros(p,r,nCollect);
Output.V = zeros(d,r,nCollect);
Output.Delta = zeros(r,nCollect);
Output.G = zeros(q,d,nCollect);

nCollectb = 1;
nCollectO = 1;

if RecordMCSampleb
    nCollectb = nCollect;
end
if RecordMCSampleO
    nCollectO = nCollect;
end
Output.b = zeros(n, d,nCollectb);
Output.O = zeros(d,d,nCollectO);

Output.PostU = zeros(p,r);
Output.PostV = zeros(d,r);
Output.PostDelta = zeros(r,1);
Output.PostB = zeros(p,d);
Output.PostG = zeros(q,d);
Output.Postb = zeros(n, d);
Output.PostO = zeros(d, d);

Output.Usd = zeros(p,r);
Output.Vsd = zeros(d,r);
Output.Deltasd = zeros(r,1);
Output.Osd = zeros(d,d);
Output.Bsd = zeros(p,d);
Output.Gsd = zeros(q, d);
Output.bsd = zeros(n, d);

Output.nBurnin = nBurnin;
N0 = nCollect;
Output.nCollect = N0;



%%
for iter = 1:(nBurnin + nCollect)
    
    %%  Sampling the B = U \Delta V
    B = U*diag(Delta)*V';
    Y = Y0 - X*B - W*G - Z*b;  % subtracting all the layers
    %     %---------------------------------------------------------------
    %     % sampling by LAYER
    %     %---------------------------------------------------------------
    for l = 1:r
        Y = Y + X*Delta(l)*U(:,l)*V(:,l)'; % adding back the l-th layer
        XY = X'*Y;
        
        %Sample U / U----------------------------------------------------------
        Dn2 = 1./(p+(Delta(l)).^2*(V(:,l)'*O*V(:,l)).*D2);
        U(:,l) = svdU*((Delta(l)).*Dn2.*(svdUt*(XY*O*V(:,l))) + sqrt(Dn2).*randn(p,1));
        
        %Vample V / V ----------------------------------------------------------
        
        sigVl = (d*Id + (Delta(l)).^2*U(:,l)'*XX*U(:,l)*O)\eye(size(Id));%inv(d*Id + (Delta(l)).^2*U(:,l)'*XX*U(:,l)*O);
        muVl = Delta(l)*(sigVl*(O*XY'*U(:,l)))';
        V(:,l) = mvnrnd(muVl, sigVl);
        
        %Sample Delta(l)
        sig_Delta = 1/(tau_del + (U(:,l)'*XX*U(:,l))*(V(:,l)'*O*V(:,l)));
        mu_Delta = sig_Delta*U(:,l)'*XY*O*V(:,l);
        Delta(l) = normrnd(mu_Delta, sqrt(sig_Delta));
        %------------------------------------------------------------------
        
        Y = Y - X*Delta(l)*U(:,l)*V(:,l)'; % subtracting the l-th layer again
    end
    B = U*diag(Delta)*V';
    
    %% Sampling spline coefficients G one columns at a time
    Y = Y + W*G;
    for k = 1:d
        Sig_Gk = (WW*O(k,k)+tau_G*Iq)\eye(size(Iq));%inv(WW*O(k,k)+tau_G*Iq);
        G_k = DropColumn(G,k);
        Y_k = DropColumn(Y,k);
        Ok = SliceColumn(O,k);
        G(:,k) = (mvnrnd(Sig_Gk*W'*(Y(:,k)*O(k,k)+(Y_k - W*G_k)*Ok), Sig_Gk));
    end;
    Y = Y - W*G;
    
    
    %% Sampling random effect b one columns at a time
    Y = Y + Z*b;
    
    if ~OrthogonalZ
        for k = 1:d
            Sig_bk = (ZZ*O(k,k) + tau_b*In)\eye(size(In));%inv(ZZ*O(k,k) + tau_b*In);
            b_k = DropColumn(b,k);
            Y_k = DropColumn(Y,k);
            Ok = SliceColumn(O,k);
            b(:,k) = (mvnrnd(Sig_bk*Z'*(Y(:,k)*O(k,k)+(Y_k - Z*b_k)*Ok), Sig_bk));
        end;
    else
        for k = 1:d
            Sig_bk_diag = 1./(ZZ2*O(k,k) + tau_b)';
            b_k = DropColumn(b,k);
            Y_k = DropColumn(Y,k);
            Ok = SliceColumn(O,k);
            vtmp1 = Z'*(Y(:,k)*O(k,k)+(Y_k - Z*b_k)*Ok);
            b(:,k) = vtmp1 .* Sig_bk_diag + randn(n,1).*sqrt(Sig_bk_diag);
        end
    end
    
    E = Y - Z*b;
    
    %% Sample the variances in the prior distributions of Delta, G, b
    tau_del = gamrnd(a0+r/2,1./(b0 + sum(Delta.^2)/2));
    tau_G = gamrnd(c0+0.5*q*d,1./(d0 + sum(sum(G.^2))/2));
    tau_b = gamrnd(e0+0.5*n*d,1./(f0 + sum(sum(b.^2))/2));
    
    
    
    %% sampling parameters in the factor model for the covariance/precision matrix
    FactorModelPar = CovFactorLoadingPar(E,FacModR,FactorModelHyperpar,FactorModelPar);
    O = FactorModelPar.O;
    
    %% Record samples
    % Collect samples
    if iter > nBurnin
        idx = ceil(iter-nBurnin);
        Output.U(:,:,idx) = U;
        Output.V(:,:,idx) = V;
        Output.Delta(:,idx) = Delta;
        Output.G(:,:,idx) = G;
        if RecordMCSampleb
            Output.b(:,:,idx) = b;
        end
        if RecordMCSampleO
            Output.O(:,:,idx) =  O;
        end
        
        SumDelta = SumDelta + Delta;
        SumSDelta = SumSDelta + Delta.^2;
        SumB = SumB + B;
        SumSB = SumSB + B.^2;
        SumU = SumU + U;
        SumSU = SumSU + U.^2;
        SumV = SumV + V;
        SumSV = SumSV + V.^2;
        SumG = SumG + G;
        SumSG = SumSG + G.^2;
        Sumb = Sumb + b;
        SumSb = SumSb + b.^2;
        
        SumO = SumO + O;
        SumSO = SumSO + O.^2;
    end
    
    %     disp(iter);
    
end

%% Posterior mean ane SE of parameters, random effects, error precision matrix

Output.PostU = SumU/N0;
Output.PostV = SumV/N0;
Output.PostDelta = SumDelta/N0;
Output.PostB = SumB/N0;
Output.PostG = SumG/N0;
Output.Postb = Sumb/N0;
Output.PostO = SumO/N0;

Output.Usd = sqrt(SumSU/N0 - Output.PostU.^2);
Output.Vsd = sqrt(SumSV/N0 - Output.PostV.^2);
Output.Deltasd = sqrt(SumSDelta/N0 - Output.PostDelta.^2);
Output.Osd = sqrt(SumSO/N0 - Output.PostO.^2);
Output.Bsd = sqrt(SumSB/N0 - Output.PostB.^2);
Output.Gsd = sqrt(SumSG/N0 - Output.PostG.^2);
Output.bsd = sqrt(SumSb/N0 - Output.Postb.^2);


end


