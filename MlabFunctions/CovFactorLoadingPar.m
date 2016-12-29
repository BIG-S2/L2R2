function FactorModelPar = CovFactorLoadingPar(Y,k,HyperParameters,FactorModelPar)
%%
[n,d] = size(Y);
% --- Define hyperparameter values --- %
as = HyperParameters.as; bs = HyperParameters.bs;    % gamma hyperparameters for diagonal elements of inv(Sigma)
df = HyperParameters.v;            % gamma hyperparameters for t_{ij}
ad1 = HyperParameters.a1; bd1 = 1;  % gamma hyperparameters for delta_1
ad2 = HyperParameters.a2; bd2 = 1;  % gamma hyperparameters delta_h, h >= 2

% --- Initial values --- %
sig = FactorModelPar.sig;                   % diagonals of sigmainv
Lambda = FactorModelPar.Lambda;                       % loadings matrix
eta = FactorModelPar.eta;        % latent factors
t = FactorModelPar.t;               % local shrinkage coefficients
delta = FactorModelPar.delta; % gobal shrinkage coefficients multilpliers
tau = FactorModelPar.tau;                      % global shrinkage coefficients
D = bsxfun(@times,t,tau');


%%
%-----Step 1: Update Lambda-----%
etasq = (eta'*eta);
%         Lambda = zeros(d,k);
for j = 1:d
    idx = D(j,:) < 1E16;
    if sum(idx)>0
        Vlam1 = diag(D(j,idx)) + sig(j)*etasq(idx,idx);
        Vlam = Vlam1\eye(size(Vlam1));%Vlam = S*S';
        Elam = Vlam*sig(j)*eta(:,idx)'*Y(:,j);
        Lambda(j,idx) = mvnrnd(Elam',Vlam);
                
    end
    Lambda(j,~idx) = zeros(1,sum(~idx));
end
%%
% -- Step 2: Update Sigma -- %
sig = gamrnd(as + n/2,1./(bs+0.5*sum((Y - eta*Lambda').^2)))';
%Sigma = diag(1./sig);

% -- Step 3: Update eta -- %
Lmsg = bsxfun(@times,Lambda,sig);
Veta1 = eye(k) + Lmsg'*Lambda;
Veta = Veta1\eye(size(Veta1));%S*S';                     % Veta = inv(Veta1)
Meta = Y*Lmsg*Veta;
eta = mvnrnd(Meta,Veta);        % update eta in a block

%------Step 4: Update tij's = phi_ijs ----------%
t = gamrnd(df/2 + 0.5,1./(df/2 + bsxfun(@times,Lambda.^2,tau')));

%------Step 5: Update delta----------%
bd = bd1 + 0.5* (1/delta(1))*sum(tau'.*(sum(t.*Lambda.^2)));
delta(1) = gamrnd(ad1 + d*k/2, 1/bd);

tau = cumprod(delta);

tLmsq = sum(t.*Lambda.^2);
for h = 2:k
    ad = ad2 + d*(k-h+1)/2;
    temp1 = tau'.*tLmsq;
    bd = bd2 + 0.5* (1/delta(h))*sum(temp1(h:k));
    delta(h) = gamrnd(ad,1/bd);
    tau  = cumprod(delta);
end


disp('Column with Max Absolute Value > 0.0001');
disp(sum(max(abs(Lambda),[],1)>0.0001));


%-------update precision parameters--------%
O = diag(sig) - Lmsg*Veta*Lmsg';
O = (O+O')/2;


FactorModelPar.sig = sig;                   % diagonals of sigmainv
FactorModelPar.Lambda = Lambda;                       % loadings matrix
FactorModelPar.eta =  eta;        % latent factors
FactorModelPar.t = t;               % local shrinkage coefficients
FactorModelPar.delta = delta; % gobal shrinkage coefficients multilpliers
FactorModelPar.tau = cumprod(delta);                      % global shrinkage coefficients
FactorModelPar.O = O;

end