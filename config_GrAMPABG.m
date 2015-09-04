function [GrampaEstimIn,GrampaEstimOut_ss,GrampaLinTrans,GrampaOptions]=config_GrAMPABG(H,y,wvar,KM,iteration_stopping_tol,maxiter)
[M,N]=size(H);
%------------------------- Configuring GrAMPA ----------------------------%
% 1D TV dictionary 
I = speye(N-1);
P = spalloc(N-1,1,0);
scaling = sqrt(norm(H,'fro')^2/M); % to equalize row norms of H and D
D = scaling/sqrt(2) *( [P I] - [I P] ); % note sparse and properly scaled
% linear transform exploiting nonsparsity of H and sparsity of D
GrampaLinTrans = LinTransConcat({MatrixLinTrans(H);MatrixLinTrans(D)});


% options that seem to work well for 1D TV dictionaries
 GrampaOptions = GampOpt;
 GrampaOptions.legacyOut = false;
 GrampaOptions.xvar0 = (norm(y)^2-wvar*M)/norm(H,'fro')^2; % guess at signal variance
 xvar0big = 100*GrampaOptions.xvar0; % much larger than signal variance 
 
 GrampaOptions.tol = sqrt(iteration_stopping_tol); % stopping tolerance: sqrt(1e-7) used by ssAMP
 GrampaOptions.nit = maxiter; % maximum number of iterations
 GrampaOptions.varNorm = false; % turn off internal normalization
 GrampaOptions.zvarToPvarMax = inf; % do not clip 
 GrampaOptions.uniformVariance = true; % off since applied externally
 GrampaOptions.xvarMin = 0; % no minimum rvar
 GrampaOptions.pvarMin = 0; % no minimum pvar
%  GrampaOptions.pvarStep = true; % apply damping to pvar 
  GrampaOptions.pvarStep = false; % apply damping to pvar 
 GrampaOptions.adaptStep = false; % leave off 
 %-----------------------------------------------------------
%  GrampaOptions.adaptStepBethe = true; % use Bethe version 
%  GrampaOptions.stepWindow = 10; % adaptive stepsize window 
%  GrampaOptions.stepIncr = 1.1; % stepsize increase rate
%  GrampaOptions.stepDecr = 0.5; % stepsize decrease rate
 %--------------------------------------------------------
 GrampaOptions.stepMax = 1.0; % maximum stepsize: 1.0 for speed, 0.5 for robustness
 GrampaOptions.stepMin = 0.25; % minimum stepsize
 GrampaOptions.step = GrampaOptions.stepMax; % initial stepsize

% trivial prior
% GrampaEstimIn = NullEstimIn(0,1);
GrampaEstimIn = AwgnEstimIn(0,xvar0big); % very mild regularization, for stability

% AWGN measurements 
MeasEstimOut = AwgnEstimOut(y,wvar);

% BG regularization 
AnaEstimOut_ss = AwbgnEstimOut(0,(scaling^2/2),KM*M/N); % spike-and-slab
GrampaEstimOut_ss = EstimOutConcat({MeasEstimOut;AnaEstimOut_ss},[M,N-1]);
end
