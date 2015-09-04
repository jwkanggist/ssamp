function opts_EFLA=EFLA_config(maxiter, lambda,tol)
%  EFLA_config.m

opts_EFLA=[];
opts_EFLA.init=2;        % starting from a zero point
% termination criterion
% opts_EFLA.tFlag=5;       % run .maxIter iterations
opts_EFLA.tFlag=3;       % iteration termination when norm( x_i - x_{i-1}, 2) <= .tol
opts_EFLA.maxIter=maxiter;   % maximum number of iterations
opts_EFLA.nFlag=0;       % without normalization
% regularization
opts_EFLA.rFlag=1;       % the input parameter 'rho' is a ratio in (0, 1)
opts_EFLA.fusedPenalty= lambda;
opts_EFLA.lFlag=0; % line search
opts_EFLA.tol=sqrt(tol); 
end
