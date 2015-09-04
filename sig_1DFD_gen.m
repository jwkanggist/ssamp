function [X,L2normX_sqr]=sig_1DFD_gen(N,q,sig_attribute)
%-------------------------------------------------------------------------
%
% filename : sig_1DFD_gen.m
%
% This file is for generation of a piecewise-constant signal having its
% signal sparsity in finite-difference.
%
% N             : the signal length
% q             : the sparsity rate K/N
% sig_attribute : the attribute of statistics of the 1d-finite-difference
%                  'Gaussian' --> the nonzero finite-differences are 
%                   Gaussian distributed.
%
%                  'Bernoulli' --> the nonzero finite-differences are
%                  Bernoulli distributed.
%
% written by Jaewook Kang 2015 Aug
% Feedback to jwkkang@gist.ac.kr 
%-------------------------------------------------------------------------
    if nargin ==2
        sig_attribute='Gaussian';
    end
    
    
    
    while(1)
        Numofsupp=length(find(rand(1,N)<q));
        Xstate=randerr2(1,N,Numofsupp);
        supp=find(Xstate==1);
        S=zeros(N,1);S(supp)=1;
        if strcmp(sig_attribute,'Gaussian')==1
            Xdiff = S.*randn(N,1);
        elseif strcmp(sig_attribute,'Bernoulli')==1
            Xdiff = S.*(2*floor(2*rand(N,1))-1);
        else
            Xdiff = S.*randn(N,1);
        end
        X = cumsum(Xdiff);
    
        if sum(diff(X))~=0
            break;
        end
    end
    L2normX_sqr=norm(X)^2;
    
end
