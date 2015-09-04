%------------------------------------------------------------------------%
% Filename: PTcurve_EFLA.m
% This file is a testbench to draw a phase transition curve of the EFLA
% algorithm.
%-------------------------------------------------------------------------
% EFLA : the corresonding paper - J. Liu, L. Yuan, and J. Ye. 
%    An efficient algorithm for a class of fused lasso problems,
%     proc of ACM SIGKDD Conference on Knowledge Discovery and Data Mining, 
%     2010.
%
% This code is used to draw Fig.6 and Fig.7 of the below paper:
%
% Jaewook Kang, Hyoyoung Jung, Heung-No Lee, Kiseon Kim,  
% "Bernoulli-Gaussian Approximate Message-Passing Algorithm for Compressed 
% Sensing with 1D-Finite-Difference Sparsity,"  submitted to IEEE TSP at SEP 2015
%-------------------------------------------------------------------------
% This testbench is only for the use of matlab version of the algorithm
%
% written by Jaewook Kang 2013 Oct., revised SEP. 2015
%------------------------------------------------------------------------%
clc
clear all
close all
% put key subdirectories in path if not already there
path(path, './solver/flsa');
path(path, './solver');
path(path, './etc_tool');
%----------------- parameters for experiment -------------------%
N=625;
exNum=1;
succ_MSE_th=1e-4;% by Schniter's comment
stepsize=0.025;
alpha=0.05:stepsize:0.99;
sparsity=0.05:stepsize:0.99;
Delta=0; wvar = max(Delta,1e-10);
stop_tol=1e-14;
maxiter=2000;% maximum number of iterations
%-----------------------------------------------------------%

lambda_array=10.^(-5:1);

EFLA_succProb=zeros(length(sparsity),length(alpha));
succProb_over_lambda=zeros(length(lambda_array),1);

for j=1:length(alpha)
    M=round(alpha(j)*N);% # of measurements
      
    for i=1:length(sparsity)
        KM=sparsity(i);% signal sparsity
        % fused penalty regularization parameter
        disp('%------------------------------------------------------------------------------------------%');
        disp(sprintf('K/M=%4.4f, M/N=%3.2f',KM,alpha(j)));
        
        for k=1:length(lambda_array)
            lambda_EFLA=lambda_array(k);
            disp(sprintf('Lambda=%1.5f',lambda_EFLA));
            EFLA_succ_cnt=0;
           
            for iter=1:exNum
                %------- measurement generation ------------
                [Xtrue,normX_sqr]=sig_1DFD_gen(N,KM*alpha(j),'Gaussian');
                H=mtxgen(M,N,'Gaussian');% measurement matrix
                y=H*Xtrue+sqrt(wvar)*randn(M,1); %Measurement generation

                % EFLA solving
                opts_EFLA=EFLA_config(maxiter, lambda_EFLA,stop_tol*normX_sqr);
                [estX_EFLA, funVal1, ValueL1]= fusedLeastR(H, y, 0, opts_EFLA);
                MSE_EFLA= norm(estX_EFLA-Xtrue)^2/normX_sqr;
            
                if MSE_EFLA < succ_MSE_th
                    EFLA_succ_cnt=EFLA_succ_cnt +1;
                end
            
            end

            succProb_over_lambda(k)=EFLA_succ_cnt/exNum;
            
        end
        
        EFLA_succProb(i,j)=max(succProb_over_lambda);
        succProb_lasso=EFLA_succProb(i,j);
        
        disp('<Resulting output>')
        disp(sprintf('Fused Lasso(EFLA): Success prob.  = %3.3f',succProb_lasso));  
        
    end
    
end

figure(1)
imagesc(alpha,sparsity,EFLA_succProb);
xlabel('M/N','fontsize',13)
ylabel('K/M','fontsize',13)

