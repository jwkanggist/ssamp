%------------------------------------------------------------------------%
% Filename: PTcurve_TVCP.m
% This file is a testbench to draw a phase transition diagramcruve
%------------------------------------------------------------------------
% Chambolle-Pock: the corresonding paper - A. Chambolle, T. Pock, 
%    “A first-order primal-dual algorithm for convex
%      problems with applications to imaging,” J. Math. Imag. Vis., vol. 40, pp.
%      120-145, May 2011.
%
% This code is used to draw Fig.6 and Fig.7 of the below paper:
%
% Jaewook Kang, Hyoyoung Jung, Heung-No Lee, Kiseon Kim,  
% "Bernoulli-Gaussian Approximate Message-Passing Algorithm for Compressed 
% Sensing with 1D-Finite-Difference Sparsity,"  submitted to IEEE TSP at SEP 2015
% -------------------------------------------------------------------------
% This testbench is only for the use of matlab version of the algorithm
%
% written by Jaewook Kang 2013 Oct., revised SEP. 2015
%------------------------------------------------------------------------%
clc
clear all
close all
% put key subdirectories in path if not already there
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
lambda_array= [0.01  0.05  0.1  0.5 1 5 10];


CP_succProb=zeros(length(sparsity),length(alpha));
succProb_over_lambda=zeros(length(lambda_array),1);


for j=1:length(alpha)
    M=round(alpha(j)*N);% # of measurements
      
    for i=1:length(sparsity)
        KM=sparsity(i);% signal sparsity
        % fused penalty regularization parameter
        disp('%------------------------------------------------------------------------------------------%');
        disp(sprintf('K/M=%4.4f, M/N=%3.2f',KM,alpha(j)));
        
        for k=1:length(lambda_array)
         
            lambda=lambda_array(k);
            disp(sprintf('TV Penalty=%1.5f',lambda));
            CP_succ_cnt=0;
           
            for iter=1:exNum
                %------- measurement generation ------------
                [Xtrue,normX_sqr]=sig_1DFD_gen(N,KM*alpha(j),'Gaussian');
                H=mtxgen(M,N,'Gaussian');% measurement matrix
                y=H*Xtrue+sqrt(wvar)*randn(M,1); %Measurement generation 
                
                % TV-CP solving
                L=max(svd(full(H)));
                TVCPout=solve_chambolle_pock_TV(H,y,zeros(N,1),lambda,maxiter,stop_tol,L,0,Xtrue);
                MSE_CP=norm(TVCPout.sol_x-Xtrue)^2/normX_sqr;
          
            
                if MSE_CP < succ_MSE_th
                    CP_succ_cnt=CP_succ_cnt +1;
                end
            
            end

            succProb_over_lambda(k)=CP_succ_cnt/exNum;
            
        end
        
        CP_succProb(i,j)=max(succProb_over_lambda);
        succProb_CP=CP_succProb(i,j);
        
        disp('<Resulting output>')
        disp(sprintf('TV-CP: Success prob.  = %3.3f',succProb_CP));  
        
    end
    
end

figure(1)
imagesc(alpha,sparsity,CP_succProb);
xlabel('M/N','fontsize',13)
ylabel('K/M','fontsize',13)


 
