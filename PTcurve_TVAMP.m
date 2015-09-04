%------------------------------------------------------------------------%
% Filename: PTcurve_TVAMP.m
% This file is a testbench to draw a phase transition curve of the TVAMP
% algorithm.
%-------------------------------------------------------------------------
% TV-AMP: the corresonding paper - D. L. Donoho, I. Johnstone, and A. Montanari, 
%     “Accurate prediction of phase transitions in compressed sensing via 
%      a connection to minimax denoising, ” IEEE Trans. Inform. Theory, 
%      vol. 59, no. 6, pp. 3396-3433, June 2013.
%
% This code is used to draw Fig.6 and Fig.7 of the below paper:
%
% Jaewook Kang, Hyoyoung Jung, Heung-No Lee, Kiseon Kim,  
% "Bernoulli-Gaussian Approximate Message-Passing Algorithm for Compressed 
% Sensing with 1D-Finite-Difference Sparsity,"  submitted to IEEE TSP at SEP 2015
%-------------------------------------------------------------------------
%
% written by Jaewook Kang 2014.Feb, the last updated at Sep. 2015
%------------------------------------------------------------------------%
clc
clear all
close all


%Handle random seed
if verLessThan('matlab','7.14')
  defaultStream = RandStream.getDefaultStream;
else
  defaultStream = RandStream.getGlobalStream;
end;

if 1
    savedState = defaultStream.State;
    save random_state.mat savedState;
else
    load random_state.mat
end
defaultStream.State = savedState;

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
TVAMP_succProb=zeros(length(sparsity),length(alpha));

for j=1:length(alpha)
        M=round(alpha(j)*N);% # of measurements
        %---------------------------------------------------------------%
    for i=1:length(sparsity)
        KM=sparsity(i);% signal sparsity

        disp('%------------------------------------------------------------------------------------------%');
        disp(sprintf('K/M=%4.4f, M/N=%3.2f',KM,alpha(j)));
        AMP_succ_cnt=0;
        
        for iter=1:exNum
            %------- measurement generation ------------
            [Xtrue,normX_sqr]=sig_1DFD_gen(N,KM*alpha(j),'Gaussian');
            H=mtxgen(M,N,'Gaussian');% measurement matrix
            y=H*Xtrue+sqrt(wvar)*randn(M,1); %Measurement generation            


            % TV-AMP solving
            TVAMP_FLSAout=solve_TVAMP_FLSA(H,y,0.5,maxiter,stop_tol,0,Xtrue,'Empirical'); 
            MSE_TVAMP_FLSA= norm(TVAMP_FLSAout.x_hat-Xtrue)^2/normX_sqr;
  
            if  MSE_TVAMP_FLSA < succ_MSE_th
                AMP_succ_cnt=AMP_succ_cnt +1;
            end
            
        end

        TVAMP_succProb(i,j)=AMP_succ_cnt/exNum;
        succProb_TVAMP=TVAMP_succProb(i,j);
        
        disp('<Resulting output>')
        disp(sprintf('TVAMP: Success prob.  = %3.3f',succProb_TVAMP));  
        
    end
    
end

figure(22)
imagesc(alpha,sparsity,TVAMP_succProb);
xlabel('M/N','fontsize',13)
ylabel('K/M','fontsize',13)

