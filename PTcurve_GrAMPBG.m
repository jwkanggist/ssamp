%------------------------------------------------------------------------%
% Filename: PTcurve_GrAMPBG.m
% This file is a testbench to draw a phase transition curve of the ssAMP-BGFD
% algorithm.
%------------------------------------------------------------------------
%  GrAMPA-BG: the corresonding paper - M. Borgerding and P. Schiniter, 
%     “Generalized approximate message passing for the cosparse analysis
%      model,” ICASSP 2015.
%
% This code is used to draw Fig.6 and Fig.7 of the below paper:
%
% Jaewook Kang, Hyoyoung Jung, Heung-No Lee, Kiseon Kim,  
% "Bernoulli-Gaussian Approximate Message-Passing Algorithm for Compressed 
% Sensing with 1D-Finite-Difference Sparsity,"  submitted to IEEE TSP at SEP 2015
%--------------------------------------------------------------------------
% This testbench is only for the use of matlab version of the algorithm
%
% written by Jaewook Kang 2013 Oct.  updated at Sep. 2015
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
path(path, './solver/GAMP_codes');
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

GrAMPABG_succProb=zeros(length(sparsity),length(alpha));

for j=1:length(alpha)
        M=round(alpha(j)*N);% # of measurements
        %---------------------------------------------------------------%
    for i=1:length(sparsity)
        KM=sparsity(i);% signal sparsity

        disp('%------------------------------------------------------------------------------------------%');
        disp(sprintf('K/M=%4.4f, M/N=%3.2f',KM,alpha(j)));
        GrAMPABG_succ_cnt=0;
        
        for iter=1:exNum
            %------- measurement generation ------------
            [Xtrue,normX_sqr]=sig_1DFD_gen(N,KM*alpha(j),'Gaussian');
            H=mtxgen(M,N,'Gaussian');% measurement matrix
            y=H*Xtrue+sqrt(wvar)*randn(M,1); %Measurement generation


            % GrAMPA estimate with bernoulli-Gaussian regularization 
            %Oracle tuning of GrAMPA-BG
            [GrampaEstimIn,GrampaEstimOut_ss,GrampaLinTrans,GrampaOptions]=...
               config_GrAMPABG(H,y,wvar,KM,stop_tol,maxiter);
            estFin1 = gampEst(GrampaEstimIn,GrampaEstimOut_ss,GrampaLinTrans,GrampaOptions);
            estX_GrAMPABG=estFin1.xhat;
            MSE_GrAMPABG=norm(estX_GrAMPABG-Xtrue)^2/normX_sqr;

            if MSE_GrAMPABG < succ_MSE_th
                GrAMPABG_succ_cnt=GrAMPABG_succ_cnt +1;
            end
            
        end

        GrAMPABG_succProb(i,j)=GrAMPABG_succ_cnt/exNum;
        succProb_GrAMPABG=GrAMPABG_succProb(i,j);
        
        disp('<Resulting output>')
        disp(sprintf('GrAMPABG: Success prob.  = %3.3f',succProb_GrAMPABG));  
        
    end


end

figure(2)
imagesc(alpha,sparsity,GrAMPABG_succProb);
xlabel('M/N','fontsize',13)
ylabel('K/M','fontsize',13)



