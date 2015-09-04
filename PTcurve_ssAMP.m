%------------------------------------------------------------------------%
% Filename: PTcurve_ssAMP.m
% This file is a testbench to draw a phase transition curve of the ssAMP-BGFD
% algorithm.
%-------------------------------------------------------------------------
% Jaewook Kang, Hyoyoung Jung, Heung-No Lee, Kiseon Kim,  
% "Bernoulli-Gaussian Approximate Message-Passing Algorithm for Compressed 
% Sensing with 1D-Finite-Difference Sparsity,"  submitted 2015
%
% This code is used to draw Fig.6, Fig.7 and Fig.8 of the above paper. 
% -------------------------------------------------------------------------
%
% written by Jaewook Kang 2014 May.  updated at Sep. 2015
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
path(path, './solver/denoiser');
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
ssAMP_succProb=zeros(length(sparsity),length(alpha));

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
            
            % ssAMP solving
            ssAMPEMout=solve_ssAMP(H,y,0.95,maxiter,stop_tol,0,Xtrue,'EM');
            MSE_EMssAMP= norm(ssAMPEMout.mu -Xtrue)^2/normX_sqr;

            if MSE_EMssAMP < succ_MSE_th
                AMP_succ_cnt=AMP_succ_cnt +1;
            end
            
        end

        ssAMP_succProb(i,j)=AMP_succ_cnt/exNum;
        succProb_AMP=ssAMP_succProb(i,j);
        
        disp('<Resulting output>')
        disp(sprintf('ssAMP: Success prob.  = %3.3f',succProb_AMP));  
        
    end
end

figure(2)
imagesc(alpha,sparsity,ssAMP_succProb);
xlabel('M/N','fontsize',16)
ylabel('K/M','fontsize',16)


