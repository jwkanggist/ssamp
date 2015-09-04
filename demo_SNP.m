%------------------------------------------------------------------------%
% Filename: demo_SNP.m
% This file is a testbench for  the CS recovery with SNP genomic data.
%
% In this testbench, we totally compare the four algorithms below:
%
%  - ssAMP-BGFD (proposed)
%  - GrAMPA-BG (http://www2.ece.ohio-state.edu/~schniter/GrAMPA/index.html), 
%  - TVAMP-Condat (Donoho et al., IEEE TIT 2013, Condat, IEEE SPL 2013)
%  - TV-CP (Chambolle and Pock,  Jounral of Math. Imag. 2011)
%
% This file is related to Fig13 of the paper 
%
%   Jaewook Kang, Hyoyoung Jung, Heung-No Lee, Kiseon Kim,  
%   "Bernoulli-Gaussian Approximate Message-Passing Algorithm for Compressed 
%   Sensing with 1D-Finite-Difference Sparsity,"  submitted 2015
%
% written by Jaewook Kang 2014 Aug. updated at Sep. 2015
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
path(path, './solver/denoiser');
path(path, './solver/condat1DTV');
path(path, './solver/GAMP_codes');
path(path, './solver');
path(path, './etc_tool');
path(path,genpath(pwd));


%-------------- Genomic sample selection ------------------------------%
load GBMpart_short_v2.mat
SNPindex=41;
X=GBMpart_short_v2(:,SNPindex);
normX_sqr=norm(X)^2;
%---------------------- Choose solvers --------------------%
usessAMP = true;
useTVAMP_Condat = true;
useCP =  true;
useGrAMPABG = true;
%--------------------- Problem dimension setting -------------------------%
alpha=0.5;% undersampling ratio M/N %0.5 / 0.05 / 0.8
N=length(X); % the signal vector dimension 
M=round (N*alpha); % the measurement vector dimension

%------------------ CS measurement generation ----------------------------%
H=mtxgen(M,N,'Gaussian');% measurement matrix
y=H*X; %Measurement generation

%------------------------ solver configulation ---------------------------
Delta=0.22; q=1e-6;
Delta_grampa=0.25; q_grampa=1e-5;
maxiter = 2000;
SNR = (mean(abs(y).^2)-Delta)/Delta; 
stop_tol = (max(1e-7,min(1e-4,(1/SNR))))^2; 

lambda_CP=10;
lambda_TVAMP=10;

if useGrAMPABG
        %Oracle tuning of GrAMPA-BG
[GrampaEstimIn,GrampaEstimOut_ss,GrampaLinTrans,GrampaOptions]=...
               config_GrAMPABG(H,y,Delta_grampa,q_grampa/alpha,stop_tol,maxiter);
end
                                         
intro_demo_1DPWC_to_SNPdata
%--------------------------- Signal recovery ----------------------------%
% ssAMP-BGFD  (Proposed)
if usessAMP
tstart=tic;
ssAMPout=solve_ssAMP(H,y,0.95,maxiter,stop_tol,0,0,'Oracle',1,q,Delta);
telapsed_ssAMP=toc(tstart);
end

% TVAMP-Condat solving 
if useTVAMP_Condat
tstart=tic;
TVAMP_Condatout=solve_TVAMP_Condat(H,y,0.5,maxiter,stop_tol,0,0,'Normal',lambda_TVAMP); 
telapsed_TVAMP_Condat=toc(tstart);
end

% chambolle-pock TV 
if useCP
L=max(svd(full(H)));
tstart=tic;
TVCPout=solve_chambolle_pock_TV(H,y,zeros(N,1),lambda_CP,maxiter,stop_tol,L);
telapsed_CP=toc(tstart);
end

% GrAMPA estimate with spike-and-slab regularization 
if useGrAMPABG
tstart=tic;
estFin1 = gampEst(GrampaEstimIn,GrampaEstimOut_ss,GrampaLinTrans,GrampaOptions);
telapsed_GrAMPABG=toc(tstart);
estX_GrAMPABG=estFin1.xhat;
end

%--------------------------- Display ------------------------------------%

figure(1); clf;
if usessAMP,   subplot(1,4,1); plot(1:N,X,'c.',1:N,ssAMPout.mu,'r-','linewidth',2);axis([1 N -2 4]); title(sprintf('CPU runtime:%2.4f sec \n # of iterations:%d',telapsed_ssAMP,ssAMPout.end_iternum),'fontsize',10);  end
if useTVAMP_Condat,   subplot(1,4,2); plot(1:N,X,'c.',1:N,TVAMP_Condatout.x_hat,'r-','linewidth',2);axis([1 N -2 4]);title(sprintf('CPU runtime:%2.4f sec\n # of iterations:%d',telapsed_TVAMP_Condat,TVAMP_Condatout.end_iternum),'fontsize',10);  end
if useCP,      subplot(1,4,3); plot(1:N,X,'c.',1:N,TVCPout.sol_x,'r-','linewidth',2);axis([1 N -2 4]);title(sprintf('CPU runtime:%2.4f sec\n # of iterations:%d',telapsed_CP,TVCPout.end_iter),'fontsize',10);  end
if useGrAMPABG, subplot(1,4,4); plot(1:N,X,'c.',1:N,estX_GrAMPABG,'r-','linewidth',2);axis([1 N -2 4]);title(sprintf('CPU runtime:%2.4f sec\n # of iterations:%d',telapsed_GrAMPABG,estFin1.nit),'fontsize',10);  end
box on;

