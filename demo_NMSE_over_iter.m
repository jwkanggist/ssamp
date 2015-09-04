%------------------------------------------------------------------------%
% filename: demo_NMSE_over_iter.m
% Objective: this script aims to compare the MSE over iteration 
% with algorithms given as 
%             1) ssAMP-BGFD (oracle) (proposed)
%             2) EFLA (Fused Lasso)
%             3) TVAMP-FLSA
%             4) TVAMP-Condat
%             5) Chambolle_Pock TV
%             6) GrAMPA-BG
%             7) ssAMP-BGFD (EM) (proposed)
% Figure plotted by this script will shows convergence behavior of
% algorithms over the number of iterations.
% This script is related to Fig.10 in the below paper.
%
% Jaewook Kang, Hyoyoung Jung, Heung-No Lee, Kiseon Kim,  
% "Bernoulli-Gaussian Approximate Message-Passing Algorithm for Compressed 
% Sensing with 1D-Finite-Difference Sparsity,"  submitted to IEEE TSP at SEP 2015
%
% written by Jaewook Kang 2014 June. finally updated at Sep. 2015
%------------------------------------------------------------------------$
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
path(path, './solver/GAMP_codes');
path(path, './solver/condat1DTV');
path(path, './solver');
path(path, './etc_tool');
%---------------------- Choose solvers --------------------%
usessAMP = true;
useEMssAMP = true;
useEFLA = true;
useTVAMP_FLSA =  true;
useTVAMP_Condat = true;
useCP =  true;
useGrAMPABG = true;
%--------------------- Problem dimension setting -------------------------%
alpha=0.1;% undersampling ratio M/N %0.5 / 0.05 / 0.8
KM=0.1;
N=1000; % the signal vector dimension 
M=round (N*alpha); % the measurement vector dimension
%  Additive noise variance
Delta=1e-10; wvar = max(Delta,1e-10);
disp('%------------------------------------------------------------------------------------------%');
disp('<Simulation codition>')
disp(sprintf('Sparsity K/M = %8.3f',KM));
disp(sprintf('Sampling ratio M/N  = %8.3f',M/N));
disp('%------------------------------------------------------------------------------------------%');
%------------------------ solver configulation ---------------------------
maxiter = 2000;
stop_tol=0;
lambda_CP=0.5;
lambda_EFLA=0.01;
Numofex=2;
%-------------------------------------------------------------------------%
MSE_iter=zeros(Numofex,maxiter,7);

accEFLA=0;accssAMP=0;acc_EMssAMP=0;accTVAMP_FLSA=0;accTVAMP_condat=0;
accCP=0;accGrAMPABG=0;


ttt=0;
Ncnt=Numofex;
for tt=1:Numofex
  
    %---------------------------signal generation-------------------------
    H=randn(M,N)/sqrt(M);
    [Xtrue,normX_sqr]=sig_1DFD_gen(N,KM*alpha,'Gaussian');
    H=mtxgen(M,N,'Gaussian');% measurement matrix
    y=H*Xtrue+sqrt(wvar)*randn(M,1);
   
    %--------------------------- Signal recovery ----------------------------%
    % ssAMP solving (Proposed)
    if usessAMP
        ssAMPout=solve_ssAMP(H,y,0.95,maxiter,stop_tol,0,Xtrue,'Oracle',3,KM*alpha,wvar);
    end

    % ssAMP solving with the EM-tuning (Proposed)
    if useEMssAMP
        ssAMPEMout=solve_ssAMP(H,y,0.95,maxiter,stop_tol,0,Xtrue,'EM');
    end

    % EFLA solving 
    if useEFLA
        opts_EFLA=EFLA_config(maxiter, lambda_EFLA,stop_tol*normX_sqr);
        [estX_EFLA, MSEiter_EFLA, funVal1, ValueL1]= fusedLeastR_MSEiter(H, y,0, opts_EFLA,Xtrue);
    end

    % TVAMP-FLSA solving 
    if useTVAMP_FLSA
        TVAMP_FLSAout=solve_TVAMP_FLSA(H,y,0.5,maxiter,stop_tol,0,Xtrue,'Empirical'); 
    end

    % TVAMP-Condat solving 
    if useTVAMP_Condat
        TVAMP_Condatout=solve_TVAMP_Condat(H,y,0.5,maxiter,stop_tol,0,Xtrue,'Empirical'); 
    end

    % chambolle-pock TV 
    if useCP
        L=max(svd(full(H)));
        TVCPout=solve_chambolle_pock_TV(H,y,H'*y,lambda_CP,maxiter,stop_tol,L,0,Xtrue);
    end

    % GrAMPABG estimate
    if useGrAMPABG
        %Oracle tuning of GrAMPA-BG
        [GrampaEstimIn,GrampaEstimOut_ss,GrampaLinTrans,GrampaOptions]=...
                   config_GrAMPABG(H,y,wvar,KM,stop_tol,maxiter);
         estFin1_GrAMPABG = gampEst_ver4(GrampaEstimIn,GrampaEstimOut_ss,GrampaLinTrans,GrampaOptions,Xtrue);
    end
    %-----------------------------------------------------------------------------------%
    if isnan(sum(ssAMPout.MSE_iter+ssAMPEMout.MSE_iter)) 
        Ncnt=Ncnt-1;
        continue
    else
        ttt=ttt+1;
    end  
  
    MSE_iter(ttt,:,1)= MSEiter_EFLA;
    MSE_iter(ttt,:,2)= ssAMPout.MSE_iter;
    MSE_iter(ttt,:,3)= TVAMP_FLSAout.MSE_iter;
    MSE_iter(ttt,:,4)= TVCPout.MSE_iter;
    MSE_iter(ttt,:,5)= estFin1_GrAMPABG.mse_iter;
    MSE_iter(ttt,:,6)= ssAMPEMout.MSE_iter;
    MSE_iter(ttt,:,7)= TVAMP_Condatout.MSE_iter;
end
MSE_iterdB=10*log10(MSE_iter(1:Ncnt,:,:));
stepsize=180;
 

figure(2)
semilogx (1:maxiter,mean(MSE_iterdB(:,:,1)),'-c',... 
         1:maxiter,mean(MSE_iterdB(:,:,2)),'-r',...
         1:maxiter,mean(MSE_iterdB(:,:,3)),'-k',...
         1:maxiter,mean(MSE_iterdB(:,:,4)),'-m',...
         1:maxiter,mean(MSE_iterdB(:,:,5)),'-g',...
         1:maxiter,mean(MSE_iterdB(:,:,6)),'-b',...
         1:maxiter,mean(MSE_iterdB(:,:,7)),'-y','linewidth',2); hold on;
axis([0 maxiter  -140  20 ]) 
xlabel('Number of iteration, t','fontsize',13)
ylabel('Normalized MSE <dB>','fontsize',12)


title(sprintf(' K/M=%1.2f, M/N=%1.2f',KM,alpha),'fontsize',12)
legend('EFLA','ssAMP-BGFD (Oracle)','TVAMP-FLSA','TV-CP','GrAMPA-BG','ssAMP-BGFD (EM)','TVAMP-Condat')


figure(2);errorbar(2:stepsize:maxiter,mean(MSE_iterdB(:,2:stepsize:maxiter,1)),std(MSE_iterdB(:,2:stepsize:maxiter,1)),'.c'); hold on
figure(2);errorbar(32:stepsize:maxiter,mean(MSE_iterdB(:,32:stepsize:maxiter,2)),std(MSE_iterdB(:,32:stepsize:maxiter,2)),'*r'); hold on
figure(2);errorbar(62:stepsize:maxiter,mean(MSE_iterdB(:,62:stepsize:maxiter,3)),std(MSE_iterdB(:,62:stepsize:maxiter,3)),'sk'); hold on
figure(2);errorbar(92:stepsize:maxiter,mean(MSE_iterdB(:,92:stepsize:maxiter,4)),std(MSE_iterdB(:,92:stepsize:maxiter,4)),'vm'); hold on
figure(2);errorbar(122:stepsize:maxiter,mean(MSE_iterdB(:,122:stepsize:maxiter,5)),std(MSE_iterdB(:,122:stepsize:maxiter,5)),'xg'); hold on
figure(2);errorbar(152:stepsize:maxiter,mean(MSE_iterdB(:,152:stepsize:maxiter,6)),std(MSE_iterdB(:,152:stepsize:maxiter,6)),'ob'); hold on
figure(2);errorbar(182:stepsize:maxiter,mean(MSE_iterdB(:,182:stepsize:maxiter,7)),std(MSE_iterdB(:,182:stepsize:maxiter,7)),'^y'); 




