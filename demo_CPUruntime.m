%-------------------------------------------------------------------------
% filename : demo_CPUruntime.m
% objective : Comparison of running time in he CS recovery with 1D
% finite-difference sparsity. 
% We consider the algorithms as given below:
%             1) ssAMP-BGFD (oracle) (proposed)
%             2) EFLA (Fused Lasso)
%             3) TVAMP-FLSA
%             4) TVAMP-Condat
%             5) Chambolle_Pock TV
%             6) GrAMPA-BG
%             7) ssAMP-BGFD (EM) (proposed)
% Figure plotted by this script will shows convergence behavior of
% algorithms over the number of iterations.
% This script is related to Fig12 in the corresponding paper.

% written by Jaewook Kang 2014 June. finally updated at Sep. 2015
%%-------------------------------------------------------------------------
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
path(path, './solver');
path(path, './etc_tool');
path(path, './solver/GAMP_codes');
path(path, './solver/Condat1DTV');
%---------------------- Choose solvers --------------------%
usessAMP = true;
useEMssAMP = true;
useEFLA = true;
useTVAMP_FLSA =  true;
useTVAMP_Condat = true;
useCP =  true;
useGrAMPABG = true;
%--------------------- Problem dimension setting -------------------------%
alpha=0.5;% undersampling ratio M/N %0.5 / 0.05 / 0.8
KM=0.1;
%  Additive noise variance
Delta=1e-10; wvar = max(Delta,1e-10);
disp('%------------------------------------------------------------------------------------------%');
disp('<Simulation codition>')
disp(sprintf('Sparsity K/M = %8.3f',KM));
disp(sprintf('Undersampling ratio M/N  = %8.3f',alpha));
disp('%------------------------------------------------------------------------------------------%');
%------------------------ solver configulation ---------------------------
maxiter = 2000;
stop_tol=0;
lambda_CP=0.5;
lambda_EFLA=0.01;
TargetMSE=1e-4;
Numofex=10;
N_array=[40*40 60*60 80*80 100*100];
%-------------------------------------------------------------------------%
algoTIME=zeros(7,length(N_array));
for i=1:length(N_array)
    
    acc_time_ssAMP=0;acc_time_EMssAMP=0;acc_time_GrAMPABG=0;acc_time_EFLA=0;
    acc_time_TVAMP_FLSA=0;acc_time_CP=0;acc_time_TVAMP_Condat=0;
    
    acc_MSE_ssAMP=0;acc_MSE_EMssAMP=0;acc_MSE_GrAMPABG=0;acc_MSE_EFLA=0;
    acc_MSE_TVAMP_FLSA=0;acc_MSE_CP=0;acc_MSE_TVAMP_Condat=0;
    
    N=N_array(i); M=round (N*alpha);
    disp(sprintf('Signal length N = %d',N));

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
            tstart=tic;
            ssAMPout=solve_ssAMP(H,y,0.95,maxiter,stop_tol,TargetMSE,Xtrue,'Oracle',3,KM*alpha,wvar);
            telapsed_ssAMP=toc(tstart);
            MSE_ssAMP= norm(ssAMPout.mu-Xtrue)^2/normX_sqr;
        end

        % ssAMP solving with the EM-tuning (Proposed)
        if useEMssAMP
            tstart=tic;
            ssAMPEMout=solve_ssAMP(H,y,0.95,maxiter,stop_tol,TargetMSE,Xtrue,'EM');
            telapsed_EMssAMP=toc(tstart);
            MSE_EMssAMP= norm(ssAMPEMout.mu -Xtrue)^2/normX_sqr;
        end

        % EFLA solving 
        if useEFLA
            opts_EFLA=EFLA_config(maxiter, lambda_EFLA,stop_tol*normX_sqr);
            opts_EFLA.target_MSE=TargetMSE;
            tstart=tic;
            [estX_EFLA, funVal1, ValueL1, end_iter_EFLA]= fusedLeastR(H, y, 0, opts_EFLA);
            telapsed_EFLA=toc(tstart);
            MSE_EFLA= norm(estX_EFLA-Xtrue)^2/normX_sqr;
        end

        % TVAMP-FLSA solving 
        if useTVAMP_FLSA
            tstart=tic;
            TVAMP_FLSAout=solve_TVAMP_FLSA(H,y,0.5,maxiter,stop_tol,TargetMSE,Xtrue,'Empirical'); 
            telapsed_TVAMP_FLSA=toc(tstart);
            MSE_TVAMP_FLSA= norm(TVAMP_FLSAout.x_hat-Xtrue)^2/normX_sqr;
        end

        % TVAMP-Condat solving 
        if useTVAMP_Condat
        tstart=tic;
        TVAMP_Condatout=solve_TVAMP_Condat(H,y,0.5,maxiter,stop_tol,TargetMSE,Xtrue,'Empirical'); 
        telapsed_TVAMP_Condat=toc(tstart);
        MSE_TVAMP_Condat= norm(TVAMP_Condatout.x_hat-Xtrue)^2/normX_sqr;
        end

        % chambolle-pock TV 
        if useCP
            L=max(svd(full(H)));
            tstart=tic;
            TVCPout=solve_chambolle_pock_TV(H,y,zeros(N,1),lambda_CP,maxiter,stop_tol,L,TargetMSE,Xtrue);
            telapsed_CP=toc(tstart);
            MSE_CP=norm(TVCPout.sol_x-Xtrue)^2/normX_sqr;
        end



        % GrAMPA estimate with bernoulli-Gaussian regularization 
        if useGrAMPABG
            [GrampaEstimIn,GrampaEstimOut_ss,GrampaLinTrans,GrampaOptions]=...
            config_GrAMPABG(H,y,wvar,KM,stop_tol,maxiter);
            %*******iteration will be stopped when the target MSE is achieved ********
            GrampaOptions.stopFcn2 = @(GS) (norm(GS.xhat-Xtrue)^2/norm(Xtrue)^2 < TargetMSE); 
            %***************************************************************************
            tstart=tic;
            estFin1 = gampEst(GrampaEstimIn,GrampaEstimOut_ss,GrampaLinTrans,GrampaOptions);
            telapsed_GrAMPABG=toc(tstart);
            estX_GrAMPABG=estFin1.xhat;
            MSE_GrAMPABG=norm(estX_GrAMPABG-Xtrue)^2/normX_sqr;
        end
        
%------------------------- divergence regulation-------------------------%
        if (MSE_CP >TargetMSE) || (MSE_TVAMP_FLSA >TargetMSE) || ...
           (MSE_ssAMP >TargetMSE) || (MSE_GrAMPABG >TargetMSE) ||...
           (MSE_EFLA >TargetMSE) || (MSE_TVAMP_Condat >TargetMSE) 
            Ncnt=Ncnt-1;
            continue;
        end
        
        if isnan(MSE_ssAMP+MSE_EMssAMP)
            Ncnt=Ncnt-1;
            continue;
        end  
        
 %------------------------------------------------------------------------%       
        % averaging runtime
        acc_time_ssAMP=acc_time_ssAMP + telapsed_ssAMP;
        acc_time_EMssAMP=acc_time_EMssAMP + telapsed_EMssAMP; 
        acc_time_GrAMPABG=acc_time_GrAMPABG + telapsed_GrAMPABG;
        acc_time_EFLA=acc_time_EFLA + telapsed_EFLA;
        acc_time_TVAMP_FLSA=acc_time_TVAMP_FLSA +  telapsed_TVAMP_FLSA;
        acc_time_CP=acc_time_CP + telapsed_CP;
        acc_time_TVAMP_Condat=acc_time_TVAMP_Condat +  telapsed_TVAMP_Condat;

        
        % averaging MSE
        acc_MSE_ssAMP = acc_MSE_ssAMP + MSE_ssAMP;
        acc_MSE_EMssAMP = acc_MSE_EMssAMP + MSE_EMssAMP;
        acc_MSE_GrAMPABG = acc_MSE_GrAMPABG + MSE_GrAMPABG;
        acc_MSE_CP = acc_MSE_CP + MSE_CP;
        acc_MSE_EFLA = acc_MSE_EFLA + MSE_EFLA;
        acc_MSE_TVAMP_FLSA = acc_MSE_TVAMP_FLSA + MSE_TVAMP_FLSA;
        acc_MSE_TVAMP_Condat = acc_MSE_TVAMP_Condat + MSE_TVAMP_Condat;
    end
    
        MSE(1,i)=acc_MSE_ssAMP/Ncnt;
        MSE(2,i)=acc_MSE_CP/Ncnt;
        MSE(3,i)=acc_MSE_EFLA/Ncnt;
        MSE(4,i)=acc_MSE_TVAMP_FLSA/Ncnt;
        MSE(5,i)=acc_MSE_GrAMPABG/Ncnt;
        MSE(6,i)=acc_MSE_EMssAMP/Ncnt;
        MSE(7,i)=acc_MSE_TVAMP_Condat/Ncnt;

    
        algoTIME(1,i)=acc_time_ssAMP/Ncnt;
        algoTIME(2,i)=acc_time_CP/Ncnt;
        algoTIME(3,i)=acc_time_EFLA/Ncnt;
        algoTIME(4,i)=acc_time_TVAMP_FLSA/Ncnt;
        algoTIME(5,i)=acc_time_GrAMPABG/Ncnt;
        algoTIME(6,i)=acc_time_EMssAMP/Ncnt;
        algoTIME(7,i)=acc_time_TVAMP_Condat/Ncnt;
        disp('%------------------------------------------------------------------------------------------%');
        disp(sprintf('Counts when all the algorithms achieve TargetMSE: %d',Ncnt));
        disp('%------------------------------------------------------------------------------------------%');
        disp('<Results of average NMSE>')
        disp(sprintf('ssAMP-BGFD (Oracle): Nomalized MSE  = %8.7f',MSE(1,i))); 
        disp(sprintf('ssAMP-BGFD (EM): Nomalized MSE  = %8.7f',MSE(6,i))); 
        disp(sprintf('Chambolle-Pock TV:Nomalized MSE   = %8.7f',MSE(2,i))); 
        disp(sprintf('EFLA: Nomalized MSE  = %8.7f',MSE(3,i))); 
        disp(sprintf('TVAMP-FLSA: Nomalized MSE   = %8.7f',MSE(4,i))); 
        disp(sprintf('TVAMP-Condat: Nomalized MSE   = %8.7f',MSE(7,i)));
        disp(sprintf('GrAMPA-BG: Nomalized MSE   = %8.7f',MSE(5,i))); 
        disp('%------------------------------------------------------------------------------------------%');
       disp('<Average runtime result>')
        disp(sprintf('Telapsed Time of ssAMP-BGFD (Oracle)  = %8.7f sec',algoTIME(1,i)));
        disp(sprintf('Telapsed Time of ssAMP-BGFD (EM)  = %8.7f sec',algoTIME(6,i)));
        disp(sprintf('Telapsed Time of Chambolle_Pock TV = %8.7f sec',algoTIME(2,i)));
        disp(sprintf('Telapsed Time of EFLA = %8.7f sec',            algoTIME(3,i)));
        disp(sprintf('Telapsed Time of TVAMP-FLSA= %8.7f sec',algoTIME(4,i)));
       disp(sprintf('Telapsed Time of TVAMP-Condat= %8.7f sec',algoTIME(7,i)));
        disp(sprintf('Telapsed Time of GrAMPA-BG     = %8.7f sec',algoTIME(5,i)));
        disp('%------------------------------------------------------------------------------------------%');
    
    
end


figure(11)
loglog(N_array,algoTIME(1,:),'X-b',N_array,algoTIME(2,:),'O-.r',N_array,algoTIME(3,:),'d--k',N_array,algoTIME(4,:),'g-s',N_array,algoTIME(5,:),'*:m',N_array,algoTIME(6,:),'c-v',N_array,algoTIME(7,:),'y-^','linewidth',2)
legend('ssAMP-BGFD (Oracle)','TV-CP','EFLA','TVAMP-FLSA','GrAMPA-BG','ssAMP-BGFD (EM-Tuning)','TVAMP-Condat')
xlabel('N','fontsize',14);
ylabel('CPU Runtime <Sec>','fontsize',14);
title(sprintf('K/M=%1.2f, M/N=%1.2f',KM,alpha),'fontsize',14)



