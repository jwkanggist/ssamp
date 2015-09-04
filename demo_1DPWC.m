%------------------------------------------------------------------------%
% Filename: demo_1DPWC.m
% This file is a testbench for an exemplary comparison among 
% algorithms for the 1D piecewise-constant (PWC) recovery problem.  
%
% Each algorithm tries to solve a Compressed sensing recovery problem 
%
% In this comparison, we consider two algorithms, 
% - ssAMP-BGFD (proposed)
% - GrAMPA-BG (http://www2.ece.ohio-state.edu/~schniter/GrAMPA/index.html), 
% which apply an MMSE method:
%
%    x_hat =arg min_x Expectation [ X | Y = y,H] over marginal posterior
%
%  We also include three solvers 
%  - TVAMP (Donoho et al., IEEE TIT 2013)
%  - EFLA (Liu et al., SLEP 4.1 package 2010)
%  - TV-CP (Chambolle and Pock,  Jounral of Math. Imag. 2011)
%, which aim to solve the TV method, defined as
%
%    x_hat = arg min_x || Y- H X ||_2^2 + lambda || DX ||_1
%
% where D is 1D-finite difference matrix.
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
path(path, './solver/GAMP_codes');
path(path, './solver/condat1DTV');
path(path, './solver');
path(path, './etc_tool');
path(path,genpath(pwd));
%---------------------- Choose solvers --------------------%
usessAMP = true;
useEMssAMP = true;
useEFLA = true;
useTVAMP_FLSA =  true;
useTVAMP_Condat = true;
useCP =  true;
useGrAMPABG = true;
%--------------------- Problem dimension setting -------------------------%
alpha=0.5;% undersampling ratio 
KM=0.1;
N=1000; % the signal vector dimension 
M=round (N*alpha); % the measurement vector dimension
%  Additive noise variance
Delta=0; wvar = max(Delta,1e-10);

%------------------ CS measurement generation ----------------------------%

[Xtrue,normX_sqr]=sig_1DFD_gen(N,KM*alpha,'Gaussian');
H=mtxgen(M,N,'Gaussian');% measurement matrix
y=H*Xtrue+sqrt(wvar)*randn(M,1); %Measurement generation

%------------------------ solver configulation ---------------------------
maxiter = 2000;
stop_tol=1e-7;
lambda_CP=0.5;
lambda_EFLA=0.1;

if useEFLA
opts_EFLA=EFLA_config(maxiter, lambda_EFLA,stop_tol*normX_sqr);
end
if useGrAMPABG
    %Oracle tuning of GrAMPA-BG
[GrampaEstimIn,GrampaEstimOut_ss,GrampaLinTrans,GrampaOptions]=...
               config_GrAMPABG(H,y,wvar,KM,stop_tol,maxiter);
end
                                         
intro_demo_1DPWC
%--------------------------- Signal recovery ----------------------------%
% ssAMP solving (Proposed)
if usessAMP
tstart=tic;
ssAMPout=solve_ssAMP(H,y,0.95,maxiter,stop_tol,0,Xtrue,'Oracle',3,KM*alpha,wvar);
telapsed_ssAMP=toc(tstart);
MSE_ssAMP= norm(ssAMPout.mu-Xtrue)^2/normX_sqr;
end

% ssAMP solving with the EM-tuning (Proposed)
if useEMssAMP
tstart=tic;
ssAMPEMout=solve_ssAMP(H,y,0.95,maxiter,stop_tol,0,Xtrue,'EM');
telapsed_EMssAMP=toc(tstart);
MSE_EMssAMP= norm(ssAMPEMout.mu -Xtrue)^2/normX_sqr;
end

% EFLA solving 
if useEFLA
tstart=tic;
[estX_EFLA, funVal1, ValueL1, end_iter_EFLA]= fusedLeastR(H, y, 0, opts_EFLA);
telapsed_fusedlasso=toc(tstart);
MSE_fusedlasso= norm(estX_EFLA-Xtrue)^2/normX_sqr;
end

% TVAMP-FLSA solving 
if useTVAMP_FLSA
tstart=tic;
TVAMP_FLSAout=solve_TVAMP_FLSA(H,y,0.5,maxiter,stop_tol,0,Xtrue,'Empirical'); 
telapsed_TVAMP_FLSA=toc(tstart);
MSE_TVAMP_FLSA= norm(TVAMP_FLSAout.x_hat-Xtrue)^2/normX_sqr;
end

% TVAMP-Condat solving 
if useTVAMP_Condat
tstart=tic;
TVAMP_Condatout=solve_TVAMP_Condat(H,y,0.5,maxiter,stop_tol,0,Xtrue,'Empirical'); 
telapsed_TVAMP_Condat=toc(tstart);
MSE_TVAMP_Condat= norm(TVAMP_Condatout.x_hat-Xtrue)^2/normX_sqr;
end

% chambolle-pock TV 
if useCP
L=max(svd(full(H)));
tstart=tic;
TVCPout=solve_chambolle_pock_TV(H,y,zeros(N,1),lambda_CP,maxiter,stop_tol,L,0,Xtrue);
telapsed_CP=toc(tstart);
MSE_CP=norm(TVCPout.sol_x-Xtrue)^2/normX_sqr;
end



% GrAMPA estimate with bernoulli-Gaussian regularization 
if useGrAMPABG
tstart=tic;
estFin1 = gampEst(GrampaEstimIn,GrampaEstimOut_ss,GrampaLinTrans,GrampaOptions);
telapsed_GrAMPABG=toc(tstart);
estX_GrAMPABG=estFin1.xhat;
MSE_GrAMPABG=norm(estX_GrAMPABG-Xtrue)^2/normX_sqr;
end

%--------------------------- Display ------------------------------------%
disp('%------------------------------------------------------------------------------------------%');
disp('<Recovery result>')
if usessAMP, disp(sprintf('ssAMP-BGFD with oracle tuning: Nomalized MSE  = %8.7f',MSE_ssAMP)); end
if useEMssAMP, disp(sprintf('ssAMP-BGFD with EM tuning: Nomalized MSE  = %8.7f',MSE_EMssAMP)); end
if useEFLA, disp(sprintf('EFLA: Nomalized MSE  = %8.7f',MSE_fusedlasso)); end
if useTVAMP_FLSA, disp(sprintf('TVAMP-FLSA: Nomalized MSE   = %8.7f',MSE_TVAMP_FLSA)); end
if useTVAMP_Condat, disp(sprintf('TVAMP-Condat: Nomalized MSE   = %8.7f',MSE_TVAMP_Condat)); end
if useCP, disp(sprintf('Chambolle-Pock L2TV:Nomalized MSE   = %8.7f',MSE_CP)); end
if useGrAMPABG, disp(sprintf('GrAMPABG: Nomalized MSE   = %8.7f',MSE_GrAMPABG)); end
disp('%------------------------------------------------------------------------------------------%');
if usessAMP, disp(sprintf('Running Time of ssAMP-BGFD with oracle tuning = %8.7f sec',telapsed_ssAMP)); end
if useEMssAMP, disp(sprintf('Running Time of ssAMP-BGFD with EM tuning = %8.7f sec',telapsed_EMssAMP)); end
if useEFLA, disp(sprintf('Running Time of EFLA = %8.7f sec',telapsed_fusedlasso)); end
if useTVAMP_FLSA, disp(sprintf('Running Time of TVAMP-FLSA= %8.7f sec',telapsed_TVAMP_FLSA)); end
if useTVAMP_Condat, disp(sprintf('Running Time of TVAMP-Condat= %8.7f sec',telapsed_TVAMP_Condat)); end
if useCP, disp(sprintf('Running Time of Chambolle-Pock TV = %8.7f sec',telapsed_CP)); end
if useGrAMPABG, disp(sprintf('Running Time of GrAMPA-BG= %8.7f sec',telapsed_GrAMPABG)); end
disp('%------------------------------------------------------------------------------------------%');


figure(1); clf;
if usessAMP,   subplot(4,2,1); plot(1:N,Xtrue,1:N,ssAMPout.mu,'rO');axis([1 N -max(abs(Xtrue))-1 max(abs(Xtrue))+1]); title('(a) ssAMP-BGFD recovery with oracle tuning','fontsize',12);xlabel('Signal index, i','fontsize',12); end
if useEMssAMP, subplot(4,2,2); plot(1:N,Xtrue,1:N,ssAMPEMout.mu,'rO');axis([1 N -max(abs(Xtrue))-1 max(abs(Xtrue))+1]); title('(b) ssAMP-BGFD recovery with EM tuning','fontsize',12);xlabel('Signal index, i','fontsize',12); end
if useEFLA,    subplot(4,2,3); plot(1:N,Xtrue,1:N,estX_EFLA,'rO');axis([1 N -max(abs(Xtrue))-1 max(abs(Xtrue))+1]); title('(c) EFLA recovery','fontsize',12);xlabel('Signal index, i','fontsize',12); end
if useTVAMP_FLSA,   subplot(4,2,4); plot(1:N,Xtrue,1:N,TVAMP_FLSAout.x_hat,'rO');axis([1 N -max(abs(Xtrue))-1 max(abs(Xtrue))+1]); title('(d) TVAMP-FLSA recovery','fontsize',12);xlabel('Signal index, i','fontsize',12); end
if useTVAMP_Condat,  subplot(4,2,5); plot(1:N,Xtrue,1:N,TVAMP_Condatout.x_hat,'rO');axis([1 N -max(abs(Xtrue))-1 max(abs(Xtrue))+1]); title('(e) TVAMP-Condat recovery','fontsize',12);xlabel('Signal index, i','fontsize',12); end
if useCP,      subplot(4,2,6); plot(1:N,Xtrue,1:N,TVCPout.sol_x,'rO');axis([1 N -max(abs(Xtrue))-1 max(abs(Xtrue))+1]); title('(f) Chambolle-Pock recovery with TV norm','fontsize',12);xlabel('Signal index, i','fontsize',13); end
if useGrAMPABG, subplot(4,2,7); plot(1:N,Xtrue,1:N,estX_GrAMPABG,'rO');axis([1 N -max(abs(Xtrue))-1 max(abs(Xtrue))+1]); title('(g) GrAMPA-BG recovery','fontsize',12);xlabel('Signal index, i','fontsize',12); end
box on;
