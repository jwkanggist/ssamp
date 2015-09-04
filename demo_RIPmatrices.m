%------------------------------------------------------------------------%
% Filename: demo_RIPmatrices.m
% This file is a testbench for ssAMP-BGFD recovery to
% the compressed sensing recovery with 1D-finite-difference sparsity for a
% variety of the measurement matrice satisfying RIP matrices:
%
%  - Standard Gaussian 
%  - Standard Bernoulli 
%  - Discrete-Cosine Transform with column sign randomization
%  - Walsh-Hadamard Transform with column sign randomization
%  - Quasi-Toeplitz with column sign randomization
%  - Amini's deterministic bipolar matrix with column sign randomization
%  (This matrix is restricted to a case of RIP order 8 and N=1024) 
%  ( We borrow Amini's Matlab code to generate this matrix)
%
%  - Sparse-Bernoulli
%
%
% In this comparison, we consider our proposed algorithms, 
%  ssAMP-BGFD (proposed), which apply an MMSE method:
%
%    x_hat =arg min_x Expectation [ X | Y = y,H] over marginal posterior
%
%   , which is given in the paper
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
%---------------------- Choose solvers --------------------%
use_Gaussian= true;
use_Bernoulli=true;
use_DCT =true;
use_WHT=true;
use_Toep=true;
use_Amini=true;
use_sparse=true;
%--------------------- Problem dimension setting -------------------------%
alpha=0.25;% undersampling ratio M/N %0.5 / 0.05 / 0.8
KM=0.1;
N=1024; % the signal vector dimension 
M=round (N*alpha); % the measurement vector dimension
%  Additive noise variance
Delta=0; wvar = max(Delta,1e-10);

%------------------ CS measurement generation ----------------------------%

[Xtrue,normX_sqr]=sig_1DFD_gen(N,KM*alpha,'Gaussian');

if use_Gaussian
    H_Gaussian=mtxgen(M,N,'Gaussian');
    y_Gaussian=H_Gaussian*Xtrue+sqrt(wvar)*randn(M,1); 
end

if use_Bernoulli
    H_Bernoulli=mtxgen(M,N,'Bernoulli');
    y_Bernoulli=H_Bernoulli*Xtrue+sqrt(wvar)*randn(M,1); 
end

if use_DCT
    H_DCT=mtxgen(M,N,'DCT');
    y_DCT=H_DCT*Xtrue+sqrt(wvar)*randn(M,1); 
end

if use_WHT
    H_WHT=mtxgen(M,N,'WHT');
    y_WHT=H_WHT*Xtrue+sqrt(wvar)*randn(M,1); 
end

if use_Toep
    H_Toep=mtxgen(M,N,'Toeplitz',0.75);
    y_Toep=H_Toep*Xtrue+sqrt(wvar)*randn(M,1); 
end

if use_Amini
    H_Amini=BCH_CS_Matrix(N ,8)*diag(2*floor(2*rand(N,1))-1);% Amini's code
    y_Amini=H_Amini*Xtrue+sqrt(wvar)*randn(M-1,1); 
end

if use_sparse
    H_sparse=mtxgen(M,N,'sparse',0.03);
    y_sparse=H_sparse*Xtrue+sqrt(wvar)*randn(M,1); 
end
%------------------------ solver configulation ---------------------------
maxiter = 2000;
stop_tol=1e-7;
%--------------------------- Signal recovery ----------------------------%


% ssAMP solving with the EM-tuning (Proposed)
if use_Gaussian
tstart=tic;
ssAMPEMout_Gaussian=solve_ssAMP(H_Gaussian,y_Gaussian,0.95,maxiter,stop_tol,0,Xtrue,'EM');
telapsed_EMssAMP_Gaussian=toc(tstart);
MSE_EMssAMP_Gaussian= norm(ssAMPEMout_Gaussian.mu -Xtrue)^2/normX_sqr;
end

if use_Bernoulli
tstart=tic;
ssAMPEMout_Bernoulli=solve_ssAMP(H_Bernoulli,y_Bernoulli,0.95,maxiter,stop_tol,0,Xtrue,'EM');
telapsed_EMssAMP_Bernoulli=toc(tstart);
MSE_EMssAMP_Bernoulli= norm(ssAMPEMout_Bernoulli.mu -Xtrue)^2/normX_sqr;
end

if use_DCT
tstart=tic;
ssAMPEMout_DCT=solve_ssAMP(H_DCT,y_DCT,0.95,maxiter,stop_tol,0,Xtrue,'EM');
telapsed_EMssAMP_DCT=toc(tstart);
MSE_EMssAMP_DCT= norm(ssAMPEMout_DCT.mu -Xtrue)^2/normX_sqr;
end

if use_WHT
tstart=tic;
ssAMPEMout_WHT=solve_ssAMP(H_WHT,y_WHT,0.95,maxiter,stop_tol,0,Xtrue,'EM');
telapsed_EMssAMP_WHT=toc(tstart);
MSE_EMssAMP_WHT= norm(ssAMPEMout_WHT.mu -Xtrue)^2/normX_sqr;
end

if use_Toep
tstart=tic;
ssAMPEMout_Toep=solve_ssAMP(H_Toep,y_Toep,0.5,maxiter,stop_tol,0,Xtrue,'EM');
telapsed_EMssAMP_Toep=toc(tstart);
MSE_EMssAMP_Toep= norm(ssAMPEMout_Toep.mu -Xtrue)^2/normX_sqr;
end

if use_Amini
tstart=tic;
ssAMPEMout_Amini=solve_ssAMP(H_Amini,y_Amini,0.95,maxiter,stop_tol,0,Xtrue,'EM');
telapsed_EMssAMP_Amini=toc(tstart);
MSE_EMssAMP_Amini= norm(ssAMPEMout_Amini.mu -Xtrue)^2/normX_sqr;
end

if use_sparse
tstart=tic;
ssAMPEMout_sparse=solve_ssAMP(H_sparse,y_sparse,0.95,maxiter,stop_tol,0,Xtrue,'EM');
telapsed_EMssAMP_sparse=toc(tstart);
MSE_EMssAMP_sparse= norm(ssAMPEMout_sparse.mu -Xtrue)^2/normX_sqr;
end

%--------------------------- Display ------------------------------------%
disp('%------------------------------------------------------------------------------------------%');
disp('<NMSE results>')
if use_Gaussian, disp(sprintf('With the standard Gaussian matrix: NMSE  = %8.7f',MSE_EMssAMP_Gaussian)); end
if use_Bernoulli, disp(sprintf('With the standard Bernoulli matrix: NMSE  = %8.7f',MSE_EMssAMP_Bernoulli)); end
if use_DCT, disp(sprintf('With the DCT matrix: NMSE  = %8.7f',MSE_EMssAMP_DCT)); end
if use_WHT, disp(sprintf('With the WHT matrix:  NMSE   = %8.7f',MSE_EMssAMP_WHT)); end
if use_Toep, disp(sprintf('With the Quasi-Toeplitz matrix:  NMSE   = %8.7f',MSE_EMssAMP_Toep)); end
if use_Amini, disp(sprintf('With the Amini''s deterministic bipolar matrix:  NMSE   = %8.7f',MSE_EMssAMP_Amini)); end
if use_sparse, disp(sprintf('With the Bernoulli-sparse matrix:  NMSE   = %8.7f',MSE_EMssAMP_sparse)); end
disp('%------------------------------------------------------------------------------------------%');
disp('<CPU runtime result>')
if use_Gaussian, disp(sprintf('With the standard Gaussian matrix: %8.7f sec',telapsed_EMssAMP_Gaussian)); end
if use_Bernoulli, disp(sprintf('With the standard Bernoulli matrix: %8.7f sec',telapsed_EMssAMP_Bernoulli)); end
if use_DCT, disp(sprintf('With the DCT matrix: %8.7f sec',telapsed_EMssAMP_DCT)); end
if use_WHT, disp(sprintf('With the WHT matrix: %8.7f sec',telapsed_EMssAMP_WHT)); end
if use_Toep, disp(sprintf('With the Quasi-Toeplitz matrix:  %8.7f sec',telapsed_EMssAMP_Toep)); end
if use_Amini, disp(sprintf('With the Amini''s deterministic bipolar matrix:  %8.7f sec',telapsed_EMssAMP_Amini)); end
if use_sparse, disp(sprintf('With the Bernoulli-sparse matrix: %8.7f sec',telapsed_EMssAMP_sparse)); end
disp('%------------------------------------------------------------------------------------------%');


figure(1); clf; colormap('Hot');
if use_Gaussian,  subplot(7,2,1); plot(1:N,Xtrue,1:N,ssAMPEMout_Gaussian.mu,'rO');axis([1 N -max(abs(Xtrue))-1 max(abs(Xtrue))+1]); title('(a) With Standard Gaussian Matrix','fontsize',12);xlabel('Signal index, i','fontsize',12); end
if use_Gaussian,  subplot(7,2,2); imagesc(H_Gaussian);title('Standard Gaussian Matrix','fontsize',12);end
if use_Bernoulli, subplot(7,2,3); plot(1:N,Xtrue,1:N,ssAMPEMout_Bernoulli.mu,'rO');axis([1 N -max(abs(Xtrue))-1 max(abs(Xtrue))+1]); title('(b) With Standard Bernoulli Matrix','fontsize',12);xlabel('Signal index, i','fontsize',12); end
if use_Bernoulli, subplot(7,2,4); imagesc(H_Bernoulli);title('Standard Bernoulli Matrix','fontsize',12);end
if use_DCT,       subplot(7,2,5); plot(1:N,Xtrue,1:N,ssAMPEMout_DCT.mu,'rO');axis([1 N -max(abs(Xtrue))-1 max(abs(Xtrue))+1]); title('(c) With DCT Matrix + column sign randomization','fontsize',12);xlabel('Signal index, i','fontsize',12); end
if use_DCT,       subplot(7,2,6);imagesc(H_DCT);title('DCT  Matrix','fontsize',12);end
if use_WHT,       subplot(7,2,7); plot(1:N,Xtrue,1:N,ssAMPEMout_WHT.mu,'rO');axis([1 N -max(abs(Xtrue))-1 max(abs(Xtrue))+1]); title('(d)  With WHT Matrix + column sign randomization','fontsize',12);xlabel('Signal index, i','fontsize',12); end
if use_WHT,       subplot(7,2,8); imagesc(H_WHT);title('WHT Matrix','fontsize',12);end
if use_Toep,      subplot(7,2,9); plot(1:N,Xtrue,1:N,ssAMPEMout_Toep.mu,'rO');axis([1 N -max(abs(Xtrue))-1 max(abs(Xtrue))+1]); title('(e) With Quasi Toeplitz + column sign randomization','fontsize',12);xlabel('Signal index, i','fontsize',12); end
if use_Toep,      subplot(7,2,10); imagesc(H_Toep);title('Quasi-Toeplitz Matrix','fontsize',12);end
if use_Amini,     subplot(7,2,11); plot(1:N,Xtrue,1:N,ssAMPEMout_Amini.mu,'rO');axis([1 N -max(abs(Xtrue))-1 max(abs(Xtrue))+1]); title('(f) With Amini''s deterministic bipolar+ column sign randomization','fontsize',12);xlabel('Signal index, i','fontsize',12); end
if use_Amini,     subplot(7,2,12); imagesc(H_Amini);title('Amini deterministic bipolar Matrix','fontsize',12);end
if use_sparse,    subplot(7,2,13); plot(1:N,Xtrue,1:N,ssAMPEMout_sparse.mu,'rO');axis([1 N -max(abs(Xtrue))-1 max(abs(Xtrue))+1]); title('(g) With Bernoulli-sparse Matrix','fontsize',12);xlabel('Signal index, i','fontsize',12); end
if use_sparse,    subplot(7,2,14); imagesc(H_sparse);title('Bernoulli-sparse Matrix','fontsize',12);end
box on;




