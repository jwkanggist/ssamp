function TVCPout=...
    solve_chambolle_pock_TV(A,b,x0,lambda,maxiter,tol,L,targetMSE,Xtrue)
%------------------------------------------------------------------------%
% filename :  solve_chambolle_pock_TV.m
%  This solver aim to solver the 1-dimensional TV regularization problem (unconstraint). 
%
%         min_x ||Ax - b ||_2^2 + lambda || x ||_TV
%
%  We refer to the papaer 
%  Sidhy et al, "Convex optimization problem portotyping for image
%  reconstruction in computed tomography with the Chambolle-Pock
%  algorithm", 2012 physics in medicine and biology.
%
%  * ----------------------------------------------------------------------
%  *Input arguments
%  - A (M by N matrix)        : measurement matrix 
%  - b (M by 1 vector)        : measurement vector
%  - x0                       : initial position of x
%  - lambda (scalar )         : The TV regularization parameter
%  - maxiter (scalar)         : the maxinum number of iterations
%  - tol (scalar)         : termination condition (stopping tolerance )
%  - L                         : maximum signgular value 
%  - targetMSE:   if targetMSE = 0, TVAMP will run until "tol" is satisfied
%                 if targetMSE > 0, TVAMP will stop when achieving the targetMSE.
%  - Xtrue  :     Given Xtrue, this solver outputs "MSE_iter" information

%  *output arguments
% - sol_x                     : solution of the optimization
% - end_iternum               : the number of iterations for algorithm termination
% - MSE_iter                  : MSE of the ssAMP estimate over iterations. If Xtrue is not
%                               given, this information is meaningless.
% -
%  * copyright@ Jaewook Kang with Gwangju Institute of Science and
%  * Technology 2014,May,  
%  * feedback: jwkkang@gist.ac.kr
%  * final update Sep, 2015
%-------------------------------------------------------------------------%
[M,N]=size(A);

if nargin < 9
    Xtrue=0;
end

if nargin < 8
    targetMSE=0;
end

if nargin < 6
    tol=1e-7;
end
if nargin < 5
    maxiter=500;
end
if nargin < 4
    lambda=0.01;
end

if nargin < 4
    x0=zeros(N,1);
end


D=sparse(N,N);

curr_x=zeros(N,1);
curr_p=zeros(M,1);
curr_q=zeros(N,1);
curr_x_bar=zeros(N,1);

next_x=zeros(N,1);
next_p=zeros(M,1);
next_q=zeros(N,1);
next_x_bar=zeros(N,1);
lambda_mul_one_vectorN=ones(N,1)*lambda;
% difference matrix (also called TV norm matrix )
I = speye(N-1);
P = spalloc(N-1,1,0);
D(1:N-1,:)=( [P I] - [I P] );
clear I P



% initialization
curr_x_bar=x0;
tau=1/L;
sigma=1/L;
theta=1;

%--for efficient computation ---%
Dt=D';
tauAt=tau*A';
sigmaD=sigma*D;
%---------------------
% iteration start
for n=1:maxiter
   next_p    = (curr_p+sigma*(A*curr_x_bar-b))/(1+sigma);
   temp1=curr_q+sigmaD*curr_x_bar ;
   next_q    = lambda* temp1 ./ max(lambda_mul_one_vectorN, abs( temp1));
   next_x     = curr_x-tauAt*next_p - tau*Dt*next_q;
   next_x_bar = next_x + theta*(next_x-curr_x);
  
  TVCPout.MSE_iter(n)=norm(next_x - Xtrue)^2/norm(Xtrue)^2;

   % stopping criterion
   if targetMSE ==0
       if norm(curr_x-next_x,2)^2/norm(curr_x,2)^2  <tol
           break;
       end
   else
       if norm(next_x-Xtrue,2)^2/norm(Xtrue,2)^2  <targetMSE
           break;
       end
       
   end
   
   curr_x=next_x;
   curr_x_bar=next_x_bar;
   curr_p=next_p;
   curr_q=next_q; 
end

TVCPout.sol_x=next_x;
TVCPout.end_iter=n;

end
