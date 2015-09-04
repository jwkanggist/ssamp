function TVAMPout=solve_TVAMP_FLSA(varargin) %(A,y,beta, maxiter,tol,targetMSE,Xtrue)
%% ======================================================================
%  * solve_TVAMP_FLSA.m: solver of TV-AMP algorithm using FLSA denoiser
%  * 
%  *  This is an AMP algorithm to solver TV minimizatio problem
%  *  We refer to Donoho et al.' paper published in IEEE inform.theory 2013
%  *  For denoiser function in TV-AMP, we have used flsa solver in the SLEP
%  *   package introduced in the paper, 
%  *   [J. Liu, L. Yuan, and J. Ye. “An efficient algorithm for a class
%  *   of fused lasso problems,” proc of ACM SIGKDD Conference on
%  *   Knowledge Discovery and Data Mining, 2010.]
%  *   And the code of the FLSA solver is available at 
%  *             http://www.public.asu.edu/?jye02/Software/SLEP
%  * 
%  * Input arguments
%  * 1) A (M by N): sparse measurement matrix 
%  * 2) y (M by 1): given measurement vector with size M
%  * 3) beta : An iteration damping factor
%  * 4) maxiter: Allowdable maximum number of AMP iterations
%  * 5) tol    : stopping tolerance - TVAMP will be terminated when the nomarlized MSE <= tol 
%  * 6) targetMSE: if targetMSE = 0, TVAMP will run until "tol" is satisfied
%                 if targetMSE > 0, TVAMP will stop when achieving the targetMSE.
%  * 7) Xtrue  :  Given Xtrue, this solver outputs "MSE_iter" information
%  * 8) mode   : the 'lambda' paramter tuning mode if mode=='Empirical', this solver will finds GrAMPABGly optimal lambda at every iteration 
%  *                                               if mode=='Normal', the lambda is given by the input argument. 
%  * 9) lambda  : control parameter for the TV regularizer
%  * 
%  * Output arguments
% * 1) x_hat: signal estimate from the AMP iteration 
% * 2) end_iternum: the number of iterations for algorithm termination
% * 3) MSE_iter: MSE of the ssAMP estimate over iterations. If Xtrue is not
%                given, this information is meaningless.
%  *
%  * copyright@ Jaewook Kang with Gwangju Institute of Science and
%  * Technology 2014,Jan,  
%  * feedback: jwkkang@gist.ac.kr
%  * final update Sep, 2015
%%  *=======================================================================*/
A=varargin{1};y=varargin{2};

if nargin <8
    mode='Empirical';
else
    mode=varargin{8};   
    if strcmp(mode,'Normal')==1
        if nargin <9
            lambda=5;
        else
            lambda=varargin{9};
        end
    end   
end

if nargin < 7
    Xtrue=0;
else
    Xtrue=varargin{7};
end

if nargin < 6
    targetMSE=0;
else
    targetMSE=varargin{6};
end

if nargin < 5
    tol=1e-7;
else
    tol=varargin{5};
end

if nargin < 4
    maxiter=500;
else
    maxiter=varargin{4};
end

if nargin < 3
    beta= 0.5;
else
    beta=varargin{3};
end


[M,N]=size(A);
x_hat=zeros(N,1);
z=zeros(M,1);

prev_x=zeros(N,1);
prev_z=zeros(M,1);

%initialization
prev_z=y;

% parameter for tvdip solver
maxiterTV=20;
lratio =10.^(-6:1);
b0=zeros(N-1,1);

for t=2:maxiter+1
    rho=A'*prev_z+prev_x;
     prev_x=x_hat;
     xMSE=1;
     
     if strcmp(mode,'Empirical')==1
         % Empirical tuning for lambda_TV at every iteration
         for tt=1:length(lratio )
            lambda=lratio(tt);
            [temp_x, b]=flsa(rho, b0,0, lambda, N,maxiterTV, 1e-5, 1, 6);
            temp_xMSE=norm(Xtrue-temp_x)^2/norm(Xtrue)^2;

            if xMSE > temp_xMSE
                xMSE=temp_xMSE;
                b0=b;
                x_hat=temp_x;
            end
         end
     else
         [x_hat, b]=flsa(rho, b0,0, lambda, N,maxiterTV, 1e-5, 1, 6);
     end
     
     
   TVAMPout.MSE_iter(t-1)=norm(x_hat-Xtrue)^2/norm(Xtrue)^2;
     
   % loop termination condition check 
   if targetMSE == 0
        if (norm(prev_x - x_hat)^2/norm(prev_x)^2 <tol) 
            break;
        end
   else
        if (norm(x_hat-Xtrue)^2/norm(Xtrue)^2 <targetMSE) 
            break;
        end
   end
   
    %--------------------------------------------------------------------%
    % Onsager term calculation
    if M/N  >=0.5
        Onsager = 1/M* (length(nonzeros(diff(prev_x))));  
    else
        Onsager = 1/M* (length(nonzeros(diff(prev_x))))-1;  % for small M/N
    end
    
     % residual message update
     z = (1-beta)*prev_z +  beta*(y - A*x_hat +  prev_z * Onsager) ;
     prev_x=x_hat;
     prev_z=z;
     
end
TVAMPout.end_iternum=t-1;
TVAMPout.x_hat=x_hat;
end


    
