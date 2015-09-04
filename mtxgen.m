function H=mtxgen(varargin)
%-------------------------------------------------------------------------

% filename: mtxgen.m

% This is a file for generating the compressed sensing matrix H according
% to the input parameters below:

% M             : the number of the rows
% N             : the number of the columns
% mtx_attribute : the matrix attribute (Default is 'Gaussian')
%                  'Gaussian'  -> the standard i.i.d. Gaussian matrix H
%                  'Bernoulli' -> the i.i.d. Bernoulli matrix H
%                  'DCT'       -> the discrete-cosine-transform matrix H
%                                 with column-sign-randomization
%                  'WHT'       -> the Walsh-Hadamard matrix H with
%                                 column-sign-randomization
%                  'Toeplitz'  -> the quasi-Toeplitz matrix H with
%                                 column-sign-randomization. In this case
%                                'para1' specifies B/N.
%                   'sparse'   -> the Bernoulli-sparse matrix whose
%                                 matrix sparsity is specified by 'para1'.
% written by Jaewook Kang 2015 Aug
% Feedback to jwkkang@gist.ac.kr 
%-------------------------------------------------------------------------
M=varargin{1};N=varargin{2};

if nargin < 3
    mxt_attribute='Gaussian';
else
    mxt_attribute=varargin{3};

    if strcmp(mxt_attribute,'Toeplitz')==1
        if nargin < 4
            para1=0.5;
        else
            para1=varargin{4};
        end
    end

    if strcmp(mxt_attribute,'sparse')==1
        if nargin < 4
            para1=0.05;
        else
            para1=varargin{4};
        end
    end
end


if strcmp(mxt_attribute,'Gaussian') ==1
    %1) i.i.d. sub-Gaussian matrices 
    H=randn(M,N)/sqrt(M);      
elseif strcmp(mxt_attribute,'Bernoulli')  ==1
    %2) i.i.d. Bernoulli matrices 
    H=(2*(rand(M,N)<0.5)-1)/sqrt(M);
elseif strcmp(mxt_attribute,'WHT')  ==1
    %3) Walsh-Hadamard matrix % only for when N is 2's power
    WHT_mat=fwht2(eye(N));
    H=WHT_mat(find(randerr(1,N,M)'),:)/sqrt(M/N)*diag(2*floor(2*rand(N,1))-1);
    H=H/norm(H(:,1));
elseif strcmp(mxt_attribute,'DCT') ==1
    %5) subsampled DCT matrix  with column sign randomization
    DCT_mat= dctmtx2(N);
    H=DCT_mat(find(randerr(1,N,M)'),:)/sqrt(M/N)*diag(2*floor(2*rand(N,1))-1);
elseif strcmp(mxt_attribute,'Toeplitz') ==1
    %6) quasi-Toeplitz with column sign randomization
    B=round(para1*N);% the number of taps
    first_row=zeros(N,1);
    first_row(1:B)=randn(B,1)/sqrt(M);
    tempH=zeros(N,M);
    for i=1:M
        tempH(:,i)=circshift(first_row,i-1);
    end
    H=tempH'*diag(2*floor(2*rand(N,1))-1);
elseif strcmp(mxt_attribute,'sparse')   ==1  
    %4) i.i.d. Bernoulli-sparse matrices 
    tempH=(2*floor(2*rand(M,N))-1)/sqrt(para1*M);% tempH \in {0,1,-1}
    sparse_H=rand(M,N)<para1;
    H=sparse(tempH.*sparse_H);
else
    H=randn(M,N)/sqrt(M); % defalut is Gaussian
end

end







                
