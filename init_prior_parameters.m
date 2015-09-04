%% Prior parameters initialization
% This initialization is based on the paper 

% J. Vila and P. Schniter, “Expectation-maximization Gaussian-mixture
% approximate message passing,” IEEE Trans. Signal Process., vol. 61, no.
% 19, pp. 4658-4672, Oct. 2013.

% This code is written by Jaewook Kang @ May 2015.

function [sigmapow0,q]=init_prior_parameters(y,H)

[M,N]=size(H);

z=1:40;

temp=(1+z.^2).*normcdf(-z,0,1) -z.*normpdf(z,0,1);
q=M/N * max (     (1- 2*N/M*temp)./(1+z.^2-2*temp)  ) ;
sigmapow0 = norm(y)^2*(1-0.01)/(q*norm(H,'fro')^2);

% disp(sprintf('The initialization of prior parameters: q_est=%1.4f,  sigmapow0_est=%1.4f\n',q,sigmapow0 ));
end
