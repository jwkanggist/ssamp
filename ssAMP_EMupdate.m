%% EM update for prior parameters in ssAMP-1D 
%% written by Jaewook Kang @ May 2015.
function [q_next,sigmapow0_next]=ssAMP_EMupdate(rho,theta,q_curr,sigmapow0_curr)

    N=length(rho);
    sum_var  = 1/sigmapow0_curr + 1/2./theta;
    diff_rho = diff(rho);
    
    pi_k    =  ( 1 +  (1-q_curr).*normpdf(0        , diff_rho, sqrt(2*theta)                 ) ./ ... 
                      (  q_curr .*normpdf(diff_rho , 0, sqrt(2*theta +sigmapow0_curr)))   ).^-1;
             

    gamma_k = ( diff_rho./2./theta) ./ sum_var;
    nu = 1/sum_var;
    
%     q_next         = 0.5*q_curr + 0.5* sum(pi_k)/(N-1);
%     sigmapow0_next = 0.5*sigmapow0_curr + 0.5*sum(pi_k.* (abs(gamma_k).^2+nu) )/q_next/(N-1);   
% 
    q_next         = sum(pi_k)./(N-1);
    sigmapow0_next = sum(pi_k.* (abs(gamma_k).^2+nu) )./q_next/(N-1);      
end
