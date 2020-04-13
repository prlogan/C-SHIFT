function [ Cov_new, Corr_new] = fCShift(Cov_obs, geneNo)
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
% DESCRIPTION: This function normalizes observed (noisy) covariance matrix
% using C-SHIFT method 
%
% INPUT:
% Cov_obs - observed covariance matrix
% geneNo - number of genes
%
% OUTPUT:
% Cov_new - normalized covariance matrix
% Corr_new  - normalized correlation matrix
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% More information about build-in fucntion used: 
% https://www.mathworks.com/help/optim/ug/fminunc.html
% https://www.mathworks.com/help/optim/ug/fmincon.html#busog7r-fun
% https://www.mathworks.com/help/optim/ug/fmincon.html#busog7r-options
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% if Cov_obs is not a full rank then we need to add additional diagonal perturbation matrix for better inversion
if rank(Cov_obs)==geneNo
    Cov_pert = Cov_obs;
else
    F = diag(rand(1,geneNo));
    Cov_pert = Cov_obs+F;
end

alpha = zeros(geneNo,1); % initial alpha

% set the current algorithm to 'trust-region' (defualt option is 'quasi-newton')
options = optimoptions('fminunc','Algorithm','trust-region','SpecifyObjectiveGradient',true,'HessianFcn','objective');
f = @ObjFunct;


    function [obj_funct,grad_obj_funct,hess] = ObjFunct(alpha)
        
        D = repmat(alpha,1,geneNo);
        M_hlp = Cov_pert+D+D';
        M_inv_hat = M_hlp\ones(geneNo,1);
        V = 1/(ones(1,geneNo)*M_inv_hat);
        % function
        obj_funct = norm(M_hlp-V*ones(geneNo,geneNo),'fro')^2;
        % gradient
        a = sum(alpha);
        c = sum(sum(Cov_pert));
        hlp1 = geneNo^2*V^2-c*V-2*geneNo*a*V;
        grad_obj_funct = 4*(geneNo*alpha + hlp1*M_inv_hat + sum(Cov_pert,2) + (a-geneNo*V)*ones(geneNo,1));
        
        H0 = M_inv_hat*ones(1,geneNo);
        H1 = H0+H0';
        H2 = M_inv_hat*M_inv_hat';
        hess = 4*(diag(geneNo*ones(1,geneNo))+ones(geneNo,geneNo)-2*geneNo*V*H1+(3*geneNo^2*V-c-2*geneNo*a)*V*H2-(geneNo^2*V-c-2*geneNo*a)*inv(M_hlp));  
        
        if norm(hess-hess','fro')<0.01
            disp('hess is symmetric')
            disp('frobnorm(Hess-hess')
            disp(norm(hess-hess','fro'))
        else
            disp('hess is not symmetric')
            disp('frobnorm(Hess-hess')
            disp(norm(hess-hess','fro'))
        end
        
    end


alpha_sol = fminunc(f,alpha,options);

% calculate the solution covariance matrix Cov_new that depends on a newly
% found alpha_sol
D = repmat(alpha_sol,1,geneNo);
M = Cov_pert+D+D';
M_inv_hat = M\ones(geneNo,1);
V = 1/(ones(1,geneNo)*M_inv_hat);
C_alpha_sol = M-V*ones(geneNo,geneNo);

Cov_new = C_alpha_sol-F;

Corr_new = zeros(geneNo,geneNo);
for m = 1:geneNo
    for n = m:geneNo
        Corr_new(m,n) = (Cov_new(m,n))/(sqrt((Cov_new(n,n))*(Cov_new(m,m))));
        Corr_new(n,m) = Corr_new(m,n);
    end
end
end


