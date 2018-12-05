function [A_hat,E_hat] = APG_IR(D, lambda, maxIter, tol, ...
            lineSearchFlag, continuationFlag, eta, mu, outputFileName)
% This code is an implementation of Proximal Gradient Algorithm for small target
% detecton in our published paper: Chenqiang Gao, Deyu Meng, Yi Yang, et al., "Infrared Patch-Image Model for Small Target Detection in a Single Image, " Image Processing, IEEE Transactions on, vol. 22, no. 12, pp. 4996-5009, 2013.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you use this code in your publications, please cite:
% Chenqiang Gao, Deyu Meng, Yi Yang, et al., "Infrared Patch-Image Model for Small Target Detection in a Single Image, " Image Processing, IEEE Transactions on, vol. 22, no. 12, pp. 4996-5009, 2013. 
% @article{Gao2013,
%    author = {Gao, Chenqiang and Meng, Deyu and Yang, Yi and Wang, Yongtao and Zhou, Xiaofang and Hauptmann, Alex},
%    title = {Infrared Patch-Image Model for Small Target Detection in a Single Image},
%    journal = {Image Processing, IEEE Transactions on},
%    volume = {22},
%    number = {12},
%    pages = {4996-5009},
%    year = {2013}
% }
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code is implemented based on the original code of Accelerated
% Proximal Gradient algorithm on the websie:
% "http://perception.csl.illinois.edu/matrix-rank/sample_code.html". So
% please also cite the following reference:
% "Robust PCA: Exact Recovery of Corrupted Low-Rank Matrices via Convex Optimization", J. Wright et al., preprint 2009.
% "An Accelerated Proximal Gradient Algorithm for Nuclear Norm Regularized Least Squares problems", K.-C. Toh and S. Yun, preprint 2009.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Paremeter Explanation
%
% D - m x n matrix of observations/data (required input)
% lambda - weight on sparse error term in the cost function (required input)
%
% tol - tolerance for stopping criterion.
%     - DEFAULT 1e-7 if omitted or -1.
% maxIter - maximum number of iterations
%         - DEFAULT 10000, if omitted or -1.
% lineSearchFlag - 1 if line search is to be done every iteration
%                - DEFAULT 0, if omitted or -1.
% continuationFlag - 1 if a continuation is to be done on the parameter mu
%                  - DEFAULT 1, if omitted or -1.
% eta - line search parameter, should be in (0,1)
%     - ignored if lineSearchFlag is 0.
%     - DEFAULT 0.9, if omitted or -1.
% mu - relaxation parameter
%    - ignored if continuationFlag is 1.
%    - DEFAULT 1e-3, if omitted or -1.
% outputFileName - Details of each iteration are dumped here, if provided.
%
% [A_hat, E_hat] - estimates for the low-rank part and error part, respectively
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If you have any questions, please contact:
% Author: Chenqiang Gao 
% Email: gaochenqiang@gmail.com or gaocq@cqupt.edu.cn
% Copyright:  Chongqing University of Posts and Telecommunications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    error('Too few arguments') ;
end

if nargin < 4
    tol = 1e-6 ;
elseif tol == -1
    tol = 1e-6 ;
end

if nargin < 3
    maxIter = 1000 ;
elseif maxIter == -1
    maxIter = 1000 ;
end

if nargin < 6
    continuationFlag = 1 ;
elseif continuationFlag == -1 ;
    continuationFlag = 1 ;
end

if ~continuationFlag
    if nargin < 8
        mu = 1e-3 ;
    elseif mu == -1
        mu = 1e-3 ;
    end
end
% mu = 1e6 ;

DISPLAY_EVERY = 20 ;
DISPLAY = 0 ;

%% Initializing optimization variables

[m, n] = size(D) ;
t_k = 1 ; % t^k
t_km1 = 1 ; % t^{k-1}

tau_0 = 2 ; % square of Lipschitz constant for the RPCA problem

X_km1_A = zeros(m,n) ; X_km1_E = zeros(m,n) ; % X^{k-1} = (A^{k-1},E^{k-1})
X_k_A = zeros(m,n) ; X_k_E = zeros(m,n) ; % X^{k} = (A^{k},E^{k})


if continuationFlag
    mu_0 = norm(D) ;
    mu_k = 0.99*mu_0 ;
    mu_bar = 1e-9 * mu_0 ;
else
    mu_k = mu ;
end
[U S V] = svd(D, 'econ');
s = diag(S);
mu_k = s(2);
mu_bar = 0.005 * s(4);

tau_k = tau_0 ;

converged = 0 ;
numIter = 0 ;

NOChange_counter = 0;
pre_rank = 0;
pre_cardE = 0;
% tol = 1e-6 ;

%% Start main loop
while ~converged
    
    Y_k_A = X_k_A + ((t_km1 - 1)/t_k)*(X_k_A-X_km1_A) ;
    Y_k_E = X_k_E + ((t_km1 - 1)/t_k)*(X_k_E-X_km1_E) ;
    
    G_k_A = Y_k_A - (1/tau_k)*(Y_k_A+Y_k_E-D) ;
    G_k_E = Y_k_E - (1/tau_k)*(Y_k_A+Y_k_E-D) ;
    
    [U S V] = svd(G_k_A, 'econ');
    diagS = diag(S) ;
%      temp = max(diagS(150), mu_k/tau_k);
    
%     X_kp1_A = U * diag(pos(diagS - temp)) * V';
    X_kp1_A = U * diag(pos(diagS- mu_k/tau_k)) * V';
            
%     X_kp1_E = sign(G_k_E) .* pos( abs(G_k_E) - lambda* temp );
    X_kp1_E = sign(G_k_E) .* pos( abs(G_k_E) - lambda* mu_k/tau_k );
            
%   rankA  = sum(diagS>temp);
    rankA  = sum(diagS>mu_k/tau_k);
            
    cardE = sum(sum(double(abs(X_kp1_E)>0)));
        
    t_kp1 = 0.5*(1+sqrt(1+4*t_k*t_k)) ;
    
    temp = X_kp1_A + X_kp1_E - Y_k_A - Y_k_E ;
    S_kp1_A = tau_k*(Y_k_A-X_kp1_A) + temp ;
    S_kp1_E = tau_k*(Y_k_E-X_kp1_E) + temp ;
    
    
    stoppingCriterion = norm([S_kp1_A,S_kp1_E],'fro') / (tau_k * max(1, norm([X_kp1_A, X_kp1_E],'fro'))) ;
    
    if stoppingCriterion <= tol
        converged = 1 ;
    end
    
    if continuationFlag
        mu_k = max(0.9*mu_k, mu_bar) ;
    end
    
    t_km1 = t_k ;
    t_k = t_kp1 ;
    X_km1_A = X_k_A ; X_km1_E = X_k_E ;
    X_k_A = X_kp1_A ; X_k_E = X_kp1_E ;
    
    numIter = numIter + 1 ;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % The iteration process can be finished if the rank of A keeps the same
    % many times;
    if pre_rank == rankA
        NOChange_counter = NOChange_counter +1;
        if NOChange_counter > 10 && abs(cardE-pre_cardE) < 20
            converged = 1 ;
        end
    else
        NOChange_counter = 0;
        pre_cardE = cardE;
    end
    pre_rank = rankA;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % In practice, the APG algorithm, sometimes, cannot get  
    % a strictly low-rank matrix A_hat after iteration process. Many
    % singular valus of the obtained matrix A_hat, however, are extremely small. This can be considered 
    % to be low-rank to a certain extent. Experimental results show that the final recoverd 
    % backgournd image and target image are good.
    % Alternatively, we can make the rank of matrix A_hat lower using the following truncation. 
    % This trick can make the APG algorithm faster and the performance of our algorithm is still satisfied. 
    % Here we set the truncated threshold as 0.3, while it can be adaptively set based on your actual scenarios.
    if rankA > 0.3 * min([m n]) 
        converged = 1 ;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    if DISPLAY && mod(numIter,DISPLAY_EVERY) == 0
        disp(['Iteration ' num2str(numIter) '  rank(A) ' num2str(rankA) ...
            ' ||E||_0 ' num2str(cardE)])
    end
    
    if nargin > 8
        fprintf(fid, '%s\n', ['Iteration ' num2str(numIter) '  rank(A)  ' num2str(rankA) ...
            '  ||E||_0  ' num2str(cardE) '  Stopping Criterion   ' ...
            num2str(stoppingCriterion)]) ;
    end
        
    if ~converged && numIter >= maxIter
        disp('Maximum iterations reached') ;
        converged = 1 ;
    end
    
end

A_hat = X_k_A ;
E_hat = X_k_E ;
