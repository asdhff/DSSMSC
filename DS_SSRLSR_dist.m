function [J,Z] =DS_SSRLSR_dist( X , Par, knn0, metric )
%X的横为N个样本数，纵为D维度
K=X'*X;
mu = 1e-4;
muMax=1e6;
iterMax=Par.maxIter;
[~,N] = size (X);
C = zeros (N, N);
J  = C;
Delta = zeros (N, N); % Y
Z=zeros(N,N);

fea{1} = X';
WW = make_distance_matrix(fea, metric);
knn = knn0 + 1;
% Construct kernel
M{1} = make_kNN_dist(WW{1}, knn);
L = laplacian(M);
%%
tol   = 1e-6;
iter    = 0;
%     err1(1) = inf; err2(1) = inf; err3(1) = inf; err4(1)=inf;
while  ( iter<iterMax )
%     Cpre = C;
%     Jpre = J;
%     Zpre = Z;
%    
 %% update C the coefficient matrix
    C=inv(2*Par.lambda_2*(Z'*Z)+mu*eye(N))*(Delta+mu*J+2*Par.lambda_2*Z'*Z);
     
    %% update J the data term matrix
    Q = (mu*C - Delta)/(Par.s*(2*Par.lambda_3+mu));
    J = SimplexProj(Q');
    J = Par.s*J';
    
    %%更新Z
    Ta=K+Par.lambda_2*eye(N);
    Tb=Par.lambda_2*(C*C'-C-C')+Par.lambda_1*L{1};
    Tc=K;
    Z=sylvester(Ta,Tb,Tc); 
    %% update Deltas the lagrange multiplier matrix
    Delta = Delta +mu * (J-C);
    mu=min(mu*Par.rho,muMax);
 
    %% computing errors
    iter = iter + 1;
%       if (  (err1(iter+1) <= tol && err1(iter+1)<=tol && err2(iter+1)<=tol && err3(iter+1)<=tol&& err4(iter+1)<=tol) || iter >= Par.maxIter )
%              terminate = true;
%             fprintf('  iter=%d\n',iter);
%       end
end

end
