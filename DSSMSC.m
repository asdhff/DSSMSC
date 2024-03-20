function [Z,C,C1,S] = DSSMSC( X , Par , metric, knn0)
%输入是p*n的
v = length(X);
fea = cell(1, v);
for i = 1:v
    fea{i} = X{i}';
end
WW = make_distance_matrix(fea, metric);
knn = knn0 + 1;

% Construct kernel
K = cell(1, v);
for i = 1:v
    K{i} = make_kNN_dist(WW{i}, knn);
end

L = laplacian(K);

n = size(X{1},2);

%mu = 1e-4;
mu = 10;
muMax=1e6;
iterMax=Par.maxIter;
tol  = 1e-6;
iter = 0;

Z = cell(1, v);
C = cell(1, v);
C1 = cell(1, v);
Y = cell(1, v);
S = zeros(n, n);
M = cell(1, v);

for i = 1:v
    M{i} = X{i}'*X{i};
    Z{i} = zeros(n, n);
    C{i} = zeros(n, n);
    C1{i} = zeros(n, n);
    Y{i} = zeros(n, n);
end
while  ( iter<iterMax )
     Z = Z_update(M, C, L, Par);
     C = C_update(Z,C1,Y,Par,mu);
     C1 = C1_update(C,S,Y,Par,mu);
     Y = Y_update(C,C1,Y,mu);
     S = S_update(C1,Par);
     mu=min(mu*Par.rho,muMax);
     iter = iter + 1;
end

end
