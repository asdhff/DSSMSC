function D = getdist(X, metric, phi)
%X是一个cell，p*n的形式， metric可选 gauss, cosine
switch metric
    case 'gauss'
        D = dist_gauss(X, phi)
    case 'cosine'
        D = dist_cos(X)
end

function D = dist_gauss(X, phi)
v = size(X, 2)
n = size(X{1},2)
data = zeros(n, n)
for i = 1:n
    for j = 1:n
        x = X[:,i]
        y = X[:,j]
        
