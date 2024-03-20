%% update C the coefficient matrix
function C = C_update(Z,C1,Y,Par,mu)
v = length(Z);
C = cell(1,v);
n = size(Z{1},2);
for i = 1:v
    C{i} = (2*Par.lambda_2*(Z{i}'*Z{i})+mu*eye(n))\(Y{i}+mu*C1{i}+2*Par.lambda_2*Z{i}'*Z{i});
end  
end
  