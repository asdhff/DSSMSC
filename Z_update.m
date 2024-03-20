%% update Z the coefficient matrix
function Z = Z_update(M, C, L, Par)
v = size(M, 2);
Z = cell(1,v);
n = size(M{1}, 2);
for i = 1:v
    Ta=M{i}+Par.lambda_2*eye(n);
    Tb=Par.lambda_2*(C{i}*C{i}'-C{i}-C{i}')+ Par.lambda_1*L{i};
    Tc=M{i};
    Z{i}=sylvester(Ta,Tb,Tc); 
end  
end