%% update C1 
function C1 = C1_update(C,S,Y,Par,mu)
v = size(C, 2);
C1 = cell(1,v);
for i = 1:v
    Q = (mu*C{i} - Y{i}+2*Par.lambda_4*S)/(Par.s*(2*Par.lambda_3+2*Par.lambda_4+mu));
    J = SimplexProj(Q');
    %J = proj_simplex(Q');
    C1{i} = Par.s*J';
end  
end