function objvalue = objectvalue(X, L, Z, C, S, Par)
v = length(X);
objvalue = 0;
for i = 1:v
    F1 = norm(X{i}-X{i}*Z{i},'fro')^2;
    F2 = Par.lambda_1*trace(Z{i}*L{i}*Z{i}');
    F3 = Par.lambda_2*norm(Z{i}-Z{i}*C{i},'fro')^2;
    F4 = Par.lambda_3*norm(C{i},'fro')^2;
    F5 = Par.lambda_4*norm(C{i}-S,'fro')^2;
    objvalue = objvalue+F1+F2+F3+F4+F5;
end
end
    
