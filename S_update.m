%% update S
function S = S_update(C1,Par)
v = size(C1, 2);
n = size(C1{1}, 2);
G = zeros(n, n);
for i = 1:v
    G = G + C1{i};
end  
M = G/(Par.s*v);
J = SimplexProj(M');
%J = proj_simplex(M');
S = Par.s*J';
end