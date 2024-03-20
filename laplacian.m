function L = laplacian(A)
% Fixed by Youwei. 2020.06.17
v = size(A, 2);
L = cell(1, v);
for i = 1:v
    degs = sum(A{i}, 2);
    D = sparse(1:size(A{i}, 1), 1:size(A{i}, 2), degs);
    L{i} = full(D - A{i});
end
end

