%% update Y
function Y = Y_update(C,C1,Y,mu)
v = size(C, 2);
B = Y;
Y = cell(1,v);
n = size(C{1}, 2);
for i = 1:v
    Y{i} = B{i} + mu * (C1{i}-C{i});
end  
end