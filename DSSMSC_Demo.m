clear;
clc;
addpath(genpath('DSSMSC'))
load('BBCSport.mat');
% load('MSRCV1.mat');
v=size(data,2);
numClust = length(unique(gt));
%%%%%%%%%%%%%   standardization %%%%%%%%%%%%%
X = cell(1, v);
for i = 1:v
    X{i}=full(data{i});
    X{i} = NormalizeFea(X{i},0); 
end
%% Subspace segmentation
% metric = 'squaredeuclidean';
% metric = 'cosine';
Par.lambda_4 = 1;
Par.rho = 1;
Par.maxIter = 50;
knn0 = 5;
metric = 'cosine';    %text data
% metric = 'squaredeuclidean';  % image data
Par.lambda_1 =  0.06;  
Par.lambda_2 = 2;    
Par.lambda_3 = 0.8;    
Par.s = 0.7;       
% MSRCV1
% Par.lambda_1 = 0.6;  
% Par.lambda_2 = 0.1;    
% Par.lambda_3 = 0.3;   
% Par.s = 0.9;   
[Z,C,C1,S] = DSSMSC( X , Par , metric, knn0);  
S= (S +  S')/2;
label = SpectralClustering(S,numClust);
score = Clustering8Measure(label, gt)
%bbcsports
%0.9458    0.9624    0.9853    0.9853    0.9693    0.9734    0.9713
%MSRCV1
%0.8784    0.8579    0.9381    0.9381    0.8818    0.8737    0.8777
