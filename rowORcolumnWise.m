close all; clear all

%% when reshape a vector, fill column first, so Matlab is column-wise
Mat = reshape(1:2:30,3,5);
disp(Mat)

%% when locationg an element from 2D matrix using linear index, Matlab is 
% also column-wise
Mat(2); % result = Mat(2,1) = 3;

