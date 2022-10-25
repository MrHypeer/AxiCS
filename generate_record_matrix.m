% this code generates the record matrix for CS algorithm
% record data format: (num_sampling, max_foci)

close all; clear all;
y_pxl = 150; z_pxl = 30;

RecordAxialY150Z3030 = randi([1,y_pxl*z_pxl],y_pxl*z_pxl,30);