%%% This code generate the full scan area that is the same as the CS 
% pattern, which can compare the full scan ability with CS performace.


close all; clear all

%% define params
stepsize = 0.4; Focinum = 1;
stepsize_zscale = 2.5;

%% load point cell from saved data
% record_name = 'RecordPointPair_X1Y100Z30_CS60Foci2Dis5.mat';
% data_path = '.\experiment20221007\yzTest2\result\';
% Pair = load(strcat(data_path,record_name));
% % Pair = load(record_name);
% Point = Pair.Point(1:Focinum,1);

%% define cell and pack data position into it
% Define scan area and pixel number
pxl_x = 1; pxl_y = 100; pxl_z = 30;
Point = cell(1,1); 

x = floor(1-pxl_x/2:1:pxl_x/2);
y = floor(1-pxl_y/2:1:pxl_y/2);
z = floor(0:pxl_z-1);

point_idx = 1;

for ii_z = 1:pxl_z
    for ii_x = 1:pxl_x
        for ii_y = 1:pxl_y
            Point{1}(point_idx,1) = x(ii_x);
            Point{1}(point_idx,2) = y(ii_y);
            Point{1}(point_idx,3) = z(ii_z);
            point_idx = point_idx+1;
        end
    end
end


%% generate hologram
for ii = 1:Focinum
    Point{ii}(:,1:2) = stepsize * Point{ii}(:,1:2);
    Point{ii}(:,3) = stepsize * stepsize_zscale * Point{ii}(:,3);
end
tic
[Uniformity,I] =  TrajectoryBGS_Multi_40gpu(50, '.\experiment20221011\yzTest\', 'Fullscan_X1Y100Z30XY04Z10', Point{1:Focinum});
toc