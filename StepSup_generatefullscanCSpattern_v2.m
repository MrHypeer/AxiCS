%%% This code generate the full scan area that is the same as the CS 
% pattern, which can compare the full scan ability with CS performace.
% 
% add portion which suits DMDgalvo scanning pattern
% this code cannot generate right DMD pattern, replace pattern generation
% code in v3

% close all;
clear all

%% define params
stepsize = 0.4; Focinum = 1;
stepsize_zscale = 2.5;

%% load point cell from saved data
record_name = 'RecordPointPair_X1Y100Z30XY04Z04_CS20Foci3Dist2.mat';
data_path = '.\experiment20221024\yzTestZ0\';
Pair = load(strcat(data_path,record_name));
% Pair = load(record_name);
Point = Pair.Point(1:Focinum,1);

%% define cell and pack data position into it
% Define scan area and pixel number
% pxl_x = 1; pxl_y = 100; pxl_z = 30;
% Point = cell(1,1); 
% 
% x = floor(1-pxl_x/2:1:pxl_x/2);
% y = floor(1-pxl_y/2:1:pxl_y/2);
% z = floor(0:pxl_z-1);
% 
% point_idx = 1;
% 
% for ii_z = 1:pxl_z
%     for ii_x = 1:pxl_x
%         for ii_y = 1:pxl_y
%             Point{1}(point_idx,1) = x(ii_x);
%             Point{1}(point_idx,2) = y(ii_y);
%             Point{1}(point_idx,3) = z(ii_z);
%             point_idx = point_idx+1;
%         end
%     end
% end


%% generate hologram
% scale Point define from original matrix into DMDgalvo position
scale = 30/100;
xyoffset = Point{1}(:,1:2)*scale;

% calculate DMDGalvo pxl positions
pxlposi = DMDgalvoFOV2pxl(xyoffset);
DMDpoints_xy = DMD2DMDgalvo(pxlposi);
Point{1} = [DMDpoints_xy, Point{1}(:,3)*0];

% for ii = 1:Focinum
%     Point{ii}(:,1:2) = stepsize * Point{ii}(:,1:2);
%     Point{ii}(:,3) = stepsize * stepsize_zscale * Point{ii}(:,3);
% end
tic
[Uniformity,I] =  TrajectoryBGS_Multi_40gpu(1, '.\experiment20221024\yzTestZ0\', 'Partialscan_X1Y100Z30XY04Z04_CS20Foci1', Point{1:Focinum});
toc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DMDpoints = DMD2DMDgalvo(DMDGalvopoints)
% this code computes affine transform that translates points from
% DMD2DMDgalvo mode.
% f: pxl@DMDGalvo --> DMDFOV(in um)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NOTE: ONLY ACCEPTS THREE POINTS EACH TIME
% input: DMDgalvo point in pxl indics
    
    % input data points obtained from imageJ, change this accordingly if
    % light path was changed
    % list FOV point in DMD path
    % (45,0,0)(45,20,0)(10,20,0)
    xFOV = [45 45 10]';
    yFOV = [0 20 20]';
    
    % list DMDGalvo point data,order: y = 20,0,-20
    uDMDGalvo = [1209.625 1313.961 1276.691]';
    vDMDGalvo = [414.437 437.75 629.606]';
    
    % list DMD point data,order: y = 20,0,-20
    uDMD = [1493.682 1499.24 1303.379]';
    vDMD = [771.236 663.923 652.896]';
    
    % flip position along horizon axis for DMDGalvo data
    % as DMD and DMDGalvo path is mirror version, so flip before transform
    uDMDGalvo_flip = repmat(1920, 3,1)-uDMDGalvo;
    
    % calculate the tranform matrix
    DMD2Galvoflip = maketform('affine',[uDMD vDMD],[uDMDGalvo_flip vDMDGalvo]);
    FOV2DMD = maketform('affine', [xFOV yFOV],[uDMD vDMD]);
    
    % flip input data for process
    x_tar_flip = repmat(1920, size(DMDGalvopoints,1),1)-DMDGalvopoints(:,1);
    y_tar = DMDGalvopoints(:,2);
    [uDMD, vDMD] = tforminv(DMD2Galvoflip, x_tar_flip, y_tar);
    DMDpoints = tforminv(FOV2DMD, uDMD, vDMD);
end

function pxlposi = DMDgalvoFOV2pxl(xyoffset)
% this code transforms brightfield postion in micrometer to pixel
% indicis in brigt field camera (1920*1200)
% input: offset from center (45,0,0) in um along camera sensor edge, dim = 3x2
%     pxl_hor = 1900; pxl_ver = 1200;
    scale = 0.1852;% um/pxl, a round value, can be precisely measured
    xCenter_DMDGalvo = 1209.555; yCenter_DMDGalvo = 414.733; % measured at (45,0,0)
    
    pxlposi = zeros(size(xyoffset,1),2);
    pxlposi(:,1) = xCenter_DMDGalvo + xyoffset(:,1)/scale;
    pxlposi(:,2) = yCenter_DMDGalvo + (-1*xyoffset(:,2))/scale; % BF saved data is upside down
end