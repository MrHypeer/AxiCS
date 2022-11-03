% this code computes affine transform that translates points from
% DMD2DMDgalvo mode.
% f: pxl@DMDGalvo --> DMDFOV
% 20221024 LIU Gan

clear all
% in BF preview window, positive y is up, positive x is right, left-hand
% cord
xyoffset = [0 15;0 0;0 -15; 0 -10; 0 5]; % in um unit
pxlposi = DMDgalvoFOV2pxl(xyoffset);
dmdpoints2 = DMD2DMDgalvo(pxlposi);

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
