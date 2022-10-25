% Revised by Chenyang, LIU Gan according to the work of Li Chengbo
%_________________________________
% This simple demo examines if TVAL3 works normally. Please try more demos
% 
% in the "Demos" directory, which would show users what TVAL3 is capable of.
% 
% I: 64x64 phantom (real, two-dimentional)
% A: random matrix without normality and orthogonality (real)
% f: observation with/without noise (real)
%
% Written by: Chengbo Li
% Advisor: Prof. Yin Zhang and Wotao Yin
% CAAM department, Rice University
% 05/21/2009
% this code should accept the data to be processed, not to modify every
% time!!! The only input should be .tiff and corresponding recoring matrix
% 09/28/2022

close all; clear all;
%% load picture and normalization, load measurement matrix
% input: target file name
% output: normalized 3D data matrix (z,x,y)
tiffName = 'pollen1_crop_150x150x30.tif';
info = imfinfo(tiffName);
numPage = length(info);
widthPage = info(1).Width;
heightPage = info(1).Height;

img3D = zeros(numPage, heightPage, widthPage); % stora data into 3D array
for ii = 1:numPage
    img3D(ii,:,:) = imread(tiffName, ii);
end
minInt = min(min(min(img3D)));
maxInt = max(max(max(img3D)));
img3D = (img3D-minInt)./(maxInt-minInt)*255; 
img3D = permute(img3D, [2,1,3]);
% !!!img3D should be in format: x_pxl*z_pxl*y_pxl
[x_length, z_length, y_length] = size(img3D);
load('RecordAxialY150Z3030_DMD.mat'); % load measurement matrix


%% add two loops for params selecting
compRatio = [10,20]; % input ratio
fociNum = [3,4];
for ii_ratio = 1:length(compRatio)
    for ii_foci = 1:length(fociNum)
        words = ['CR:',num2str(compRatio(ii_ratio)),'; Foci:',num2str(fociNum(ii_foci)),' is running...'];
        disp(words)
        
        %% loop operation to process multiple z or other requirement
        % input: loopNum
        
        loopNum = x_length; % all data is looped along x-direction
        imgRaw = zeros(x_length,z_length,y_length); % format: x_length * z_length * y_length
        imgCS = zeros(x_length,z_length,y_length);
        t = cputime;
        for ii_loop = 1:loopNum 
%             fprintf('%dth loop is running...', ii_loop)
        %     zIdx = 1; xIdx = randi(130); yIdx = randi(18);
            I = squeeze(img3D(ii_loop,:,:)); 
            I = im2double(I(:,:));
            I = I/max(max(I))*256;
            I = I-10;
            I(I<0) = 0;
            % nrmI = norm(I,'fro');
            imgRaw(ii_loop,:,:) = I;
            nrmI = norm(I,'fro');
            [p,q] = size(I);
            m = int16(compRatio(ii_ratio)/100*y_length*z_length);
            
            %% define measurement vector
            ratio = m/y_length*z_length;
            A = zeros(m,p*q);
            AA = eye(p*q);
            for i = 1:1:fociNum(ii_foci)
            % for i = 1:1:size(Record,2)
                A = AA(RecordAxialY150Z3030_DMD(1:m,i),:) + A; % change with filename
            end
            % clear AA;
            
            f = A*I(:);
            %% Run TVAL3
            clear opts
            opts.mu = 2^12;
            opts.beta = 2^7;
            opts.upsilon = 2^8;
            opts.mu0 = 2^8;
            opts.beta0 = 2^2;
            opts.upsilon0 = 2^2;
            opts.Ohm = 2;
            opts.tol = 1E-4;
            opts.maxit = 400;
            opts.TVnorm = 1;
            opts.nonneg = true;
            opts.rate_ctn = 1.2;
            
            [U, out] =  ftvcs_alp(A,f,p,q,opts);
            imgCS(ii_loop,:,:) = U;
            % imwrite(uint8(I),'original.bmp','bmp');
            % imwrite(uint8(U),'rebuild.bmp','bmp');
            % imwrite(uint8(I3),'bilinear.bmp','bmp');
        end
        
        %% save generated data
        splitName = strsplit(tiffName,'.');
        saveName = [cell2mat(splitName(1)),'_AxiCS','Rtio',num2str(compRatio(ii_ratio)),'Foci',num2str(fociNum(ii_foci))];
        save(saveName,"imgCS","imgRaw")
        disp([saveName,' saved!'])
        t = cputime - t;
        fprintf('Time spent: %f.3 \n', t);
    end
end

