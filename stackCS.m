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

close all; clear all;
%% load picture and normalization, load measurement matrix
% input: target file name
% output: normalized 3D data matrix (z,x,y)
tiffName = 'neuron-stack-330um_scale01.tif';
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

load('RecordZ100Y201.mat'); % load measurement matrix

%% select the region of interest manually
% input: vIdx/hIdx 
% output: satisfying zIdx, xIdx, yIdx
zIdx = 50; xIdx = 2.2; yIdx = 2.1; % for neuron cell
% zIdx = 1; xIdx = randi(130); yIdx = randi(18);
xidx = int16(xIdx*100); yidx = int16(yIdx*100);
I = squeeze(img3D(zIdx,xidx:xidx+200, yidx:yidx+200)); % select a region cover 200x200 pixel
I = im2double(I(:,:));
I = I/max(max(I))*256;
I = I -10;
I(I<0) = 0;
figure; imagesc(I); title('region selection')


%% loop operation to process multiple z or other requirement
% input: loopNum

loopNum = 201;
imgRaw = zeros(loopNum,100,201);
imgCS = zeros(loopNum,100,201);
t = cputime;
for ii_loop = 1:loopNum 
    fprintf('%dth loop is running...', ii_loop)
%     zIdx = 1; xIdx = randi(130); yIdx = randi(18);
    xidx = int16(xIdx*100); yidx = int16(yIdx*100);
    I = squeeze(img3D(1:100,xidx+ii_loop, yidx:yidx+200)); 
    I = im2double(I(:,:));
    I = I/max(max(I))*256;
    I = I -10;
    I(I<0) = 0;
    % nrmI = norm(I,'fro');
    imgRaw(ii_loop,:,:) = I;
    nrmI = norm(I,'fro');
    [p,q] = size(I);
    m = 5000;
    
    %% define measurement vector
    ratio = m/20100;
    A = zeros(m,p*q);
    AA = eye(p*q);
    for i = 1:1:10
    % for i = 1:1:size(Record,2)
        A = AA(RecordZ100Y201(1:m,i),:) + A;
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
saveName = [cell2mat(splitName(1)),'_AxiCS'];
save(saveName,"imgCS","imgRaw")
t = cputime - t;
fprintf('Time spent: %f.3', t);


