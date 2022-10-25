%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  This demo show the ability of CS along axial direction
%%%  The key issue here is the requirement for pixel number 
%%%  along axial direction.
%%%  Refer to LI Chengbo and WEN Chenyang for help
%%%  working flow of this code is:
%%%  1. Load raw data and normalize the data
%%%  2. select the ROI where contains most structures manually
%%%  3. slice the raw data and apply CS algorithm on y-z plane along x-axis


close all; clear all;
%% load rawdata and normalization
tiffName = 'pollen1_crop_150x150x30.tif';
info = imfinfo(tiffName);
numPage = length(info); % get the stack number of stack 
widthPage = info(1).Width; % get single frame's width and height information
heightPage = info(1).Height;

img3D = zeros(numPage, widthPage, heightPage); % stora data into 3D array
for ii = 1:numPage
    img3D(ii,:,:) = imread(tiffName, ii);
end

minInt = min(min(min(img3D)));
maxInt = max(max(max(img3D)));
img3D_Norm = uint8((img3D-minInt)./(maxInt-minInt).*255); % normalize the image
[zPxl, xPxl, yPxl] = size(img3D);

%% select the ROI manually
% sliceLat = squeeze(img3D_Norm(15,:,:));
% sliceAxi = squeeze(img3D_Norm(:,100,:));
% 
% % prview the slice
% figure; 
% subplot(1,2,1)
% imshow(sliceLat); title('Lateral cross section')
% subplot(1,2,2)
% imshow(sliceAxi); title('Axial cross section')

%% apply axial CS algorithm along x-axis

% slice a data volume with (150*150*30)
img3D_CS = zeros(zPxl,xPxl,yPxl); % save CS-recovered data
for ii_x = 1:150
    %% choose the CS ROI
    I = double(squeeze(img3D_Norm(:, ii_x, :)));% dim: 33*512 pxls
    I = I - 40;
    I(I<0) = 0;
    nrmI = norm(I,'fro');
    [p,q] = size(I);


    %% generate the measurement matrix
    % RecordAxial = randi(33*512,5000,30);
    record_file = 'RecordAxialY150Z3030_DMD.mat';
    load(record_file); % 4500 samples and 30 points every time
    m = 450; % sample number
    foci = 3;
    ratio = m/(p*q);
    A = zeros(m,p*q);
    AA = eye(p*q);
    for i = 1:1:foci
    % for i = 1:1:size(Record,2)
        A = AA(RecordAxialY150Z3030_DMD(1:m,i),:) + A;
    end
    % clear AA;

    %% reconstrcut
    t = cputime;

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
    t = cputime - t;

    img3D_CS(:,ii_x,:) = U;
    disp(t)
end
%% preview and compare the result
% figure(2);
% I2 = imresize(I, sqrt(ratio));
% I3 = imresize(I2,1/sqrt(ratio), 'bilinear');
% subplot(121); imshow(I3,[]);
% subplot(122); imshow(U,[]);

% figure(3);
% plot(I(30,:),'b');
% hold on;
% plot(U(30,:),'r');
% hold on;
% plot(I3(30,:),'g');
% legend('original','rebuild','bilinear')
% imwrite(uint8(I),'original.bmp','bmp');
% imwrite(uint8(U),'rebuild.bmp','bmp');
% imwrite(uint8(I3),'bilinear.bmp','bmp');


%% statistics
nrmI = norm(I,'fro');
deviImg = norm(U-I,'fro')/nrmI*100;

%% plot the result for a single slice
% figure(1)
% set(gcf, 'PaperSize', [4 2]);
% subplot(1,2,1)
% imshow(squeeze(img3D_Norm(:, ii_x, :)))
% title('original')
% subplot(1,2,2)
% imshow(U,[])
% title('reconstructed')


