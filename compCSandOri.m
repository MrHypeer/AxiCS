close all;clear all

matName = 'pollen1_crop_150x150x30_AxiCSRtio20Foci4.mat';

file = load(matName);
tiffCS = file.imgCS; %CS result
tiffOri = file.imgRaw; % original pic for CS reconstruction
[xPxl, zPxl, yPxl] = size(tiffCS);


%% calculate error
errorCs = zeros(xPxl,1);
for ii_x = 1:xPxl
    I = squeeze(tiffOri(ii_x,:,:));
    U = squeeze(tiffCS(ii_x,:,:));
    nrmI = norm(I,'fro');
    errorCs(ii_x) = norm(U-I,'fro')/nrmI;
end

% figure; plot(errorCs);title('error statistic')
meanError = mean(errorCs)
varError = var(errorCs);


%% plot the original and CS iamge
% [~, zIdx] = max(errorCs,[],'linear');
% figure
% subplot(121)
% tiffCS2 = uint8(squeeze(tiffCS(zIdx,:,:)));
% imshow(tiffCS2)
% title('CS pic')
% subplot(122)
% tiffOri2 = uint8(squeeze(tiffOri(zIdx,:,:)));
% imshow(tiffOri2)
% title('Original pic')

%% save CS file into tiff format, 
tiffCS = uint8(permute(tiffCS,[2,3,1]));
tiffOri = uint8(permute(tiffOri,[2,3,1]));
splitName = strsplit(matName,'.');

% save CS file
targetfile = [cell2mat(splitName(1)),'_CS.tif'];
t = Tiff(targetfile,'w');
tagstruct.ImageLength = size(tiffCS,1);
tagstruct.ImageWidth = size(tiffCS,2);
tagstruct.SampleFormat = 1; % uint
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 8;
tagstruct.SamplesPerPixel = 1;
tagstruct.Compression = Tiff.Compression.Deflate;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
for ii=1:size(tiffCS,3)
   setTag(t,tagstruct);
   write(t,tiffCS(:,:,ii));
   writeDirectory(t);
end
close(t)

% save original NormSub image 
targetfile = [cell2mat(splitName(1)),'_NormSub.tif'];
t = Tiff(targetfile,'w');
tagstruct.ImageLength = size(tiffOri,1);
tagstruct.ImageWidth = size(tiffOri,2);
tagstruct.SampleFormat = 1; % uint
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = 8;
tagstruct.SamplesPerPixel = 1;
tagstruct.Compression = Tiff.Compression.Deflate;
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
for ii=1:size(tiffOri,3)
   setTag(t,tagstruct);
   write(t,tiffOri(:,:,ii));
   writeDirectory(t);
end
close(t)