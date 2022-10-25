%%%%%%%%%%%
% This code use differnet metrics to evaluate the performance of CS 
% metrics including: 
% 1. PSNR

close all; clear all
%% load data
filename = [".\data\pollen2_crop_150x150x30_AxiCSRtio20Foci4.mat",".\data\pollen2_crop_150x150x30_AxiCSRtio20Foci3.mat"...
    ,".\data\pollen2_crop_150x150x30_AxiCSRtio10Foci4.mat",".\data\pollen2_crop_150x150x30_AxiCSRtio10Foci3.mat"];
file_num = length(filename);
for ii_file = 1:file_num
    %% calculate PSNR and relative error
    load(filename(ii_file))
    [z_length, x_length, y_length] = size(imgRaw);
    if ii_file == 1
        psnr_x = zeros(x_length,file_num);
        relativeError_x = psnr_x;
        psnr_z = zeros(z_length,file_num);
    end

    % calculate the PSNR along x-axis
    for ii_x = 1: size(imgCS,2)
        psnr_x(ii_x,ii_file) = psnr(squeeze(imgCS(:,ii_x,:)),squeeze(imgRaw(:,ii_x,:)),255);
        relativeError_x(ii_x,ii_file) = norm(squeeze(imgCS(:,ii_x,:) - imgRaw(:,ii_x,:)),"fro")...
            /norm(squeeze(imgRaw(:,ii_x,:)))*100; % in percentage
    end
    % calculate the PSNR along z-axis
    for ii_z = 1: size(imgCS,1)
        psnr_z(ii_z, ii_file) = psnr(squeeze(imgCS(ii_z,:,:)),squeeze(imgRaw(ii_z,:,:)),255);
    end
    
    
    %% print result
    % deal with Inf 
    temp = squeeze(psnr_x(:,ii_file));
    inf_falg = isinf(temp);
    one_flag = find(~inf_falg);
    fprintf('mean PSNR along x: %4.2f dB\n',sum(temp(one_flag))/length(one_flag))
    % deal with Nan
    temp = squeeze(relativeError_x(:,ii_file));
    nan_falg = isnan(temp);
    nan_flag = find(~nan_falg);
    fprintf('mean relative Error along x: %4.2f percent \n', sum(temp(nan_flag))/length(nan_flag))
    % deal with Inf 
    temp = squeeze(psnr_z(:,ii_file));
    inf_falg = isinf(temp);
    one_flag = find(~inf_falg);
    fprintf('mean PSNR along z: %4.2f dB \n',sum(temp(one_flag))/length(one_flag))
end

%% plot result
figure(1)
subplot(1,2,1)
plot(psnr_x)
title('Ratio and Foci vs. PSNR')
xlabel('x plane number'); ylabel('PSNR/dB')
legend({'CS20Foci4','CS20Foci3','CS10Foci4','CS10Foci3'})
subplot(1,2,2)
plot(relativeError_x)
title('Ratio and Foci vs. relative error')
xlabel('x plane number'); ylabel('relative error/%')
legend({'CS20Foci4','CS20Foci3','CS10Foci4','CS10Foci3'})
sgtitle('pollen2-tulip pollen')
