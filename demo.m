% Revised by Chenyang according to the work of Li Chengbo
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

clear; close all;
% path(path,genpath(pwd));
% i = 89;
% ObjDir = 'E:\Dropbox\Dropbox\TVAL3\Simu_diff_p\DFT\chosen\';
% bgFile = [ObjDir,int2str(i),'.bmp'];
% I = im2double(imread(bgFile));

fullscreen = get(0,'ScreenSize');

%% size of I
I = imread('9.png');
I = im2double( I(220:320,250:350,1));
I = I/max(max(I))*256;
I = I -70;
I(I<0) = 0;
% nrmI = norm(I,'fro');
nrmI = norm(I,'fro');
[p,q] = size(I);
m = 3000;

%% A
load('Record.mat');
ratio = m/10000;
A = rand(m,p*q)-0.5;
% % 
A = zeros(m,p*q);
AA = eye(p*q);
for i = 1:1:10
% for i = 1:1:size(Record,2)
A = AA(Record(1:m,i),:) + A;
end
% clear AA;

f = A*I(:);
figure(2);
figure('Name','TVAL3','Position',...
    [fullscreen(1) fullscreen(2) fullscreen(3) fullscreen(4)]);
subplot(121); imshow(I,[]);
title('Original phantom','fontsize',18); drawnow;


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

t = cputime;
[U, out] =  ftvcs_alp(A,f,p,q,opts);
t = cputime - t;


subplot(122); 
imshow(U,[]);
title('Recovered by TVAL3','fontsize',18);
xlabel(sprintf(' %2d%% measurements \n Rel-Err: %4.2f%%, CPU: %4.2fs ',ratio*100,norm(U-I,'fro')/nrmI*100,t),'fontsize',16);

figure(2);
I2 = imresize(I, sqrt(ratio));
I3 = imresize(I2,1/sqrt(ratio), 'bilinear');
subplot(121); imshow(I3,[]);
subplot(122); imshow(U,[]);

figure(3);
plot(I(30,:),'b');
hold on;
plot(U(30,:),'r');
hold on;
plot(I3(30,:),'g');
legend('original','rebuild','bilinear')
% imwrite(uint8(I),'original.bmp','bmp');
% imwrite(uint8(U),'rebuild.bmp','bmp');
% imwrite(uint8(I3),'bilinear.bmp','bmp');
