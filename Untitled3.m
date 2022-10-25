%% 利用LEE hologram中的方法三控制强度，生成多点
%% 生成点的均匀度不如binary rounding 的方法

% clc;
% clear all;
m = 1024/2; n = 768;

wl = 780e-9;
px_size = 13.68e-6;
fl = 200e-3;
Mobj = 40.509;
[row,col] = meshgrid(1:m,1:n);
[row2,col2] = meshgrid(1:2*m,1:n);

t = 16;
Points = zeros(4,3);
Points(1,:) = Point{1}(t,:);
Points(2,:) = Point{2}(t,:);
Points(3,:) = Point{3}(t,:);
Points(4,:) = Point{4}(t,:);
% 
% Points = [3.9,-3.8,0;...
%          -1.2,15,0;...
%          -11.4,-10.8,0;...
%          -9.6,-6.9,0;];
FocusNum = size(Points,1);

theta = zeros(n,m);
theta_test = zeros(n,m*2,FocusNum);
V = zeros(FocusNum,1);

% original point where x=0, y=0, z=0
% u, v, P are x and y diretion's spatial frequency of DMD pattern, and  optical power, respectively
u0 = 0.1;
v0 = 0;
P0 = 0;

for t=1:1:FocusNum
            x = Points(t,1)+ 38;
            y = Points(t,2);

           
 % accordingly parameters under optical settng: 54 mm scan lens; 175 mm tube
% lens, 25X water immersion
              u = x/199.05;
              v = y/199.05;

            % desired grating spatial frequency for X and Y direction
            FreX = u+u0;
            FreY = v+v0;

            % calculate grating period for X0 and Y0 direction
            FreX0 = (FreX-FreY)/2;
            FreY0 = (FreX+FreY)/2;

            % tilted phase 
            X0 = col*FreX0;
            Y0 = 2*row*FreY0;

%             XY = (X0+ Y0);

            % computer spherical wavefront
%             Wf = fi((row-(m+1)/2)*pSize,(col-(n+1)/2)*pSize,P+P0);

            % add titled phase
            theta = exp(1i*(X0 + Y0)*2*pi) + theta;
            theta_test(:,:,t) = (col2*FreX0 + row2*FreY0)*2*pi;
%             C(:,:,NF) = Wf + XY*2*pi;
 end
C1 = asin(abs(theta)/max(max(abs(theta))));
C2 = angle(theta);

fi1 = C1 + C2;
fi2 = C1 - C2 + pi;
threshold = 0;
M1 = cos(fi1);
M1(M1 > threshold) = 1;
M1(M1<= threshold) = 0;
M2 = cos(fi2);
M2(M2 > threshold) = 1;
M2(M2<= threshold) = 0;

M = zeros(n,2*m);
M(:,1:2:2*m) = M1;
M(:,2:2:2*m) = M2;

R = 1-M;


[x,y,F_DMD]=Cal_Fourier_Surface(R,px_size,4,Mobj,wl,fl);
% toc
figure(2);
imagesc(x,y,abs(F_DMD));
axis square;
axis xy;
colorbar;
colorbar;
colormap gray;
xlim([-100e-6 100e-6]);
ylim([-100e-6 100e-6]);


%%
 for i=1:1:FocusNum
    V(i) = abs(sum(sum(R.*exp(1i*(-theta_test(:,:,i)))))/(2*m*n));
 end
 V = V.^2;
 Uniformity=1-(max(abs(V))-min(abs(V)))/(max(abs(V))+min(abs(V)));
 I(t) = sum(V)/size(V,1);