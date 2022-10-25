function [Uniformity,I] = TrajectoryBGS_Multi_40gpu(iterations, Path_File, suffixChar, varargin)
%%%%%%
% 'Points' is an array with a dimension of Num*3. All the points need to be
% 'Path_File' is a path to write image,with '/' in the end, e.g. Path_File = 'D:\DLP_figures\Dien2\20170817\A\';
% fabricated are contained in the varargin
% 
% change paramters fitting into COCHE TPM. @20220929 LIU Gan


%% Initial parameters

% wavelength, nm
wl = 780*10^-9; 

% spherical wavefront, for beam shaping, you can replace it by other functions
% P, optical power, equals to 1/f, unit m^-1
fi = @(x,y,P) pi*(x.^2+y.^2).*(-P)./wl;

% control the width of the fringes
% q = 1;

% pixels of the DMD
m = 1024; n = 768;

% size of a pixel, um
pSize = 13.68*10^-6;

[row,col] = meshgrid(1:m,1:n);
row = gpuArray(row);
col = gpuArray(col);

% original point where x=0, y=0, z=0
% u, v, P are x and y diretion's spatial frequency of DMD pattern, and  optical power, respectively
% u0 = 0.1;
% v0 = 0;
% P0 = 0;

% on COCHE TPM
u0 = 0;
v0 = 0;
P0 = 0;
%% input trajactories
NumFocus = nargin - 3; % use number of input to calculate foci number
Points = varargin;
clear varargin;
NumPoints = zeros(NumFocus, 1); %Points Number of each trajactory
V = gpuArray(zeros(NumFocus,1));
for NF = 1:NumFocus
    [NumPoints(NF),~] = size(Points{NF});
end
Num = max(NumPoints);

%% scanning trajactories
% under the system in lab303, the x scan range is (0 um ,66 um), y scan range is (-82 um ,82 um),z scan range is (-30 um ,30 um)
% The unit is micron for the positions.


disp('---PICTURE WRITING START---'); 
stack_num = 40;
bigPic = zeros(768 * stack_num, 128, 'uint8');

pic_cnt_total = 0;
pic_cnt_pack = 0;
pack_num = 1;
last_pack_num = mod(Num, stack_num);
total_pack = ceil(Num / stack_num);


Uniformity = gpuArray(zeros(Num,1));
I  = gpuArray(zeros(Num,1));
C = gpuArray(zeros(n,m,NumFocus));

for t=1:Num
    N_effective = 0;
%     FocusDel = [];
    for NF=1:NumFocus
%         if t>NumPoints(NF)
%             FocusDel = [FocusDel,NF];
%         else
            N_effective = N_effective + 1;
            x = Points{NF}(t,1) + 30; % change this parameter according
            y = Points{NF}(t,2);
            z = Points{NF}(t,3); % uncomment this line when generating axial pattern
            
            % accordingly parameters
%             u = x/165.43;
%             v = y/165.43;
%             P = z/24.24;
% accordingly parameters under optical settng: 54 mm scan lens; 175 mm tube
% lens, 25X water immersion
%               u = x/199.05;
%               v = y/199.05;
%               P = z/32.50;

% For COCHE TPM
            u = x/109;
            v = y/109;
            P = z/9.7370;
            
            % desired grating spatial frequency for X and Y direction
            FreX = u+u0;
            FreY = v+v0;

            % calculate grating period for X0 and Y0 direction
            FreX0 = gpuArray((FreX-FreY)/2);
            FreY0 = gpuArray((FreX+FreY)/2);

            % tilted phase 
            X0 = col*FreX0;
            Y0 = row*FreY0;

            XY = (X0+ Y0);
            
            % add titled phase
%             C(:,:,NF) = (X0+ Y0)*2*pi; % calculate x-y only

            % computer spherical wavefront for z-depth and add tilt phase
            Wf = fi((row-(m+1)/2)*pSize,(col-(n+1)/2)*pSize,P+P0);
            C(:,:,NF) = Wf + XY*2*pi;
%         end
    end
    % according to the Lee Holography
%     Ef_Focus = 1:NumFocus;
%     Ef_Focus(FocusDel) = [];
	weights = gpuArray.ones(NumFocus,1);
    theta = gpuArray.zeros(NumFocus,1);
% 	theta = 2*pi*rand(length(Ef_Focus),1);
 
    num_it = 0;
	for iter=1:iterations
         SupField = gpuArray.zeros(n,m);
         for k=1:NumFocus
             SupField = SupField + weights(k)*exp(1i*(theta(k)+C(:,:,k)));
         end
        M = cos(angle(SupField));
        M(M > 0) = 1;
        M(M <= 0) = 0;
        R = 1 - M;
        for k=1:NumFocus
           V(k) = sum(sum(R.*exp(1i*(-C(:,:,k)))))/(m*n);
        end
        VV = abs(V);
        VV = VV.^2;
        I(t) = sum(VV)/size(V,1);
        Uniformity(t)=1-(max(abs(V))-min(abs(V)))/(max(abs(V))+min(abs(V)));
        if(Uniformity(t) > 0.99) % change threshold may can help converge but low fidelity
            break;
        else
            weights = weights*(sum(abs(V))/length(V))./abs(V);
            theta = angle(V);
        end
%       num_it = num_it+1
%       Uniformity(t)   
	end
      R = gather(R);
%     
%     if(mod(t,10)==0)
%         disp(['point ',num2str(t)]);
%     end

% transfer to DLP4100 required mode
       stride=128;
       im = uint8(zeros(n,stride));
        for i=1:stride
               im(:,i)=128*R(:,8*i-7)+64*R(:,8*i-6)+32*R(:,8*i-5)+16*R(:,8*i-4)+8*R(:,8*i-3)+4*R(:,8*i-2)+2*R(:,8*i-1)+R(:,8*i);
        end
    %% save single pic for EasyProj
    if t<21
        imwrite(R, strcat(Path_File, num2str(t), suffixChar, '_large.bmp')); 
    end

% save each binary patterns as bmp file, each file is around 96 KB
        pic_cnt_total = pic_cnt_total + 1;
        
        pic_cnt_pack = pic_cnt_pack + 1;
        
        bigPic(((pic_cnt_pack - 1) * 768 + 1) : pic_cnt_pack * 768, :) = im;
        
        if pic_cnt_pack == stack_num
            
            imwrite(bigPic, strcat(Path_File, num2str(pack_num), suffixChar, '.bmp')); % save stacked pic
            
            disp(['The ', num2str(pack_num), ' pack has been finished, with ', num2str(total_pack - pack_num), ' packs remained.']);
            
            pic_cnt_pack = 0;
            
            pack_num = pack_num + 1;
            
            bigPic = zeros(768 * stack_num, 128, 'uint8');
            
        end
        
        
%       imwrite(R,strcat(num2str(10000+t),'.bmp'));   
end      

bigPic((last_pack_num * 768 + 1) : ((last_pack_num + 1) * 768), :) = zeros(768, 128, 'uint8');

bigPic = bigPic(1 : ((last_pack_num + 1)* 768), :);

imwrite(bigPic, strcat(Path_File, num2str(pack_num), suffixChar, '.bmp'));

Uniformity = gather(Uniformity);
I = gather(Uniformity);
disp('---PICTURE WRITING FINISHED---');
end