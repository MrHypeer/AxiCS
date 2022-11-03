function pattern = TrajectoryByPoint_withcor_stack40Scan(Point,Amp,Cali,Path_File,suffixChar) 
%%
% 'Points' is an array with a dimension of Num*3. All the points need to be
% fabricated are contained in the 'Points' array.

%% Initial parameters
% Amp = gpuArray(Amp);
%radius of light spot on DMD
wf_r = 4*10^-3;

% wavelength, nm
wl = 780*10^-9; 

% spherical wavefront, for beam shaping, you can replace it by other functions
% P, optical power, equals to 1/f, unit m^-1
fi = @(x,y,P) pi*(x.^2+y.^2).*(P)./wl;

% control the width of the fringes
q = 1;

% pixels of the DMD
m = 1024; n = 768;

% size of a pixel, um
pSize = 13.68*10^-6;

[row,col] = meshgrid(1:m,1:n);
% original point where x=0, y=0, z=0
% u, v, P are x and y diretion's spatial frequency of DMD pattern, and  optical power, respectively
u0 = 0;
v0 = 0;
P0 = 0;

%% add zernike coefficient and generate stack bmp file

[Num, Useless] = size(Point);

disp('---PICTURE WRITING START---'); 

stack_num = 40;

bigPic = zeros(n * stack_num, m/8, 'uint8');

pic_cnt_total = 0;
pic_cnt_pack = 0;
onepackflag = (Num<stack_num);
pack_num = 1;
last_pack_num = mod(Num, stack_num);
total_pack = ceil(Num / stack_num);
for ii_sample = 1:Num
%     disp([num2str(ii_sample),'th sample is running...']);
    K=zeros(n,m);
    % a circular aperture
    for i=1:m
        for j=1:n
           if ((i-((m+1)/2))*pSize)^2+((j-((n+1)/2))*pSize)^2 <= wf_r^2
               K(j,i)=1;
           end
        end
    end
    
    %Zernike mode
    [theta,r] = cart2pol((row-(m+1)/2)*pSize/wf_r,(col-(n+1)/2)*pSize/wf_r);
    Zer = [2,-2;2,2;3,-3;3,-1;3,1;3,3;4,-4;4,-2;4,0;4,2;4,4;5,-5;5,-3;5,-1;5,1;5,3;5,5];%17
    cor = zeros(n,m);
    Amp_sig = Amp;
    for i = 1:size(Amp_sig')
        if Amp_sig(i) ~= 0
               %cor = cor + Amp(i)*zernfun_GPU(Zer(i,1),Zer(i,2),r,theta,'norm')+Cali(i)*Amp(i)*zernfun_GPU(2,0,r,theta,'norm');
            cor = cor + Amp_sig(i)*zernfun(Zer(i,1),Zer(i,2),r,theta,'norm')+Cali(i)*Amp_sig(i)*zernfun(2,0,r,theta,'norm');
        end
    end

% z = Points(1,3);
% P = z/32.50;
% A = fi((row-(m+1)/2)*pSize,(col-(n+1)/2)*pSize,P+P0);

        x = Point(ii_sample,1);
        y = Point(ii_sample,2);
        z = Point(ii_sample,3);
% accordingly parameters, COCHE TPM

         u = x/109;
         v = y/109;
         P = z/9.7370;



% desired grating spatial frequency for X and Y direction
        FreX = u+u0;
        FreY = v+v0;

% calculate grating period for X0 and Y0 direction
        FreX0 = (FreX-FreY)/2;
        FreY0 = (FreX+FreY)/2;

% tilted phase 
        X0 = col*FreX0;
        Y0 = row*FreY0;

        XY = (X0+ Y0);

% computer spherical wavefront
       A = fi((row-(m+1)/2)*pSize,(col-(n+1)/2)*pSize,P+P0);
%         A = 0;
  
    
% add titled phase
        C = A./(2*pi)+cor+ XY;%+RNN_phase.RNN_phase(n*(RNN_idx-1)+1:n*RNN_idx,:);
%         C= A+XY;
        M = C-floor(C);
        Mf = abs(M);

% according to the Lee Holography
        Mf(Mf < q/2) = 0;
        Mf(Mf >= q/2) = 1;
        
        Mf = Mf.*K;
        
        R = 1-Mf;

       
%         R=R.*K;
        count_large = 0;
        if count_large < 10
            PatternName = [Path_File,'X',num2str(x),'_Y',num2str(y),'_Z',num2str(z),'.bmp'];
%             imwrite(logical(R),PatternName); % write large image
            count_large = count_large+1;
        end

% transfer to DLP4100 required mode
        stride=m/8;

%        im = gpuArray(uint8(zeros(n,stride)));
        im = uint8(zeros(n,stride));
        for ij=1:stride
            im(:,ij)=128*R(:,8*ij-7)+64*R(:,8*ij-6)+32*R(:,8*ij-5)+16*R(:,8*ij-4)+8*R(:,8*ij-3)+4*R(:,8*ij-2)+2*R(:,8*ij-1)+R(:,8*ij);
        end
         
% save each binary patterns as bmp file, each file is around 96 KB
        pic_cnt_total = pic_cnt_total + 1;
        pic_cnt_pack = pic_cnt_pack + 1;
        bigPic(((pic_cnt_pack - 1) * n + 1) : pic_cnt_pack * n, :) = im;
        if (onepackflag == 1) && (pic_cnt_total == Num)
        % writing to bmp file
            imwrite(bigPic, strcat(Path_File,num2str(pack_num,'%4d'),suffixChar, '.bmp'));


            disp(['The ', num2str(Num), ' image has been finished, with ', num2str(total_pack - pack_num), ' packs remained.']);
            bigPic = zeros(n * stack_num, m/8, 'uint8');
        end
        if pic_cnt_pack == stack_num
            % writing to bmp file
            imwrite(bigPic, strcat(Path_File,num2str(pack_num,'%4d'),suffixChar,'.bmp'));
            disp(['The ', num2str(pack_num), ' pack has been finished, with ', num2str(total_pack - pack_num), ' packs remained.']);
            pic_cnt_pack = 0;
            pack_num = pack_num + 1;
            bigPic = zeros(n * stack_num, m/8, 'uint8');
        end
end
pattern = 'ok';
bigPic = bigPic(1 : (last_pack_num * n), :);
imwrite(bigPic, strcat(Path_File, num2str(pack_num), suffixChar, '.bmp'));

%% Mr.Geng's code
%% Modified by Mindan