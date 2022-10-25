%%%%%%%%%%%%%%%%%%%%
% This code generates a record matrxi covering all points in the sample
% plane and fill the scanning pattern to desired compression ratio.
% It has two outputs:
% 1. record matrix: sample_num * foci_num, 
% 2. hologram: both stacked and single frame is save to default folder
%
% change the code, to make sure the z-position is right for hologram
% generation. But the sampling distance is the same.@20221007 LIU Gan
% 
% change sampling distance and allow anisotropy scaling. @20221011 LIU Gan
% 

clear all; clc;

%% define compression ratio and scanning point

p = 100; q = 1; z = 30;  %define the size of the image, p-y axis, q-x, z-z
stepsize  = 0.4;  %um
stepsize_zscale = 2.5;
mindis = 5; % um, how to change this value and why to change???
                % --> maybe use this to control sparsity
ratio = 0.6; FocusNum = 2;
record_name = 'RecordPointPair_X1Y100Z30_CS60Foci2Dis5.mat';
suffix = 'Particalscan_X1Y100Z30_CS60Foci2Dis5';
pattern_path = '.\experiment20221007\yzTest2\';

m = round(ratio*p*q*z);
Point = cell(FocusNum,1);
Record = zeros(1,FocusNum);

% define kk as the record matrix, size = (num_sample, foci)
kk = randi(p*q*z,10*p*q*z,FocusNum);
% kk  = zeros(p*q, FocusNum);

% for i = 1:1:FocusNum
%   k = randperm(p*q);
%   k = k(:);
%   kk(:,i) = k;
% end

point_x = zeros(FocusNum,1);
point_y = zeros(FocusNum,1);
point_z = zeros(FocusNum,1);
w = 1; % counter for valid CS pattern(multi-foci generation), should have (w-1)*foci>p*q

%% get about partial points that are not sampled twice
for i = 1:1:size(kk,1)
    for t = 1:1:FocusNum
         % calculate the col and row of designed linear index
         % in Matlab, y is column, x is row
         point_z(t) = floor((kk(i,t)-1)/p/q);
         point_y(t) = ceil((kk(i,t)-point_z(t)*p*q)/q); % y is slow moving axis
         point_x(t) = kk(i,t) - point_z(t)*p*q - (point_y(t)-1)*q;
         % move corordinate to center
         point_x(t) = point_x(t) - ceil(q/2);
         point_y(t) = point_y(t) - ceil(p/2);
    end
    % use stepsize_zscale to adjust distance, because X is only used for
    % distance judgement, so dont need to change value in {Point}
    X = [point_x,point_y,point_z*stepsize_zscale]; 
    % calculate pairwise distance
    distmat = pdist(X);
    
    % do following operation only when distance > 20 pxl
    if stepsize*min(distmat)> mindis 
        overlap = 0;  % indicator whether CS scanning contain duplicate points
        for t = 1:1:FocusNum
            % update overlap only when (~isempty(find(Record==kk(i,t)))) = 1
            % one or more than one points can flip overlap
            overlap = ~isempty(find(Record==kk(i,t))) || overlap; 
        end
        % store point information into effectve CS pattern
        if overlap == 0
            for t = 1:1:FocusNum
                Point{t}(w,1) = point_x(t);
                Point{t}(w,2) = point_y(t);
                Point{t}(w,3) = point_z(t);
                Record(w,t) = kk(i,t);
            end
            w = w + 1; % this operation updates {Point} during every loop
        end
    end
end

%% find the point that is not sampled and re-generate
for cc = 1:1:20 % adjust the loop numner to make sure the effect
    if p*q*z-(w-1)*4 > 0
        % calculate number of points never been scanned
        Nonsampledarray = zeros(p*q*z-(w-1)*4,1);
        t = 1; % counter for points never been scanned
        for i = 1:1:p*q*z
          if isempty(find(Record==i))
              Nonsampledarray(t) = i;
              t = t+1;
          end
        end

        kk = randi(length(Nonsampledarray),4*p*q*z,FocusNum); % a new array to store point information
        for i = 1:1:size(kk,1)
            for t = 1:1:FocusNum
                point_z(t) = floor((Nonsampledarray(kk(i,t))-1)/p/q);
                point_y(t) = ceil((Nonsampledarray(kk(i,t))-p*q*point_z(t))/q);
                point_x(t) = Nonsampledarray(kk(i,t)) - p*q*point_z(t)-(point_y(t) - 1)*q;
                point_x(t) = point_x(t) - ceil(q/2);
                point_y(t) = point_y(t) - ceil(p/2);
            end
            X = [point_x,point_y, point_z*stepsize_zscale];
            distmat = pdist(X);
            if stepsize*min(distmat)> mindis
                overlap = 0;
                for t = 1:1:FocusNum
                    overlap = ~isempty(find(Record==Nonsampledarray(kk(i,t)))) || overlap;
                end

                if overlap == 0
                    for t = 1:1:FocusNum
                    	Point{t}(w,1) = point_x(t);
                        Point{t}(w,2) = point_y(t);
                        Point{t}(w,3) = point_z(t);
                        Record(w,t) = Nonsampledarray(kk(i,t));
                    end
                    w = w + 1;
                end
            end
        end
    end
end 

%% most pattern is fixed, but some are unsampled. Sampling these point with overlaped point
if p*q*z-(w-1)*4 > 0
    Nonsampledarray = zeros(p*q*z-(w-1)*4,1);
    t = 1;
    for i = 1:1:p*q*z
      if isempty(find(Record==i))
          Nonsampledarray(t) = i;
          t = t+1;
      end
    end

    kk = randi(p*q*z,p*q*z,FocusNum);
    cc = 1;
    for i = 1:1:size(kk,1)
        if cc > length(Nonsampledarray)
            break;
        end

%         point_y(1) = ceil(Nonsampledarray(cc)/q);
%         point_x(1) = Nonsampledarray(cc) - (point_y(1) - 1)*q;
        
        point_z(1) = floor((Nonsampledarray(cc)-1)/p/q);
        point_y(1) = ceil((Nonsampledarray(cc)-p*q*point_z(1))/q);
        point_x(1) = Nonsampledarray(cc) - p*q*point_z(1)-(point_y(1) - 1)*q;
        point_x(1) = point_x(1) - ceil(q/2);
        point_y(1) = point_y(1) - ceil(p/2);
        for t = 2:1:FocusNum
            point_z(t) = floor((kk(i,t)-1)/p/q);
            point_y(t) = ceil((kk(i,t)-p*q*point_z(t))/q);
            point_x(t) = kk(i,t) - p*q*point_z(t)-(point_y(t) - 1)*q;
            point_x(t) = point_x(t) - ceil(q/2);
            point_y(t) = point_y(t) - ceil(p/2);
        end
        X = [point_x,point_y, point_z*stepsize_zscale];
        distmat = pdist(X);
        if stepsize*min(distmat)> mindis
            overlap = length(find(Record==Nonsampledarray(cc)))>1;
            for t = 2:1:FocusNum
                overlap = length(find(Record==kk(i,t)))>1 || overlap;
            end
    %         overlap
            if overlap == 0
                Point{1}(w,1) = point_x(1);
                Point{1}(w,2) = point_y(1);
                Point{1}(w,3) = point_z(1);
                Record(w,1) = Nonsampledarray(cc);
                for t = 2:1:FocusNum
                    Point{t}(w,1) = point_x(t);
                    Point{t}(w,2) = point_y(t);
                    Point{t}(w,3) = point_z(t);
                    Record(w,t) = kk(i,t);
                end
                w = w + 1;
                cc = cc + 1 ;
            end
        end
    end
end
%% all the point are sampled, then generate more point to the satisify the undersampling ratio
% all points are scanned, but the CS ratio should be also met
for ii = 1:3
    kk = randi(p*q*z,10*p*q*z,FocusNum);
    for i = 1:1:size(kk,1)
        if w > m
            break;
        end
        for t = 1:1:FocusNum
            point_z(t) = floor((kk(i,t)-1)/p/q);
            point_y(t) = ceil((kk(i,t)-p*q*point_z(t))/q);
            point_x(t) = kk(i,t) - p*q*point_z(t)-(point_y(t) - 1)*q;
            point_x(t) = point_x(t) - ceil(q/2);
            point_y(t) = point_y(t) - ceil(p/2);
        end
        X = [point_x,point_y, point_z*stepsize_zscale];
        distmat = pdist(X);
        if stepsize*min(distmat)> mindis
            overlap = 0;
            for t = 1:1:FocusNum
                overlap = length(find(Record==kk(i,t)))>1 || overlap;
            end
    %         overlap
            if overlap == 0
                for t = 1:1:FocusNum
                 Point{t}(w,1) = point_x(t);
                 Point{t}(w,2) = point_y(t);
                 Point{t}(w,3) = point_z(t);
                 Record(w,t) = kk(i,t);
                end
                w = w + 1;
            end
        end
    end
end
for t = 1:1:FocusNum
    Point{t} = Point{t}(1:m,:);
end
Record = Record(1:m,:);

%%

%  t = zeros(p,q);
%  t1 = zeros(p,q);
%  cy  = 10;
% for k = 1:1:FocusNum
%     t(Point{k}(cy,2)+ceil(q/2),Point{k}(cy,1)+ ceil(p/2))=1;
%     
%     point_y = ceil(Record(cy,k)/q);
%     point_x = Record(cy,k) - (point_y - 1)*q ;
%     t1(point_y,point_x) = 1;
% end
% 
% subplot(2,1,1);imshow(t,[]);
% subplot(2,1,2);imshow(t1,[]);
% 
% 
% point_y = ceil(Record/q);
% point_x = Record - (point_y - 1)
%% check generated matrix is right or not
k = zeros(size(Record));
 for t = 1:1:FocusNum
       k(:,t) = (Point{t}(:,1) +  ceil(q/2)) + (Point{t}(:,2) - 1 + ceil(p/2))*q ...
           + Point{t}(:,3)*p*q;
 end

debug = Record - k; % good generation should have all-zero matrix

%% plot one of the focus for a single focus
point = Point{1}(:,1:2);
t = zeros(p,q);
for k = 1:1:m
    t(point(k,2)+ceil(p/2), point(k,1)+ ceil(q/2))=1;
end
imshow(t,[]);


%% save generated measurement matrix
% save data here will not change generated scanning points
save(record_name, 'Record','Point');

%% generate hologram
for t = 1:1:FocusNum
    Point{t}(1:m,1:2) = stepsize * Point{t}(1:m,1:2);
    Point{t}(1:m,3) = stepsize * stepsize_zscale * Point{t}(1:m,3);
%     Point{t} = 1 * Point{t}(1:m,:);
end
tic
% [Uniformity,I] =  TrajectoryBGS_Multi_40gpu(50, pattern_path, suffix, Point{1:FocusNum});
toc
% save('I', 'I');
% save('Uniformity', 'Uniformity');