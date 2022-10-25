%%%%%%%%%%%%%%%%%%%%
% This code generates a record matrxi covering all points in the sample
% plane and fill the scanning pattern to desired compression ratio.
% It has two outputs:
% 1. record matrix: sample_num * foci_num, 
% 2. hologram: both stacked and single frame is save to default folder

clear all; clc;

p = 100; q = 30;  %define the size of the image, p-y axis, q-x/z axis
stepsize  = 0.4;  %um
mindis = 2; % um, how to change this value and why to change???
                % --> maybe use this to control sparsity
z_flag = 1; % indicate this code to generate axial CS pattern or not
            % if z_flag = 1. swap p with third column in {Point} and add x
            % coordinate to first column in {Point}

% define compression ratio and scanning point
ratio = 0.3; 
m = round(ratio*p*q);

% define number of focus and store data
FocusNum = 4;
Point = cell(FocusNum,1);
Record = zeros(1,FocusNum);

% define kk as the record matrix, size = (num_sample, foci)
kk = randi(p*q,10*p*q,FocusNum);
% kk  = zeros(p*q, FocusNum);

% for i = 1:1:FocusNum
%   k = randperm(p*q);
%   k = k(:);
%   kk(:,i) = k;
% end

point_x = zeros(FocusNum,1);
point_y = zeros(FocusNum,1);
w = 1; % counter for valid CS pattern(multi-foci generation), should have (w-1)*foci>p*q

%% get about partial points that are not sampled twice
for i = 1:1:size(kk,1)
    for t = 1:1:FocusNum
         % calculate the col and row of designed linear index
         % in Matlab, y is column, x is row
         point_y(t) = ceil(kk(i,t)/q); % y is slow moving axis
         point_x(t) = kk(i,t) - (point_y(t) - 1)*q;
         % move corordinate to center
         point_x(t) = point_x(t) - ceil(q/2);
         point_y(t) = point_y(t) - ceil(p/2);
    end
    X = [point_x,point_y];
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
                Point{t}(w,3) = 0;
                Record(w,t) = kk(i,t);
            end
            w = w + 1; % this operation updates {Point} during every loop
        end
    end
end

%% find the point that is not sampled and re-generate
for cc = 1:1:20 % adjust the loop numner to make sure the effect
    if p*q-(w-1)*4 > 0
        % calculate number of points never been scanned
        Nonsampledarray = zeros(p*q-(w-1)*4,1);
        t = 1; % counter for points never been scanned
        for i = 1:1:p*q
          if isempty(find(Record==i))
              Nonsampledarray(t) = i;
              t = t+1;
          end
        end

        kk = randi(length(Nonsampledarray),4*p*q,4); % a new array to store point information
        for i = 1:1:size(kk,1)
            for t = 1:1:FocusNum
                 point_y(t) = ceil(Nonsampledarray(kk(i,t))/q);
                 point_x(t) = Nonsampledarray(kk(i,t)) - (point_y(t) - 1)*q;
                 point_x(t) = point_x(t) - ceil(q/2);
                 point_y(t) = point_y(t) - ceil(p/2);
            end
            X = [point_x,point_y];
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
                     Point{t}(w,3) = 0;
                     Record(w,t) = Nonsampledarray(kk(i,t));
                    end
                    w = w + 1;
                end
            end
        end
    end
end 

%% most pattern is fixed, but some are unsampled. Sampling these point with overlaped point
if p*q-(w-1)*4 > 0
    Nonsampledarray = zeros(p*q-(w-1)*4,1);
    t = 1;
    for i = 1:1:p*q
      if isempty(find(Record==i))
          Nonsampledarray(t) = i;
          t = t+1;
      end
    end

    kk = randi(p*q,p*q,4);
    cc = 1;
    for i = 1:1:size(kk,1)
        if cc > length(Nonsampledarray)
            break;
        end

        point_y(1) = ceil(Nonsampledarray(cc)/q);
        point_x(1) = Nonsampledarray(cc) - (point_y(1) - 1)*q;
        point_x(1) = point_x(1) - ceil(q/2);
        point_y(1) = point_y(1) - ceil(p/2);
        for t = 2:1:FocusNum
             point_y(t) = ceil(kk(i,t)/q);
             point_x(t) = kk(i,t) - (point_y(t) - 1)*q;
             point_x(t) = point_x(t) - ceil(q/2);
             point_y(t) = point_y(t) - ceil(p/2);
        end
        X = [point_x,point_y];
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
                 Point{1}(w,3) = 0;
                 Record(w,1) = Nonsampledarray(cc);
                for t = 2:1:FocusNum
                 Point{t}(w,1) = point_x(t);
                 Point{t}(w,2) = point_y(t);
                 Point{t}(w,3) = 0;
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
    kk = randi(p*q,10*p*q,4);
    for i = 1:1:size(kk,1)
        if w > m
            break;
        end
        for t = 1:1:FocusNum
             point_y(t) = ceil(kk(i,t)/q);
             point_x(t) = kk(i,t) - (point_y(t) - 1)*q;
             point_x(t) = point_x(t) - ceil(q/2);
             point_y(t) = point_y(t) - ceil(p/2);
        end
        X = [point_x,point_y];
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
                 Point{t}(w,3) = 0;
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
       k(:,t) = (Point{t}(:,1) +  ceil(q/2)) + (Point{t}(:,2) - 1 + ceil(p/2))*q;
 end

debug = Record - k; % good generation should have all-zero matrix

%% plot one of the focus for a single focus
point = Point{1}(:,1:2);
t = zeros(p,q);
for k = 1:1:m
    t(point(k,2)+ceil(p/2), point(k,1)+ ceil(q/2))=1;
end
imshow(t,[]);

%% scale the matrix using stepsize, consider axial generation
if z_flag == 1
    % add q/2 to third column make sure z > 0
    % swap first and third column in {Point}
    for i_cell = 1:1:size(Point,1)
        min_z = min(Point{i_cell}(:,1));
        Point{i_cell}(:,1) = Point{i_cell}(:,1) - min_z + 1;
        Point{i_cell}(:,[1 3]) = Point{i_cell}(:,[3 1]);
    end
end

%% save generated measurement matrix
% save data here will not change generated scanning points
save('RecordPointPair_X1Y100Z30_CS30Foci4.mat', 'Record','Point');

%% generate hologram
for t = 1:1:FocusNum
    Point{t} = stepsize * Point{t}(1:m,:);
%     Point{t} = 1 * Point{t}(1:m,:);
end
tic
% [Uniformity,I] =  TrajectoryBGS_Multi_40gpu(50, '.\experiment20221007\yzTest\', 'Particalscan_X1Y100Z30_CS60Foci2', Point{1:FocusNum});
toc
% save('I', 'I');
% save('Uniformity', 'Uniformity');