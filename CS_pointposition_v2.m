clear all; clc;

%CS array



p = 100; q = 100;  %define the size of the image 
stepsize  = 0.5;  %um
mindis = 3;

ratio = 0.4;
m = round(ratio*p*q);

FocusNum = 20;
Point = cell(FocusNum,1);
Record = zeros(m,FocusNum);

kk = randi(p*q,10*p*q,FocusNum);

point_x = zeros(FocusNum,1);
point_y = zeros(FocusNum,1);
w = 1;
for i = 1:1:size(kk,1)
    for t = 1:1:FocusNum
         point_y(t) = ceil(kk(i,t)/q);
         point_x(t) = kk(i,t) - (point_y(t) - 1)*q;
         point_x(t) = point_x(t) - ceil(q/2);
         point_y(t) = point_y(t) - ceil(p/2);
    end
    X = [point_x,point_y];
    distmat = pdist(X);
    if stepsize*min(distmat)> mindis
        for t = 1:1:FocusNum
         Point{t}(w,1) = point_x(t);
         Point{t}(w,2) = point_y(t);
         Point{t}(w,3) = 0;
         Record(w,t) = kk(i,t);
        end
        w = w + 1;
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
% point_x = Record - (point_y - 1)*q ;
%%
    
save('Record', 'Record');

%%
k = zeros(size(Record));
 for t = 1:1:FocusNum
       k(:,t) = (Point{t}(:,1) +  ceil(q/2)) + (Point{t}(:,2) - 1 + ceil(p/2))*q;
 end

debug = Record - k;

%% plot one of the focus 
point = Point{3}(:,1:2);
t = zeros(p,q);
for k = 1:1:m
    t(point(k,1)+ ceil(p/2), point(k,2)+ceil(q/2))=1;
end
imshow(t,[]);

%%
for t = 1:1:FocusNum
    Point{t} = stepsize * Point{t}(1:m,:);
%     Point{t} = 1 * Point{t}(1:m,:);
end

tic
% [Uniformity,I] =  TrajectoryBGS_Multi_40gpu(15, 'F:\CS pattern\Pattern0408\', 'a', Point{1:FocusNum});
toc
save('I', 'I');
save('Uniformity', 'Uniformity');