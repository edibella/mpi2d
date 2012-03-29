function [epi, endo, LV,Roi_start_angle] = Best_Segmentation(cinemri1,show,RV,LVpoint)
if(~exist('show'))
    show = 0;
end
[sx sy st] = size(cinemri1);
myvals = zeros(sx*sy,st);
side = 5;
mymean = mean(cinemri1(:,:,1:5),3);
mycroppedcine1 = zeros(sx,sy,st);
for t=1:st
    mycroppedcine1(:,:,t) = cinemri1(:,:,t) - mymean;
end
%smooth through time
padded = cat(3,repmat(mean(mycroppedcine1(:,:,1:side),3),[1 1 side]), mycroppedcine1);
h = fspecial('gaussian',[1 10],2);
temp = zeros(sx,sy,st);
for i=1:sx
    for j=1:sy
        curve = squeeze(padded(i,j,1:end));
        curve = filter2(h,curve);
        myvals((j-1)*sx+i,:) = curve((1+side):end);
        temp(i,j,:) = curve((1+side):end);
        %convert to deltaSI
        myvals((j-1)*sx+i,:) = myvals((j-1)*sx+i,:) - mean(myvals((j-1)*sx+i,1:5));
    end
end
if(~exist('RV'))
    [RV, LVpoint] = FindLVRV(cinemri1);
end
p1 = squeeze(padded(RV(1),RV(2),(1+side):end))';
p2 = squeeze(padded(LVpoint(1),LVpoint(2),(1+side):end))';
temp = round(mean(vertcat(RV,LVpoint),1));
p3 = squeeze(padded(temp(1),temp(2),(1+side):end))';
p4 = squeeze(padded(1,1,(1+side):end))';

IDX = kmeans(myvals,4,'distance','correlation','start',vertcat(p1,p2,p3,p4),'emptyaction','singleton');
master = reshape(IDX,[sx sy]);
if(show)
    figure(1), 
    imagesc(master), colorbar, title('Master regions');
end

se = strel('diamond',1);
RVvalue = master(RV(1), RV(2));
RVMask = imfill((master == RVvalue),RV);
LVvalue = master(LVpoint(1), LVpoint(2));
LVMask = imdilate(imfill(imerode(master==LVvalue,se),LVpoint),se);
LVMask = ~imfill(~LVMask,[1,1]);

if(show)
    figure(2), imagesc(LVMask); title('LV mask');
end
se = strel('diamond',3);
myofinder = imdilate(LVMask,se) - LVMask;
[r, c] = find(myofinder);
clear myoValue
myoValue = [];
myoPoint = [r(1) c(1)];
for i=1:length(r)
    myoValue = [myoValue master(r(i), c(i))];
end
myoValue = median(myoValue);
myoMask = (master == myoValue);


if(show)
    figure(3), imagesc(myoMask);title('initial myocardial mask');
end
myperim = bwperim(LVMask);
[r,c] = find(myperim);
[X, Y] = meshgrid(1:sx, 1:sy);
temp = zeros(sx,sy,1);
for i=1:length(r)
    mymask = sqrt((X-c(i)).^2 + (Y-r(i)).^2);
    mymask = 1./(1+(mymask/10).^2);
    mymask = reshape(mymask,sx,sy,1);
    temp = max(cat(3,temp,mymask),[],3);
end
myoMask = filter2(fspecial('gaussian',[sx,sy],1),myoMask);
newimg  = myoMask.* temp;
if(show)
    figure(4), imagesc(newimg)
end
goodRegion = newimg > .5;
myoMask = (myoMask .* goodRegion)>0;

%kill islands but preserve true mycardial islands
StuffToDo = myofinder.* myoMask;
runningMyo = zeros(sx,sy);
while(any(StuffToDo(:)>0))
    [r,c] = find(StuffToDo>0);
    myo = imfill(myoMask,[r(1), c(1)]);
    runningMyo = runningMyo + myo;
    StuffToDo = StuffToDo - StuffToDo.* myo;
end
myoMask = runningMyo;
myoMask = imdilate(myoMask,se);
myoMask = imerode(myoMask,se);

if(show)
    figure(6), imagesc(myoMask),title('final myocardial mask');
end

%Find which piece is the true myocardium
%find holes in the LV
temp = imfill(~LVMask,[1 1]);
LVMask = ~(temp);
myoMask = myoMask - (LVMask .* myoMask);
myoMask = imfill((myoMask + LVMask) > 0,LVpoint);
myoMask = myoMask - LVMask;

    

se = strel('diamond',8);
newLV = imerode(LVMask,se);
[r,c] = find(newLV);
LV = bwtraceboundary(newLV,[r(1) c(1)],'W',8);
se = strel('diamond',1);
epimask = imdilate(LVMask,se);

[r,c] = find(epimask);
epi = bwtraceboundary(epimask,[r(1) c(1)],'W',8);

temp = myoMask+imdilate(epimask,se);
x = find(temp(:,LVpoint(2)));
endo = bwtraceboundary(temp,[x(1),LVpoint(2)],'W',8);

%reconstruct myoMask based on epi and endo
outter = roipoly(myoMask,endo(:,2), endo(:,1));
inner = roipoly(myoMask,epi(:,2), epi(:,1));
mynewmyo = (outter - inner)>0;
[r, c] = find(mynewmyo);

myRGB = zeros(sx,sy,3);
temp = mean(cinemri1,3);
temp = temp - min(temp(:));
temp = temp / max(temp(:));
myRGB(:,:,1) = temp;
myRGB(:,:,2) = myRGB(:,:,1);
myRGB(:,:,3) = myRGB(:,:,1);
[r,c] = find(mynewmyo);
for i=1:length(r)
    myRGB(r(i),c(i),1) = min(1,myRGB(r(i),c(i),1)*2);
end

figure, imagesc(myRGB);title('Segmentation');


if(any(RVMask == LVMask))
    %the FindRVLV function couldn't find the RV so we'll hard code an
    %arbitrary angle
    Roi_start_angle = 45;
else
    
    [x,y] = find(RVMask);
    [theta,r] = cart2pol(x-LVpoint(1), y-LVpoint(2));
    [theta,IX] = sort(theta);
    r = r(IX);
    [apexPoint(1,:),apexPoint(2,:)] = pol2cart(theta,r);
    apexPoint = apexPoint(:,1)' + LVpoint;
    if(show)
        figure(5), hold on, plot(apexPoint(2), apexPoint(1), 'o');
    end
    Roi_start_angle = theta(1)/pi*360;
end
