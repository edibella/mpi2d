function [X,Y,xCenter, yCenter, start_angle,bw1,bw2,bw3]= mpi_loadContours(img,sliceNum, studyNum, outpath, endoFlag)  


if ~exist('endoFlag')
  endoFlag=0;
end

[nrows, ncols, nFrames]=size(img);
h = figure
clf
tmpImg=mean(img(:,:,:),3);  % hard-coded instead of passing in referenceFrame
imagesc(tmpImg), colormap gray;
hold on
axis('image')

% add error checking since did have different filename!!!
try
	tmp=load(strcat(outpath,'endo_polyCoords.study',int2str(studyNum),'.slice',int2str(sliceNum)));
catch
	disp('Dropped rest or stress part of filename March 2002. Will now attempt to read old filename')	
	tmp=load(strcat(outpath,'endo_polyCoords.rest.study',int2str(studyNum),'.slice',int2str(sliceNum)));
end
x1poly=tmp(:,1); y1poly=tmp(:,2);
bw1=mpi_roipoly(tmpImg,x1poly,y1poly);
%h1=line(x1poly,y1poly);
%set(h1,'Color',[0 1 0]); 
[I, J]=find(bw1);
center1=[mean(x1poly) mean(y1poly)];

try
	tmp=load(strcat(outpath,'epi_polyCoords.study',int2str(studyNum),'.slice',int2str(sliceNum)));
catch
	disp('Dropped rest or stress part of filename March 2002. Will now attempt to read old filename')	
	tmp=load(strcat(outpath,'epi_polyCoords.rest.study',int2str(studyNum),'.slice',int2str(sliceNum)));
end
x2poly=tmp(:,1); y2poly=tmp(:,2);
bw2=mpi_roipoly(tmpImg,x2poly,y2poly);
%h2=line(x2poly,y2poly);
%set(h2,'Color',[1 1 0]); 
center2=[mean(x2poly) mean(y2poly)];

center=(center1+center2)/2;
xCenter = center(1);
yCenter = center(2);

try
   tmp = load(strcat(outpath,'Roi_start_angle.study',int2str(studyNum),'.slice',int2str(sliceNum)));
catch
   disp('Dropped rest or stress part of filename March 2002. Will now attempt to read old filename')	
   tmp = load(strcat(outpath,'Roi_start_angle.rest.study',int2str(studyNum),'.slice',int2str(sliceNum)));
end
start_angle = tmp;

%-load blood Roi-
try
   tmp=load(strcat(outpath,'blood_polyCoords.study',int2str(studyNum),'.slice',int2str(sliceNum)));
catch
   disp('Dropped rest or stress part of filename March 2002. Will now attempt to read old filename')	
   tmp=load(strcat(outpath,'blood_polyCoords.rest.study',int2str(studyNum),'.slice',int2str(sliceNum)));
end
x3 = tmp(:,1); y3 = tmp(:,2);
bw3 = mpi_roipoly(tmpImg,x3,y3);
%line(x3, y3, 'Color', [0.5 0.3 0.9]);


myo = double(bw2)-double(bw1);   % check all are 0 or 1 ... (2nd mask encompasses first!)
[Y, X] = find(myo);



if(endoFlag~=0)
try
        tmp=load(strcat(outpath,'midwall_polyCoords.study',int2str(studyNum),'.slice',int2str(sliceNum)));
catch
        disp('Not finding midwall contour - may not have yet donw stage 3.25 to raw it')
        disp('DO NOT USE THIS DATA - CURVES WRITTEN HERE NOT VALID WITHOUT MIDWALL!')
        return;
end
x2poly=tmp(:,1); y2poly=tmp(:,2);
bw5=mpi_roipoly(tmpImg,x2poly,y2poly);
%h2=line(x2poly,y2poly);
%set(h2,'Color',[1 1 1]);
%center2=[mean(x2poly) mean(y2poly)];
%center=(center1+center3)/2;
clear myo, X, Y
if(endoFlag==1)
   myo = double(bw5)-double(bw1);
   [Y, X] = find(myo);
end
if(endoFlag==2)
   myo = double(bw2)-double(bw5); 
   [Y, X] = find(myo);
end

end  % if

% To replace with fits, make it seem myo is all image
%clear X
%clear Y
%[a,b]=size(tmpImg);
%new_img=zeros(a,b);
%new_img((1:a),(1:b))=1;
%[Y,X]=find(new_img);
%%keyboard
%%save ('Y.mat','Y')
%% save ('X.mat','X')
close(h);
return;
