function mpi_pickStartAngle(img,sliceNum, studyNum, outpath, referenceFrame)

% to get initial set of ROIs:

% read in data
[nrows, ncols, nFrames]=size(img)


% display endo already chosen:
tmpname=strcat(outpath,'endo_polyCoords.study',int2str(studyNum),'.slice',int2str(sliceNum));
% to get back do e.g. 
ss=load(tmpname);
x1poly=ss(:,1); y1poly=ss(:,2);
center1=[mean(x1poly) mean(y1poly)];   % is this a reasonable way to pick? probably if enough points chosen

%polyCoords=[x2poly y2poly];
tmpname=strcat(outpath,'epi_polyCoords.study',int2str(studyNum),'.slice',int2str(sliceNum));
ss2=load(tmpname);
x2poly=ss2(:,1); y2poly=ss2(:,2);
center2=[mean(x2poly) mean(y2poly)];

center=(center1+center2)/2;
xCenter = center(1);
yCenter = center(2);

if ~exist('referenceFrame')
 % then use image and graphic already displayed
else
   figure(1); clf;
   frameShown=referenceFrame;
   imagesc(img(:,:,frameShown))
   axis('image')
   hold
%fill(x1poly,y1poly,'r');
   h1=line(x1poly,y1poly);
   set(h1,'Color',[0 1 0]); 
   h2=line(x2poly,y2poly);
   set(h2,'Color',[1 1 0]); 
end


disp('...')
disp('Choose the starting direction for circumferential segmentation of the LV by clicking the right mouse button once the cursor is placed near the RV-LV insertion site')
disp('...')
disp('A line will be drawn from this point to the centroid of the LV, to then divide the LV myocardium into discrete equi-angular regions')

[x,y] = getline(gcf);

xStart = x(end);  
yStart = y(end);

dx   = xStart - xCenter;
dy   = -(yStart - yCenter);

h3 = line([xStart, xCenter],[yStart, yCenter],'Color', 'm');

start_angle = atan2(dy, dx)*180/pi + 180;
start_angle = atan2(dy, dx)*180/pi + 180 - 135 ;  % 11/29/04
start_angle = atan2(dy, dx)*180/pi + 180 + (360-120) ;  % 3/15/05
% for 6 regions 17 seg. model
% need to worry about how to identify angle in apical slice (4 regions)
%keyboard
% worry about if quad. II or III for sign?
tmpname = strcat(outpath,'Roi_start_angle.study',int2str(studyNum),'.slice',int2str(sliceNum));
save(tmpname, 'start_angle','-ascii');
