function mpi_drawEpi(img,sliceNum, studyNum, outpath, referenceFrame)

% to get initial set of ROIs:

% read in data
[nrows, ncols, nFrames]=size(img)

clf


figure(1)
clf
imagesc(img(:,:,1))
%spect=load('/home/mirl/ed/spect.cmap');
%m=max(max(spect));
%spect=spect/m;
%colormap(spect)

axis('square')
frameShown=referenceFrame
imagesc(img(:,:,frameShown))
axis('square')
% get some user input if this frame is ok to draw contours on...

% -- Added 04/03/2003 --  SV%
flag = input('Enter "1" if this frame is okay to draw contours else enter "0": ');
while (flag == 0)
        clf;
        frameShown = input('Enter which other frame you would like to try: ');
        imagesc(img(:,:,1))
%        axis('square')
        imagesc(img(:,:,frameShown))
%        axis('square')
        flag = input('Enter "1" if this frame is okay to draw contours else enter "0": ');
end;
%s = fix(clock);   % date;
disp('Ref. frame chosen being written to log file');
% -- To log the changes made -- %
if (flag == 1)
%    tmpname=strcat(outpath,'log.study',int2str(studyNum),'.slice',int2str(sliceNum),'.txt');
    f = fopen('log.txt','a+');
    fprintf(f,'Date: ');
    fprintf(f,' %d,\n',fix(clock));
    fprintf(f,'study: ');
    fprintf(f,'%d,\n',studyNum);
    fprintf(f,'slice: ');
    fprintf(f,'%d,\n',sliceNum);
    fprintf(f,'frameUsedforContours = frameShown =');
    fprintf(f,'%d,\n\n',frameShown);
end;


% display endo already chosen:

tmpname=strcat(outpath,'endo_polyCoords.study',int2str(studyNum),'.slice',int2str(sliceNum));
% to get back do e.g. 
ss=load(tmpname);
x1poly=ss(:,1); y1poly=ss(:,2);
bw1=mpi_roipoly(img(:,:,17),x1poly,y1poly);
hold
%fill(x1poly,y1poly,'r');
h1=line(x1poly,y1poly);
set(h1,'Color',[0 1 0]); 

[I, J]=find(bw1);

center1=[mean(x1poly) mean(y1poly)];   % is this a reasonable way to pick? probably if enough points chosen

[bw2,x2poly,y2poly]=mpi_roipoly;   %  polygon vertices
polyCoords=[x2poly y2poly];
tmpname=strcat(outpath,'epi_polyCoords.study',int2str(studyNum),'.slice',int2str(sliceNum));
save(tmpname, 'polyCoords','-ascii');
%save epi_polyCoords polyCoords -ascii;
hold
h2=line(x2poly,y2poly);
set(h2,'Color',[1 1 0]); 
center2=[mean(x2poly) mean(y2poly)];

center=(center1+center2)/2;
xCenter = center(1);
yCenter = center(2);
% now that have center write out blood ROI (move edges in a little)
for i=1:length(I)
    dist(i)=sqrt((J(i)-center(1))*(J(i)-center(1))+(I(i)-center(2)*(I(i)-center(2))));
end
maxdist=max(dist);
%fid=fopen(strcat(outpath,'bld.slice',int2str(sliceNum),'.dat'),'w');
%for i=1:length(I)
%    % check distance from center - move in 40% to just get blood
%    if( sqrt((J(i)-center(1))*(J(i)-center(1))+(I(i)-center(2)*(I(i)-center(2)))) < 0.6*maxdist)
%    	fprintf(fid,'%d	%d	%d \n',J(i),I(i),1);
%    end
%end
%fclose(fid);


% added from Dmitri 6/22/02
disp('Choose the starting direction by hiting right mouse button')
[x,y] = getline(gcf);

xStart = x(end);  
yStart = y(end);

dx   = xStart - xCenter;
dy   = -(yStart - yCenter);

h3 = line([xStart, xCenter],[yStart, yCenter],'Color', 'm');

start_angle = atan2(dy, dx)*180/pi + 180 - 135;
tmpname = strcat(outpath,'Roi_start_angle.study',int2str(studyNum),'.slice',int2str(sliceNum));
save(tmpname, 'start_angle','-ascii');

disp('Use  mouse to choose Blood contour. Backspace or delete to redo last polygon vertex point')
[bw3, x3, y3] = mpi_roipoly;   %  polygon vertices
polyCoords = [x3 y3];
tmpname = strcat(outpath,'blood_polyCoords.study',int2str(studyNum),'.slice',int2str(sliceNum));
save(tmpname, 'polyCoords','-ascii');
hold
line(x3,y3,'Color', [0.3 0.8 0.5]); 
%---blood curve------------------- 
N = sum(sum(bw3));
for j = 1 : nFrames
   a = img(:, :, j) .* bw3;
   bldcurve(j) = sum(sum(a)) / N;
end
%------------------------------------------
%figure(3);
%plot(times, bldcurve, '-', 'Color', s(i).col);
