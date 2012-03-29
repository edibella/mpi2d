function mpi_drawEndoAndEpi(img,sliceNum, studyNum, outpath, referenceFrame)

% to get initial set of ROIs:

% read in data
[nrows, ncols, nFrames]=size(img)

close all
% varimg = zeros(size(img,1),size(img,2));
% for i=1:size(img,1)
%     for j=1:size(img,2)
%         varimg(i,j) = var(img(i,j,:));
%     end
% end
% varimg = varimg/(max(varimg(:)))*max(img(:));
% rgbimg = zeros(size(img,1),size(img,2),3);
% rgbimg(:,:,1) = img(:,:,referenceFrame);
% rgbimg(:,:,2) = img(:,:,referenceFrame);
% rgbimg(:,:,3) = img(:,:,referenceFrame);
% q1 = mean(varimg(:));%+.5*std(varimg(:));
% for i=1:size(img,1)
%     for j=1:size(img,2)
%         if varimg(i,j)>.03
%             %rgbimg(i,j,1) = max(max(min(varimg(i,j)/q1,1),0),rgbimg(i,j,1));
%         end
%     end
% end
%rgbimg(varimg(:,:)>.03) = max(min(varimg(:,:)/q1,1),0);


figure(1)
clf
%spect=load('/home/mirl/ed/spect.cmap');
%m=max(max(spect));
%spect=spect/m;
%colormap(spect)

axis('image')
%imagesc(maxImg) % maxImg not very representative! big motion as bolus hit lv!
%disp('shoswing maxImg in case want to choose epi and endo from it. Hit key to go on ') 
frameShown=referenceFrame
%imagesc(img(:,:,frameShown))
imagesc(img(:,:,frameShown)),colormap gray
axis('image')
% get some user input if this frame is ok to draw contours on...

% -- Added 04/03/2003 --  SV%
chooseRefFrame = input('Enter "1" if this frame is okay to draw contours on else enter another frame to try ');
while (chooseRefFrame ~= 1)
        frameShown=chooseRefFrame;
        imagesc(img(:,:,chooseRefFrame)); axis('image')
        chooseRefFrame = input('Enter "1" if this frame is okay to draw contours on else enter another frame to try ');
end;
%s = fix(clock);   % date;
disp('Ref. frame chosen being written to log file');
% -- To log the changes made -- %
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

    
disp('...')
disp('Use the mouse to trace the boundary contour of the endocardium first (try to draw contour away from the LV blood pool), then press "Enter" or double click.')
set(gcf,'Renderer','zbuffer');   % to stop the blinking. Trying 3/22/05
%[bw1,x1poly,y1poly]=mpi_roipoly;   %  polygon vertices.See Notes 10/03, modified mpi_roipoly so you don't need to click mouse at each vertices.
[bw1,x1poly,y1poly]=mpi_roipoly;
polyCoords=[x1poly y1poly];
tmpname=strcat(outpath,'endo_polyCoords.study',int2str(studyNum),'.slice',int2str(sliceNum));
save(tmpname, 'polyCoords','-ascii');
% to get back do e.g. 
%ss=load(tmpname);
%x1poly=ss(:,1); y1poly=ss(:,2);
%bw1=mpi_roipoly(img(:,:,17),x1poly,y1poly);
hold
%fill(x1poly,y1poly,'r');
h1=line(x1poly,y1poly);
set(h1,'Color',[0 1 0]); 

[I, J]=find(bw1);

center1=[mean(x1poly) mean(y1poly)];   % is this a reasonable way to pick? probably if enough points chosen


disp('...')
disp('Use the mouse to trace the boundary contour of the epicardium (try to draw contour away from the pericardial fat layer), then press "Enter" or double click.')
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


% added  6/22/02

% Next to choose a starting angle for divididing the LV in equi-angular regions circumferentially 
mpi_pickStartAngle(img,sliceNum, studyNum, outpath)
%[x,y] = getline(gcf);
%xStart = x(end);  
%yStart = y(end);
%dx   = xStart - xCenter;
%dy   = -(yStart - yCenter);
%h3 = line([xStart, xCenter],[yStart, yCenter],'Color', 'm');
%start_angle = atan2(dy, dx)*180/pi + 180 - 135;
%tmpname = strcat(outpath,'Roi_start_angle.study',int2str(studyNum),'.slice',int2str(sliceNum));
%save(tmpname, 'start_angle','-ascii');

disp('...')
disp('Use the mouse to draw a blood contour within the LV blood pool to generate an arterial input function (AIF) for this slice')
[bw3, x3, y3] = mpi_roipoly;   %  polygon vertices

%has it put an rbg picture into a 2d matrix
if(size(bw3,2)/size(img,2) == 3)
    bw3 = reshape(bw3,[size(img,1) size(img,2) 3]);
    bw3 = (bw3(:,:,1) + bw3(:,:,2)+bw3(:,:,3))>0;
end
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
