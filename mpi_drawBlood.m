function mpi_drawBlood(img, sliceNum, studyNum, outpath)

% read in data 
[nrows, ncols, nFrames]=size(img)
figure(1)
tmpImg=img(:,:,15);
imagesc(tmpImg);
colormap('gray')
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
h1=line(x1poly,y1poly);
set(h1,'Color',[0 1 0]); 
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
h2=line(x2poly,y2poly);
set(h2,'Color',[1 1 0]); 
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


chooseRefFrame = input('Enter "1" if this frame is okay to draw contours on elseenter another frame to try ');
while (chooseRefFrame ~= 1)
        frameShown=chooseRefFrame;
        imagesc(img(:,:,chooseRefFrame)); axis('image')
        chooseRefFrame = input('Enter "1" if this frame is okay to draw contours on else enter another frame to try ');
end;
%s = fix(clock);   % date;
disp('Ref. frame chosen being written to log file');
% -- To log the changes made -- %
    f = fopen('log.txt','a+');
    fprintf(f,'Date: ');
    fprintf(f,' %d,\n',fix(clock));
    fprintf(f,'study: ');
    fprintf(f,'%d,\n',studyNum);
    fprintf(f,'slice: ');
    fprintf(f,'%d,\n',sliceNum);
    fprintf(f,'frameUsedforAIFin3.105 = frameShown =');
    fprintf(f,'%d,\n\n',frameShown);
set(gcf,'Renderer','zbuffer');   % to stop the blinking. Trying 3/22/0

disp('Use  mouse to choose Blood contour. Backspace or delete to redo last polygon vertex point')
[bw3, x3, y3] = mpi_roipoly;   %  polygon vertices
polyCoords = [x3 y3];
tmpname = strcat(outpath,'blood_polyCoords.study',int2str(studyNum),'.slice',int2str(sliceNum));
save(tmpname, 'polyCoords','-ascii');
hold
line(x3,y3,'Color', [0.3 0.8 0.5],'Linewidth',1); 
