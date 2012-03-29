% From function mpi_showRegions(img,sliceNum, studyNum, outpath, referenceFrame)
%function mpi_calcRegionsScales(img,sliceNum, studyNum, outpath, referenceFrame)
function [regionXcoords, regionYcoords] = mpi_RegionCoords(img,sliceNum, studyNum, outpath, referenceFrame)

%This is  to be called by mpi_display*.m and return vectors of pixel coords of each region
% where xcoords is 2D vector, each col. is a region, each row has x-coord of pixels in that region.  Or should it be a struct

% read in data 
[nrows, ncols, nFrames]=size(img)

clf

tmpImg=img(:,:,referenceFrame);
figure(3); clf;
imagesc(tmpImg);   %axis('square')
figure(1)
imagesc(tmpImg);
%axis('square')

% add error checking since did have different filename!!!
try
	tmp=load(strcat(outpath,'endo_polyCoords.study',int2str(studyNum),'.slice',int2str(sliceNum)));
catch
	disp('Dropped rest or stress part of filename March 2002. Will now attempt to read old filename')	
	tmp=load(strcat(outpath,'endo_polyCoords.rest.study',int2str(studyNum),'.slice',int2str(sliceNum)));
end
x1poly=tmp(:,1); y1poly=tmp(:,2);
bw1=mpi_roipoly(tmpImg,x1poly,y1poly);
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
center2=[mean(x2poly) mean(y2poly)];

try
   tmp = load(strcat(outpath,'Roi_start_angle.study',int2str(studyNum),'.slice',int2str(sliceNum)));
catch
   disp('Dropped rest or stress part of filename March 2002. Will now attempt to read old filename')
   tmp = load(strcat(outpath,'Roi_start_angle.rest.study',int2str(studyNum),'.slice',int2str(sliceNum)));
end
start_angle = tmp;


%-load blood Roi-
   tmp=load(strcat(outpath,'blood_polyCoords.study',int2str(studyNum),'.slice',int2str(sliceNum)));
x3 = tmp(:,1); y3 = tmp(:,2);
bw3 = mpi_roipoly(tmpImg,x3,y3);
line(x3, y3, 'Color', [0.5 0.3 0.9]);



center=(center1+center2)/2;
xCenter = center(1);
yCenter = center(2);


myo = double(bw2)-double(bw1);   % check all are 0 or 1 ... (2nd mask encompasses first!)
[I, J]=find(myo);
[Y, X] = find(myo);
nX     = length(X);

nSect = 8;
dAng  = 360 / nSect;
mask = zeros(nrows, ncols, nSect);
Angles =  atan2(-(Y - yCenter), X - xCenter) * 180/pi + 180;
for i = 1 : nSect
   clear RoiX RoiY;
   cn = 0;
   Ang0 = (i - 1) * dAng + start_angle
   Ang1 = i * dAng       + start_angle
   for j = 1 : nX
      if( ((Angles(j) > Ang0) & (Angles(j) < Ang1)) | ((Angles(j) + 360 > Ang0) & (Angles(j) + 360 < Ang1)))
         cn = cn + 1;
         RoiX(cn) = X(j);
         RoiY(cn) = Y(j);
      end
   end

   for n = 1 : cn
      mask(RoiY(n), RoiX(n), i) = 1;
      regionXcoords(i, n)=RoiY(n); 
      regionYcoords(i, n)=RoiX(n); 
   end

   rr = nrows/2;
   xx = [xCenter,    -rr * cos(Ang1*pi/180) + xCenter] ;
   yy = [yCenter,     rr * sin(Ang1*pi/180) + yCenter] ;
%   line(xx, yy,'Color', 'm');
	hh=line(xx,yy);
% find intersection with contour:
%        find xx in myo
	set(hh,'Color',[0 0 1])
	set(hh,'Linewidth',3)
        xx=[xCenter-(rr/1.5)*cos(Ang1*pi/180+pi/nSect) ] ;
        yy = [yCenter+(rr/1.3)*sin(Ang1*pi/180+pi/nSect)];
        if i==nSect
           labelString=sprintf('%d',1)
        else
           labelString=sprintf('%d',i+1)
        end
        ht=text(xx,yy,labelString);
        set(ht,'FontSize',18,'Color',[1 1 1])

   t = mask(:,:,i);
   N = sum(sum(t));
end

return;
