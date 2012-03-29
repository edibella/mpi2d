function mpi_showRegions(img,sliceNum, studyNum, outpath, numAzimuthalRegions, numRadialRegions, referenceFrame)

% read in data 
[nrows, ncols, nFrames]=size(img)

clf

tmpImg=img(:,:,referenceFrame);
%figure(3); clf;
%imagesc(tmpImg);   axis('equal')
figure(1)
imagesc(tmpImg);
axis('image')
colormap gray

%nrows=2*(2*nrows-1)-1;
%ncols=2*(2*ncols-1)-1; % if odd same? yes.
%newImg=interp2(tmpImg,2,'cubic');  so size is 2*(2*nrows-1)-1 x 
% same as xi,yi=meshgrid(1:.25:60,...

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
set(h1,'Linewidth',2.5)
%set(h1,'Color',[1 1 0]); set(h1,'Linestyle','x'); set (h1,'Linestyle','-'); 
set(h1,'Color',[1 1 0]);
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
set(h2,'Linewidth',2.5)
%set(h2,'Color',[1 1 0]); set(h2,'Linestyle','pentagram'); set (h2,'Linestyle','-'); 
%set(h2,'Color',[0 1 0]); set(h2,'Linestyle','pentagram'); set (h2,'Linestyle','-'); 
center2=[mean(x2poly) mean(y2poly)];

%poly_length = min(length(x1poly),length(x2poly))
%xpoly_midline=(x1poly(1:poly_length)+x2poly(1:poly_length))/2;
%ypoly_midline=(y1poly(1:poly_length)+y2poly(1:poly_length))/2;
%h3=line(xpoly_midline,ypoly_midline);
%set(h3,'Color',[1 1 0]); 
%disp('look!')
%pause


tmp = load(strcat(outpath,'Roi_start_angle.study',int2str(studyNum),'.slice', int2str(sliceNum)));
start_angle = tmp;

%-load blood Roi-
   tmp=load(strcat(outpath,'blood_polyCoords.study',int2str(studyNum),'.slice',int2str(sliceNum)));
x3 = tmp(:,1); y3 = tmp(:,2);
bw3 = mpi_roipoly(tmpImg,x3,y3);
line(x3, y3, 'Color', [0.5 0.3 0.9]);

try
   if numRadialRegions>1
      tmp=load(strcat(outpath,'midwall_polyCoords.study',int2str(studyNum),'.slice',int2str(sliceNum)));
   end
catch
    disp('cannot find midwall  need to make this more graceful')	
end
x4 = tmp(:,1); y4 = tmp(:,2);
bw4 = mpi_roipoly(tmpImg,x4,y4);
line(x4, y4, 'Color', [0.7 0.2 0.4]);
%line(x4, y4, 'Color', [1 1 0]);
h4=line(x4,y4);
set(h4,'Color',[0 1 1]); set(h4,'Linestyle','pentagram'); set (h4,'Linestyle','-'); 


center=(center1+center2)/2;
xCenter = center(1);
yCenter = center(2);

for i=1:length(I)
    dist(i)=sqrt((J(i)-center(1))*(J(i)-center(1))+(I(i)-center(2)*(I(i)-center(2))));
end
maxdist=max(dist);
j=1;
clear bld_xCoords bld_yCoords;
for i=1:length(I)
    % check distance from center - move in 40% to just get blood
    if( sqrt((J(i)-center(1))*(J(i)-center(1))+(I(i)-center(2))*(I(i)-center(2))) < 0.4*maxdist)   % was .5
        bld_xCoords(j)=I(i);
        bld_yCoords(j)=J(i);
%        if( mod(I(i)+J(i),2) )
%	   tmpImg(I(i),J(i))=0;    % do some kind of overlay here! - matlab supports!
%        end
        j=j+1;
    end
end
imagesc(tmpImg);
%axis('square')
axis('image')
%colormap gray
h1=line(x1poly,y1poly);
set(h1,'Color',[0 1 0]); 
set(h1,'Linewidth',2.5)
h2=line(x2poly,y2poly);
set(h2,'Color',[1 1 0]); 
set(h2,'Linewidth',2.5)

line(x3, y3, 'Color', [0.5 0.3 0.9]);
h4=line(x4,y4);
set(h4,'Color',[1 1 0]); 
set(h4,'Linewidth',2.5)

% create a mask to overlay where blood coming from - better a contour:
% hmm, should convert to polar for this!
 
clear bld_xCoords bld_yCoords;


myo = double(bw2)-double(bw1);   % check all are 0 or 1 ... (2nd mask encompasses first!)
[I, J]=find(myo);
[Y, X] = find(myo);
nX     = length(X);
delta_t=1.8;
times   = (0 : nFrames - 1) * delta_t;




dAng  = 360 / numAzimuthalRegions;
mask = zeros(nrows, ncols, numAzimuthalRegions);
Angles =  atan2(-(Y - yCenter), X - xCenter) * 180/pi + 180;
for i = 1 : numAzimuthalRegions
   clear RoiX RoiY;
   cn = 0;
   Ang0 = (i - 1) * dAng + start_angle
   Ang1 = i * dAng       + start_angle
   for j = 1 : nX
      if( ((Angles(j) > Ang0) & (Angles(j) <= Ang1)) | ((Angles(j) + 360 > Ang0) & (Angles(j) + 360 <= Ang1)))
         cn = cn + 1;
         RoiX(cn) = X(j);
         RoiY(cn) = Y(j);
% want pixel with closest angle to angle 1 to be epicardium??? no - already have the contours! Do this another way!
%         if (Angles(j)==Ang1)
%            radDistToEpi(i)=sqrt((X(j)-xCenter)*(X(j)-xCenter)+(Y(j)-yCenter)*(Y(j)-yCenter));
%         end
      end
   end

   for n = 1 : cn
      mask(RoiY(n), RoiX(n), i) = 1;
   end

%   plot(RoiX, RoiY, '.', 'Color', [0.5, 0.7, 0.9]);
   rr = nrows/2;
   xx = [xCenter,    -rr * cos(Ang1*pi/180) + xCenter] ;
   yy = [yCenter,     rr * sin(Ang1*pi/180) + yCenter] ;
%   xx=[xCenter, xCenter-(radDistToEpi(i)*cos(Ang1*pi/180+pi/numAzimuthalRegions)) ] ;
%   yy = [yCenter, yCenter+(radDistToEpi(i)*sin(Ang1*pi/180+pi/numAzimuthalRegions))];
%   line(xx, yy,'Color', 'm');
%	hh=line(xx,yy);
% find intersection with contour:

%that is, given angle, what point out of x2poly, y2poly is closest to that ray? 
%   tmpp=atan2(x2poly-xCenter,y2poly-yCenter)*180/pi+180;
   tmpp=-atan2(y2poly-yCenter,x2poly-xCenter)*180/pi+180;
%   [closestAngle, indexIntersect]=sort(abs(tmpp-Ang1-90));
   [closestAngle, indexIntersect]=sort(abs(tmpp-Ang1));
   radDistToEpi=sqrt((x2poly(indexIntersect(1))-xCenter)*(x2poly(indexIntersect(1))-xCenter)+(y2poly(indexIntersect(1))-yCenter)*(y2poly(indexIntersect(1))-yCenter));
% changed to match applyRegions 8/17/05
xx=[xCenter, xCenter-1.2*(radDistToEpi*cos(Ang1*pi/180)) ] ;
yy = [yCenter, yCenter+1.2*(radDistToEpi*sin(Ang1*pi/180))];
	hh=line(xx,yy);

%        find xx in myo
	set(hh,'Color',[0 0 1])
	set(hh,'Linewidth',3)
%        xx=[xCenter-(rr/1.5)*cos(Ang1*pi/180+pi/numAzimuthalRegions) ] ;
%        yy = [yCenter+(rr/1.3)*sin(Ang1*pi/180+pi/numAzimuthalRegions)];
xx=[xCenter-1.4*(radDistToEpi*cos(Ang1*pi/180+pi/numAzimuthalRegions)) ] ;
yy = [yCenter+1.4*(radDistToEpi*sin(Ang1*pi/180+pi/numAzimuthalRegions))];
        if i==numAzimuthalRegions
           labelString=sprintf('%d',1)
        else
           labelString=sprintf('%d',i+1)
        end
        ht=text(xx,yy,labelString);
        set(ht,'FontSize',18,'Color',[1 1 1])

   t = mask(:,:,i);
   N = sum(sum(t));
end

figure(2)
for nframe=1:nFrames            % to display overlay of contours in each frame
   disp('Frame number '),nframe
   tmpImg=img(:,:,nframe);
   imagesc(tmpImg);
   h1=line(x1poly,y1poly);
   set(h1,'Color',[0 1 0]); 
   h2=line(x2poly,y2poly);
   set(h2,'Color',[1 1 0]); 
   h4=line(x4,y4);
   set(h4,'Color',[0.7 0.2 0.4]); 
   axis('image'); colormap gray
   pause   % if want to look at contours each time
end

%figure(2); clf; imagesc(sectorImage);

figure(1)

return;
