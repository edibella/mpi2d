function curves=mpi_applyRegions(img,sliceNum, studyNum, outpath, referenceFrame, numAzimuthalRegions, flagPixelwise,numRadialRegions, pixelSizeRatio)
% read in epi, endo, bld contours and get out regional curves

% read in data 
[nrows, ncols, nFrames]=size(img)

clf
figure(1)
tmpImg=img(:,:,referenceFrame); % not sure if this is needed 
imagesc(tmpImg);
hold on
axis('image')


[X,Y,xCenter, yCenter, start_angle,bw1, bw2, bw3]= mpi_loadContours(img,sliceNum, studyNum, outpath);

%---blood curve------------------- 
N = sum(sum(bw3));
for j = 1 : nFrames
   a = img(:, :, j) .* bw3;
   bldcurve(j) = sum(sum(a)) / N;
end
%------------------------------------------


nX     = length(X);
if (flagPixelwise)
   for j = 1 : nX
      tisscurves(j,:) = img(Y(j), X(j), :);
   end
else


%% new 5/11/04 not implemented
%tisscurves=mpi_regions(numAzimuthalRegions,endoFlag,xCenter,yCenter,X,Y)

dAng  = 360 / numAzimuthalRegions;
mask = zeros(nrows, ncols, numAzimuthalRegions);
Angles =  atan2(-(Y - yCenter), X - xCenter) * 180/pi + 180; 
% if worry about pixel sizes
%pixelSizeRatio=pixelsizeX/pixelsizeY; 
%pixelSizeRatio=1.8229/2.983
%pixelSizeRatio=1.0
Angles =  atan2(-(Y - yCenter), pixelSizeRatio*(X - xCenter)) * 180/pi + 180; 
%%% NOTE - IF GET ERROR WITH RoiX, IS FROM START_ANGLE!!!  NOT ROBUST TO PICKING IN ANY  QUADRANT!!!
for i = 1 : numAzimuthalRegions
   clear RoiX RoiY;
   cn = 0;
   Ang0 = (i - 1) * dAng + start_angle;
   Ang1 = i * dAng       + start_angle;
   for j = 1 : nX
%      if( ((Angles(j) > Ang0) & (Angles(j) < Ang1)) | ((Angles(j) + 360 > Ang0) & (Angles(j) + 360 < Ang1)))
      if( ((Angles(j) > Ang0) & (Angles(j) < Ang1)) | ((Angles(j) + 360 > Ang0) & (Angles(j) + 360 < Ang1)) | ((Angles(j)+720 > Ang0) & (Angles(j)+720 < Ang1)))
         cn = cn + 1
         RoiX(cn) = X(j); 
         RoiY(cn) = Y(j);
      end
   end
   
   for n = 1 : cn
      mask(RoiY(n), RoiX(n), i) = 1;
   end
   
   
%   for j = 1 : nFrames
%      t = img(RoiY, RoiX, j);
%      curve1(i,j) = mean(mean(t));
%   end
   
   figure(1);
   plot(RoiX, RoiY, '.', 'Color', [0.5, 0.7, 0.9]);
   rr = nrows/2;
   rrx = pixelSizeRatio*nrows/2;
   rry = nrows/2;
   xx = [xCenter,    -rrx * cos(Ang1*pi/180) + xCenter] ;
   yy = [yCenter,     rry * sin(Ang1*pi/180) + yCenter] ;
   line(xx, yy,'Color', 'm');

%        find xx in myo
        xx=[xCenter-(rr/1.5)*cos(Ang1*pi/180+pi/numAzimuthalRegions) ] ;
        yy = [yCenter+(rr/1.3)*sin(Ang1*pi/180+pi/numAzimuthalRegions)];
        if i==numAzimuthalRegions
           labelString=sprintf('%d',1)
        else
           labelString=sprintf('%d',i+1)
        end
        ht=text(xx,yy,labelString);
        set(ht,'FontSize',18,'Color',[1 1 1])

   
%   figure(3);
%   plot(times, curve1(i,:), '--', 'Color', s(i).col);   
   
   t = mask(:,:,i);
   N = sum(sum(t));
   for j = 1 : nFrames
      a = img(:, :, j) .* t;
      tisscurves(i,j) = sum(sum(a)) / N;
   end
end
end   % if


% stick in 5/12/04 for model with 32 regions:
% if (numAzimuthalRegions==32 & ~flagPixelwise)
 if (numRadialRegions==2 & ~flagPixelwise)
disp('doing endo/epi regions ')
tmpname=strcat(outpath,'midwall_polyCoords.study',int2str(studyNum),'.slice',int2str(sliceNum));
if(~exist(tmpname))
   disp('Creating midwall automatically since no file found')
   mpi_computeMidwall(img,sliceNum,studyNum,outpath,referenceFrame);   % first compute and write out midwall contour
end

maskEndo = zeros(nrows, ncols, numAzimuthalRegions);
maskEpi = zeros(nrows, ncols, numAzimuthalRegions);
   dAng  = 360 / numAzimuthalRegions;
   endoFlag=1;
   [Xendo,Yendo,xCenter, yCenter, start_angle,bw1, bw2, bw3 ]= mpi_loadContours(img,sliceNum, studyNum,outpath,endoFlag);
   endoFlag=2;
   [Xepi,Yepi,xCenter, yCenter, start_angle,bw1, bw2, bw3]= mpi_loadContours(img,sliceNum, studyNum,outpath, endoFlag);

for i = 1 : numAzimuthalRegions
   mask=mpi_regions(Xcoords, Ycoords, numAzimuthalRegions, start_angle, pixelSizeRatio, totalRegion, xCenter,yCenter)

   nX=length(Xendo);
   clear RoiX RoiY;
   cn = 0;
   Ang0 = (i - 1) * dAng + start_angle;
   Ang1 = i * dAng       + start_angle;
%   Angles =  atan2(-(Yendo - yCenter), Xendo - xCenter) * 180/pi + 180;
% modify 5/24/04 evrd
   Angles =  atan2(-(Yendo - yCenter), pixelSizeRatio*(Xendo - xCenter)) * 180/pi + 180; 
   for j = 1 : nX
%      if( ((Angles(j) > Ang0) & (Angles(j) < Ang1)) | ((Angles(j) + 360 > Ang0)& (Angles(j) + 360 < Ang1)))
% modify 5/24/05 evrd
      if( ((Angles(j) > Ang0) & (Angles(j) < Ang1)) | ((Angles(j) + 360 > Ang0) & (Angles(j) + 360 < Ang1)) | ((Angles(j)+720 > Ang0) & (Angles(j)+720 < Ang1)))
         cn = cn + 1;
         endoRoiX(cn) = Xendo(j);
         endoRoiY(cn) = Yendo(j);
      end
   end
   for n = 1 : cn
      maskEndo(endoRoiY(n), endoRoiX(n),i) = 1;
   end
   t = maskEndo(:,:,i);
   N = sum(sum(t));
   for j = 1 : nFrames
      a = img(:, :, j) .* t;
      tisscurves(i,j) = sum(sum(a)) / N;
   end


   nX=length(Xepi);
   clear RoiX RoiY;
   cn = 0;
   Angles =  atan2(-(Yepi - yCenter), Xepi - xCenter) * 180/pi + 180;
   Angles =  atan2(-(Yepi - yCenter), pixelSizeRatio*(Xepi - xCenter)) * 180/pi + 180; 
   for j = 1 : nX
%      if( ((Angles(j) > Ang0) & (Angles(j) < Ang1)) | ((Angles(j) + 360 > Ang0)& (Angles(j) + 360 < Ang1)))
      if( ((Angles(j) > Ang0) & (Angles(j) < Ang1)) | ((Angles(j) + 360 > Ang0) & (Angles(j) + 360 < Ang1)) | ((Angles(j)+720 > Ang0) & (Angles(j)+720 < Ang1)))
         cn = cn + 1;
         epiRoiX(cn) = Xepi(j);
         epiRoiY(cn) = Yepi(j);
      end
   end
   for n = 1 : cn
      maskEpi(epiRoiY(n),epiRoiX(n),i) = 1;
   end

%  keyboard

   t = maskEpi(:,:,i);
   N = sum(sum(t));
   for j = 1 : nFrames
      a = img(:, :, j) .* t;
      tisscurves(i+numAzimuthalRegions,j) = sum(sum(a)) / N;
   end
end
end  % if 32

%end  % for i=1:numAzimuthalRegions

curves=[bldcurve;tisscurves];
return;
