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
   ind=find(a>0); 
   bldcurve(j) = median(a(ind));
end
%------------------------------------------


nX     = length(X);
if (flagPixelwise)
   for j = 1 : nX
      tisscurves(j,:) = img(Y(j), X(j), :);
   end
% add 9/6/05 to model all pixels
%   X=1:size(img,1); Y=1:size(img,2); nX=length(X)
%numSkip=2;
%numPreContrast_bld=6;
%numPreContrast_tiss=7;
%   for j = 1 : size(img,1)
%    for k = 1 : size(img,2)
%      tisscurves((j-1)*size(img,2)+k,:) = img(j,k, :);
%      init_img(j,k)=mean(img(j,k, numSkip+1:numPreContrast_tiss+numSkip));
%    end
%   end
%outfilename=strcat(outpath,'initImg.study',int2str(studyNum),'.slice',int2str(sliceNum),'.prebld',int2str(numPreContrast_bld),'.pretiss',int2str(numPreContrast_tiss),'.float');
%ff=fopen(outfilename,'w','ieee-be');
%fwrite(ff,init_img','float');
%fclose(ff);


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
         cn = cn + 1;
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
% try median filter 9/5/05
      ind=find(a>0); 
      tisscurves(i,j) = median(a(ind));
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
   [endoRoiX, endoRoiY]= mpi_returnRegionCoords(i, Xendo, Yendo, numAzimuthalRegions, start_angle, pixelSizeRatio, xCenter,yCenter)

   for n = 1 : length(endoRoiX)
      maskEndo(endoRoiY(n), endoRoiX(n),i) = 1;
   end
   t = maskEndo(:,:,i);
   N = sum(sum(t));
   for j = 1 : nFrames
      a = img(:, :, j) .* t;
      tisscurves(i,j) = sum(sum(a)) / N;
% try median filter 9/5/05
      ind=find(a>0); 
      tisscurves(i,j) = median(a(ind));
   end

   [epiRoiX, epiRoiY]= mpi_returnRegionCoords(i, Xepi, Yepi, numAzimuthalRegions, start_angle, pixelSizeRatio, xCenter,yCenter)
   for n = 1 : length(epiRoiX)
      maskEpi(epiRoiY(n),epiRoiX(n),i) = 1;
   end
   t = maskEpi(:,:,i);
   N = sum(sum(t));
   for j = 1 : nFrames
      a = img(:, :, j) .* t;
      tisscurves(i+numAzimuthalRegions,j) = sum(sum(a)) / N;
% try median filter 9/5/05
      ind=find(a>0); 
      tisscurves(i+numAzimuthalRegions,j) = median(a(ind));
   end
end
end  % if 32


curves=[bldcurve;tisscurves];



% write out file for each region:

nRegs=numAzimuthalRegions*numRadialRegions;
for ii=1:nRegs
   outfilename=strcat(outpath,'tisscurve.study',int2str(studyNum),'.slice',int2str(sliceNum),'.',int2str(numAzimuthalRegions),'.',int2str(numRadialRegions),'_',int2str(ii),'.txt')
  ff=fopen(outfilename,'wt');
  fprintf(ff,'%6.4f\n',tisscurves(ii,:)); 
  fclose(ff);
end

return;
