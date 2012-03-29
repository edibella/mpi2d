function mpi_writeParamImages(img,sliceNum, studyNum, outpath,numAzimuthalRegions, numRadialRegions, referenceFrame, flagPixelwise, numPreContrast_bld, numPreContrast_tiss)


if nargin==0,
    error('arguments needed for mpi_writeOParamImgaes');
 end
if ~exist('numPreContrast_bld'),
   numPreContrast_bld=7;
end
if ~exist('numPreContrast_tiss'),
   numPreContrast_tiss=9;
end


% read in data
[nrows, ncols, nFrames]=size(img)
clf
figure(1)

tmpImg=img(:,:,referenceFrame);
imagesc(tmpImg);
hold on
axis('square')


[X,Y,xCenter, yCenter, start_angle,bw1, bw2, bw3]= mpi_loadContours(img,sliceNum, studyNum,outpath);


%---blood curve------------------- 
N = sum(sum(bw3));
for j = 1 : nFrames
   a = img(:, :, j) .* bw3;
   bldcurve(j) = sum(sum(a)) / N;
end
%------------------------------------------

filenamePart1='deltaSIcurves'
if ~exist('curvesFromClusters')
   curvefilename=strcat(outpath,filenamePart1,'.study',int2str(studyNum),'.slice',int2str(sliceNum),'.mat');
   load(curvefilename);
   curves = eval(filenamePart1);   % may want to re-do this so not so tricky
else
   curves=curvesFromClusters;
end

bldcurve=curves(1,:)';   % this is deltaSI bldcurve



%-------------------------------------------



% read 1st line
if (flagPixelwise==1)
   filename=strcat(outpath,'fitparams.study',int2str(studyNum),'.slice',int2str(sliceNum),'.pixelwise.txt');
else
   filename=strcat(outpath,'fitparams.study',int2str(studyNum),'.slice',int2str(sliceNum),'.',int2str(numAzimuthalRegions),'.',int2str(numRadialRegions),'.txt');
end

fid=fopen(filename ,'r' );
a =fscanf(fid,'%s',3);
a =fscanf(fid,'%s',1);
delta_t=str2num(a(9:14));

indexkwi=1; indexkwo=1; indexfv=1; indext0=1;
while ~isempty(a)
  a =fscanf(fid,'%s',1);
  switch(a)
     case {'flow'}
        a =fscanf(fid,'%s',1);
        a =fscanf(fid,'%s',1);   % argh, must be a cleaner way!
        aa=str2num(a);
        kwi(indexkwi)=aa; indexkwi=indexkwi+1;
     case {'ve'}
        a =fscanf(fid,'%s',1);
        a =fscanf(fid,'%s',1);
        aa=str2num(a);
        kwo(indexkwo)=aa;  indexkwo=indexkwo+1;
     case {'fv'}
        a =fscanf(fid,'%s',1);
        a =fscanf(fid,'%s',1);
        aa=str2num(a);
        est_fb(indexfv)=aa;  indexfv=indexfv+1;
     case {'t0'}
        a =fscanf(fid,'%s',1);
        a =fscanf(fid,'%s',1);
        aa=str2num(a);
        t_delay(indext0)=aa;  indext0=indext0+1;
     otherwise
  end

end

kwi=kwi/2;
Hct=0.45;
kwo=kwi./((1-Hct)*kwo);


nSect=length(kwi);
nX     = length(X);
imgParams=zeros(nrows,ncols);

if (flagPixelwise~=1)
   dAng  = 360 / nSect;
   Angles =  atan2(-(Y - yCenter), X - xCenter) * 180/pi + 180; 
end

for i = 1 : nSect
  [fittedcurve, tmpbldcurve]=mpi_fwdModel(bldcurve, delta_t, kwi(i), kwo(i), t_delay(i), est_fb(i));


%%%%%%%%%%%%%%%%%%
  if (flagPixelwise)
      img(Y(i), X(i), :) = fittedcurve;
      imgParams(Y(i), X(i)) = kwi(i);
  else
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
      img(RoiY(n), RoiX(n), :) = fittedcurve;
      imgParams(RoiY(n), RoiX(n)) = kwi(i);
   end
  end
end  % for i=1:nSect

 if (nSect==32)
disp('doing 32 model')
   imgParams=zeros(nrows,ncols);
   nSect=16; 
   dAng  = 360 / nSect;
   endoFlag=1;
   [Xendo,Yendo,xCenter, yCenter, start_angle,bw1,bw2,  bw3]= mpi_loadContours(img,sliceNum, studyNum,outpath,endoFlag);
   endoFlag=2;
   [Xepi,Yepi,xCenter, yCenter, start_angle,bw1,bw2,  bw3]= mpi_loadContours(img,sliceNum, studyNum,outpath, endoFlag);

   for i = 1 : nSect
      [fittedcurve, tmpbldcurve]=mpi_fwdModel(bldcurve, delta_t, kwi(i), kwo(i), t_delay(i), est_fb(i));
      [fittedcurveEpi, tmpbldcurve]=mpi_fwdModel(bldcurve, delta_t, kwi(i+nSect), kwo(i+nSect), t_delay(i+nSect), est_fb(i+nSect));
   nX=length(Xendo);
   clear RoiX RoiY;
   cn = 0;
   Ang0 = (i - 1) * dAng + start_angle
   Ang1 = i * dAng       + start_angle
   Angles =  atan2(-(Yendo - yCenter), Xendo - xCenter) * 180/pi + 180; 
   for j = 1 : nX
      if( ((Angles(j) > Ang0) & (Angles(j) < Ang1)) | ((Angles(j) + 360 > Ang0) & (Angles(j) + 360 < Ang1)))
         cn = cn + 1;
         endoRoiX(cn) = Xendo(j); 
         endoRoiY(cn) = Yendo(j);
      end
   end
   for n = 1 : cn
      img(endoRoiY(n), endoRoiX(n), :) = fittedcurve;
      imgParams(endoRoiY(n), endoRoiX(n)) = kwi(i);
   end
  

   nX=length(Xepi);
   clear RoiX RoiY;
   cn = 0;
   Angles =  atan2(-(Yepi - yCenter), Xepi - xCenter) * 180/pi + 180; 
   for j = 1 : nX
      if( ((Angles(j) > Ang0) & (Angles(j) < Ang1)) | ((Angles(j) + 360 > Ang0) & (Angles(j) + 360 < Ang1)))
         cn = cn + 1;
         epiRoiX(cn) = Xepi(j); 
         epiRoiY(cn) = Yepi(j);
      end
   end
   for n = 1 : cn
      img(epiRoiY(n), epiRoiX(n), :) = fittedcurveEpi;  %(nSect-i)*10;
      imgParams(epiRoiY(n), epiRoiX(n)) = kwi(i+nSect);
   end
  end  % if 32

end  % for i=1:nSect

[Ybld, Xbld] = find(bw1); 
for j=1:length(Ybld)
%   img(Ybld(j),Xbld(j),:) = bldcurve;
end


   % doing some testing 5/10/04 evrd
%%     curves=[bldcurve, allfittedcurves]';
%%     fitparams2 = fit(sliceNum,studyNum,outpath,'deltaSIcurves',0,curves); 
%% keyboard



imgPreScale=img;

figure(2); clf; imagesc(img(:,:,referenceFrame));

infilename=strcat(outpath,'scaleImg.study',int2str(studyNum),'.slice',int2str(sliceNum),'.prebld',int2str(numPreContrast_bld),'.pretiss',int2str(numPreContrast_tiss),'.float');
ff=fopen(infilename,'r');
scale_img=fread(ff,[ncols, nrows],'float');
scale_img=scale_img';
fclose(ff);
infilename=strcat(outpath,'initImg.study',int2str(studyNum),'.slice',int2str(sliceNum),'.prebld',int2str(numPreContrast_bld),'.pretiss',int2str(numPreContrast_tiss),'.float');
ff=fopen(infilename,'r');
init_img=fread(ff,[ncols, nrows],'float');
init_img=init_img';
fclose(ff);

for j=1:length(Ybld)    % add this to try to leave original LV blood 
   init_img(Ybld(j),Xbld(j)) = 0;
end

for ii=1:nFrames   % reversed 4/30/04. Go without scaling  5/10/04
   img(:,:,ii) = img(:,:,ii).*scale_img;
   img(:,:,ii) = init_img + img(:,:,ii);
end

%disp('now hit key to see after scale and add ')

figure(2); clf; imagesc(img(:,:,referenceFrame));

% crop like do for clusters:
% find max and min X and Y to set crop boundaries:
   maxX=max(X)+1; minX=min(X)-1;
   maxY=max(Y)+1; minY=min(Y)-1;
   imgParams=imgParams(minY:maxY,minX:maxX,:);
[nrows ncols]=size(imgParams)


if (flagPixelwise)
   outfilename=strcat(outpath,'img.study',int2str(studyNum),'.slice',int2str(sliceNum),'newOriginal.pixelwise.washins.',int2str(ncols),'x',int2str(nrows),'.flt');
   outfilename2=strcat(outpath,'img.study',int2str(studyNum),'.slice',int2str(sliceNum),'newOriginal.pixelwise.float');
else
   %outfilename=strcat(outpath,'img.study',int2str(studyNum),'.slice',int2str(sliceNum),'PreScale.newOriginal.float');
   outfilename=strcat(outpath,'img.study',int2str(studyNum),'.slice',int2str(sliceNum),'newOriginal.washins.float');
   outfilename2=strcat(outpath,'img.study',int2str(studyNum),'.slice',int2str(sliceNum),'newOriginal.float');

end

ff=fopen(outfilename,'w')
ff2=fopen(outfilename2,'w')
for i=1:nFrames
   fwrite(ff2,img(:,:,i)','float');
end
fwrite(ff,imgParams','float');
fclose(ff);
%fwrite(ff2,imgParams','float');  % tack on parametric image to end
fclose(ff2);

return;
