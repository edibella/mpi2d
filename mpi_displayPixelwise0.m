% re_createRegions - hacked from drawEndoAndEpi.m - 6/01. 
% 
% Probably should integrate getting curves into drawEndo... and
% then this could be just to visualize regions over time

% add option to write out curve for each pixel in myocardium  9/03
function mpi_displayPixelwise(img,sliceNum, studyNum, outpath, curves,referenceFrame, numPreContrast_bld, numPreContrast_tiss)

% read in data 
[nrows, ncols, nFrames]=size(img)
nTimes=nFrames;
bldcurve=curves(1,:)';


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



clf
figure(1)
tmpImg=img(:,:,referenceFrame);
imagesc(tmpImg);
hold on
axis('square')

[X,Y,xCenter, yCenter, start_angle,bw1, bw3]= mpi_loadContours(img,sliceNum, studyNum,outpath);
nX     = length(X);


% read 1st line
filename=strcat(outpath,'fitparams.study',int2str(studyNum),'.slice',int2str(sliceNum),'.pixelwise','.txt')
%filename=strcat('fitparams.study',int2str(studyNum),'.slice',int2str(sliceNum),'.clusters','.txt');
fid=fopen(filename ,'r' );
a =fscanf(fid,'%s',3);
a =fscanf(fid,'%s',1);
delta_t=str2num(a(9:14));

indexkwi=1; indexkwo=1; indexfv=1; indext0=1; indexfval=1;
while ~isempty(a)
  a =fscanf(fid,'%s',1);
  switch(a)
     case {'kwi'}
        a =fscanf(fid,'%s',1);
        aa=str2num(a);
        kwi(indexkwi)=aa; indexkwi=indexkwi+1;
     case {'kwo'}
        a =fscanf(fid,'%s',1);
        aa=str2num(a);
        kwo(indexkwo)=aa;  indexkwo=indexkwo+1;
     case {'fv'}
        a =fscanf(fid,'%s',1);
        aa=str2num(a);
        est_fb(indexfv)=aa;  indexfv=indexfv+1;
     case {'t0'}
        a =fscanf(fid,'%s',1);
        aa=str2num(a);
        t_delay(indext0)=aa;  indext0=indext0+1;
     case {'fval'}
        a =fscanf(fid,'%s',1);
        aa=str2num(a);
        fval(indexfval)=aa;  indexfval=indexfval+1;
     otherwise
        disp('Ending loop'),a
  end

end

kwi;
kwo;
est_fb;
t_delay;
delta_t;
tmp=0:nTimes-1;
ttimes = delta_t * tmp';
ttimes = tmp';

imp=zeros(nTimes,1);
imp(1)=1.0;


imgParams=zeros(size(tmpImg));
%imgParams=0.005*tmpImg;   % good for P09 not PCA
%imgParams=0.000005*tmpImg;
imgParams=0.00005*tmpImg;
%clf; hold on;
for j = 1 : nX
   [fittedcurve, tmpbldcurve]=mpi_fwdModel(bldcurve, delta_t, kwi(j), kwo(j), t_delay(j), est_fb(j));
%    tmpbldcurve = interp1(1:nTimes,bldcurve,1-t_delay(j):nTimes-t_delay(j),'splines',0);
%    fittedcurve = conv(tmpbldcurve,((1-est_fb(j))*delta_t/60*kwi(j)*exp(-kwo(j)*delta_t/60*(ttimes-t_delay(j)*ones(nTimes,1))) + est_fb(j)*imp) );
%    fittedcurve = conv(tmpbldcurve, delta_t/60*((1-est_fb(j))*kwi(j)*exp(-kwo(j)*delta_t/60*(ttimes)) + est_fb(j)*imp) );
%kwi(j), kwo(j), est_fb(j), , t_delay(j)
%plot(fittedcurve(1:nFrames),'Linewidth',2)
%plot(squeeze(img(Y(j), X(j),:)),':r','Linewidth',2);
    imgParams(Y(j), X(j)) = kwi(j);
    img(Y(j), X(j),:) = fittedcurve(1:nFrames);
%pause
end

[Ybld, Xbld] = find(bw1);   % need to pass this!  fixed 4/23/04
for j=1:length(Ybld)
   img(Ybld(j),Xbld(j),:) = tmpbldcurve;
end


figure(2)
clf
imagesc(imgParams)
   h1=line(x1poly,y1poly);
   set(h1,'Color',[0 1 0]); 
   h2=line(x2poly,y2poly);
   set(h2,'Color',[1 1 0]); 
axis('square')
disp('adjust and scale imgParams')
disp('may also want to smooth parametric iamge')
%keyboard
figure(1)

for j=1:length(Ybld)    % add this to try to leave original LV blood
   init_img(Ybld(j),Xbld(j)) = 0;
end
for ii=1:nFrames
   img(:,:,ii) = init_img + img(:,:,ii);
   img(:,:,ii) = img(:,:,ii).*scale_img;
end

outfilename=strcat(outpath,'pixelwiseFits.study',int2str(studyNum),'.slice',int2str(sliceNum),'.float');
  ff=fopen(outfilename,'w');
  for i=1:nFrames         % transpose for display purposes
      fwrite(ff,img(:,:,i)','float');  % get back origina orientation
  end
  fwrite(ff,1e4*imgParams','float'); % last image is just parameters
  fclose(ff);




for nframe=1:nFrames            % to display overlay of contours in each frame
%   disp('Frame number '),nframe
   tmpImg=img(:,:,nframe);
   imagesc(tmpImg);
   h1=line(x1poly,y1poly);
   set(h1,'Color',[0 1 0]); 
   h2=line(x2poly,y2poly);
   set(h2,'Color',[1 1 0]); 
   axis('square')
%   pause   % if want to look at contours each time
end
%clf; hold on
%plot(1:19,tisscurves,'Linewidth',2)  % this chooses different colors
%ylabel('SI')
%xlabel('Frame number')

%figure(2); clf; imagesc(sectorImage);
   
%curves=[bldcurve;tisscurves];


[regionXcoords, regionYcoords] = mpi_RegionCoords(img,sliceNum, studyNum, outpath, referenceFrame);
outfilename=strcat(outpath,'regionCoords.study',int2str(studyNum),'.slice',int2str(sliceNum),'.short');
save outfilename regionXcoords regionYcoords



return;
