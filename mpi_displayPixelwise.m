% re_createRegions - hacked from drawEndoAndEpi.m - 6/01. 
% 
% Probably should integrate getting curves into drawEndo... and
% then this could be just to visualize regions over time

% add option to write out curve for each pixel in myocardium  9/03
function mpi_displayPixelwise(img,sliceNum, studyNum, outpath, curves,referenceFrame, numPreContrast_bld, numPreContrast_tiss, useIntegralLinearFit,fixedDelay, fixedVp, seriesNumAIF, sliceNumAIF, scaleAIF,modelType)

if ~exist('useIntegralLinearFit')
   useIntegralLinearFit=0;
end
if ~exist('fixedDelay')
    fixedDelay=0;
end
if ~exist('fixedVp')
    fixedVp=0;
end
    
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




figure(1),clf
tmpImg=img(:,:,referenceFrame);
imagesc(tmpImg);
hold on
axis('square')

[X,Y,xCenter, yCenter, start_angle,bw1, bw3]= mpi_loadContours(img,sliceNum, studyNum,outpath);
nX     = length(X);


flagPixelwise=1;
numAzimuthalRegions=8; numRadialRegions=1;  % doesn't matter for pixelwise but needs to be defined 
clusterFlag=0;


if exist('seriesNumAIF')
   filename=mpi_getFilename(sliceNum, studyNum, outpath, numAzimuthalRegions, numRadialRegions,  clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp, seriesNumAIF, sliceNumAIF, scaleAIF)
else
   filename=mpi_getFilename(sliceNum, studyNum, outpath, numAzimuthalRegions, numRadialRegions,  clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp)
end

filename=strcat(filename,'.',modelType,'.txt')



[flow ve est_fb t_delay fval delta_t]=readFits(filename);
%%% Or you can comment the line above and manually load a flows file processed from
%%% other analysis methods (ie. Model-Independent analysis) to view here too.

% %%% ie. stress flows
%  load /v/raid6/ed/npack/SignalEqnProcessing/PixelwiseFlows/FLOW_P110906_17_2_pixelwise.mat
%  flow=FLOW_P110906_17_2_pixelwise(3,:);

% %%% ie. rest flows
% load /v/raid6/ed/npack/SignalEqnProcessing/PixelwiseFlows/FLOW_P110906_13_2_pixelwise.mat
% flow=FLOW_P110906_13_2_pixelwise(3,:);





times   = (0 : nFrames - 1) * delta_t;
kwi=flow*0.5;
Hct=0.45;
kwo=(1-Hct)*kwi./ve;   % divide by zero...
imgParams=zeros(size(tmpImg,1),size(tmpImg,2),5);

tmp=0:nTimes-1;
ttimes = delta_t * tmp';
ttimes = tmp';

imp=zeros(nTimes,1);
imp(1)=1.0;

ve_Threshold=1.0; % (ie. 1.0 represents 100% dist volume)


tmp=0:nTimes-1;
ttimes = delta_t * tmp';
ttimes = tmp';

imp=zeros(nTimes,1);
imp(1)=1.0;


%imgParams=zeros(size(tmpImg));
%imgParams=0.005*tmpImg;   % good for P09 not PCA
%imgParams=0.000005*tmpImg;
%imgParams=0.00005*tmpImg;
%clf; hold on;

for j = 1 : nX
%   [fittedcurve, tmpbldcurve]=mpi_fwdModel(bldcurve, delta_t, kwi(j), kwo(j), t_delay(j), est_fb(j));
%    imgParams(Y(j), X(j),1) = kwi(j);
    imgParams(Y(j), X(j),1) = flow(j);
    if ve(j) < ve_Threshold
       imgParams(Y(j), X(j),2) = ve(j); 
    else
       imgParams(Y(j), X(j),2) = 0;
    end
    imgParams(Y(j), X(j),3) = est_fb(j);
    imgParams(Y(j), X(j),4) = -t_delay(j);
%    imgParams(Y(j), X(j),5) = fval(j);
%    img(Y(j), X(j),:) = fittedcurve(1:nFrames);
end

%[Ybld, Xbld] = find(bw1);   % need to pass this!  fixed 4/23/04
%for j=1:length(Ybld)
%   img(Ybld(j),Xbld(j),:) = tmpbldcurve;
%end


figure(2); clf; % set(gca,'FontSize',16); set(gcf,'Color',[1 1 1]);
imagesc(imgParams(:,:,1)); axis('image')
%imagesc(imgParams(22:130,12:130,1)./2); axis('image') % these bounds were specific to the Ve Paper figures
%imagesc(imgParams(:,:,1)./2,[0 2]); axis('image') %this was used a grant image and was saved as a 1200dpi .tif file using 'print -dtiff -r1200 P101706_K1img'
title('Flow parameters in ml/min/g','FontSize',16)
spect=load('/home/mirl/ed/spect.cmap');
m=max(max(spect));
spect=spect/m;
colormap(spect)
hh=colorbar;
set(hh,'FontSize',11); set(hh,'FontWeight','demi'); 
%colormap('gray')
figure(3); clf
imagesc(imgParams(:,:,2)); axis('image')
%imagesc(imgParams(22:130,12:130,2),[0 0.35]); axis('image') % these bounds were specific to the Ve Paper figures
title('Ve in 1/min','FontSize',16)
hh=colorbar;
set(hh,'FontSize',11); set(hh,'FontWeight','demi'); 
figure(4); clf
imagesc(imgParams(:,:,3)); axis('image')
title('Vp (unitless)','FontSize',16)
colorbar;
figure(5); clf
imagesc(imgParams(:,:,4)); axis('image')
title('Tdelay in time frames, negative','FontSize',16)
colorbar;
figure(6)

disp('adjust and scale imgParams')
disp('may also want to smooth parametric iamge')
%keyboard


%for j=1:length(Ybld)    % add this to try to leave original LV blood
%   init_img(Ybld(j),Xbld(j)) = 0;
%end
%for ii=1:nFrames
%   img(:,:,ii) = init_img + img(:,:,ii);
%   img(:,:,ii) = img(:,:,ii).*scale_img;
%end

outfilename=strcat(outpath,'pixelwiseFits.study',int2str(studyNum),'.slice',int2str(sliceNum),'.',int2str(ncols),'x',int2str(nrows),'x5.flt');
  ff=fopen(outfilename,'w');
%  for i=1:nFrames         % transpose for display purposes
%      fwrite(ff,img(:,:,i)','float');  % get back origina orientation
%  end
%  fwrite(ff,1e4*imgParams','float'); % last image is just parameters
  for ii=1:4
      tmpp=imgParams(:,:,ii);
      fwrite(ff,tmpp','float');  % get back origina orientation
  end
  fclose(ff);




%for nframe=1:nFrames            % to display overlay of contours in each frame
%   tmpImg=img(:,:,nframe);
%   imagesc(tmpImg);
%   h1=line(x1poly,y1poly);
%   set(h1,'Color',[0 1 0]); 
%   h2=line(x2poly,y2poly);
%   set(h2,'Color',[1 1 0]); 
%   axis('square')
%end

   
%[regionXcoords, regionYcoords] = mpi_RegionCoords(img,sliceNum, studyNum, outpath, referenceFrame);
%outfilename=strcat(outpath,'regionCoords.study',int2str(studyNum),'.slice',int2str(sliceNum),'.short');
%save outfilename regionXcoords regionYcoords



return;
