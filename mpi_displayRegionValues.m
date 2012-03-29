function mpi_displayRegionValues(img,sliceNum, studyNum, outpath, numAzimuthalRegions, numRadialRegions, referenceFrame, clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp, seriesNumAIF, sliceNumAIF, scaleAIF, modelType, flagTimeStamps)

%function mpi_displayPixelwise(img,sliceNum, studyNum, outpath, curves,referenceFrame, numPreContrast_bld, numPreContrast_tiss, useIntegralLinearFit)

if ~exist('modelType')
   modelType='full';
end
if ~exist('flagTimeStamps')
    flagTimeStamps=1;
end
if ~exist('useIntegralLinearFit')
   useIntegralLinearFit=0;
end


% read in data 
[nrows, ncols, nFrames]=size(img);
nTimes=nFrames;

 clf

tmpImg=img(:,:,referenceFrame);
figure(1); clf; colormap gray;
imagesc(tmpImg);
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

%flagTimeStamps=0;
%flagTimeStamps=1;
if exist('seriesNumAIF') 
   filename=mpi_getFilename(sliceNum, studyNum, outpath, numAzimuthalRegions, numRadialRegions,  clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp, seriesNumAIF, sliceNumAIF, scaleAIF, flagTimeStamps);
%   filename=mpi_getFilename(sliceNum, studyNum, outpath, numAzimuthalRegions, numRadialRegions,  clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp, seriesNumAIF, sliceNumAIF, scaleAIF, flagTimeStamps)
else
   filename=mpi_getFilename(sliceNum, studyNum, outpath, numAzimuthalRegions, numRadialRegions,  clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp);
end
%filename=strcat(filename,'.no_ts.txt')

Hct=0.45;
if strcmp(modelType,'fermi')
   filename=strcat(filename,'.fermi.txt')
   [flow ve est_fb t_delay fval delta_t]=readFitsFermi(filename);
   kwi=flow;
   kwo=ve;
elseif strcmp(modelType,'fermiFull')
   filename=strcat(filename,'.fermiFull.txt')
   [flow ve est_fb t_delay fval delta_t]=readFitsFermi(filename);
   kwi=flow;
   kwo=ve;
elseif strcmp(modelType,'fermi2')
   filename=strcat(filename,'.fermi2.txt')
   [flow ve est_fb t_delay fval delta_t]=readFitsFermi(filename)
   kwi=flow;
   kwo=ve;
elseif strcmp(modelType,'fermiFull')
   filename=strcat(filename,'.fermiFull.txt')
   [flow ve est_fb t_delay fval delta_t]=readFitsFermi(filename)
   kwi=flow;
   kwo=ve;
elseif strcmp(modelType,'fermiFullTrunc')
   filename=strcat(filename,'.fermiFullTrunc.txt')
   [flow ve est_fb t_delay fval delta_t]=readFitsFermi(filename)
   kwi=flow;
   kwo=ve;
elseif strcmp(modelType,'globaldelay')
   filename=strcat(filename,'.globaldelay.txt')
   [flow ve est_fb t_delay fval delta_t]=readFits(filename)
   kwi=flow*0.5;
   kwo=(1-Hct)*kwi./ve;
elseif strcmp(modelType,'noBlood')
   filename=strcat(filename,'.noBlood.txt')
   [flow ve est_fb t_delay fval delta_t spillover]=readFits(filename)
   kwi=flow*0.5;
   kwo=(1-Hct)*kwi./ve;
elseif strcmp(modelType,'noBlood1')
   filename=strcat(filename,'.noBlood1.txt')
   [flow ve est_fb t_delay fval delta_t spillover]=readFits(filename)                                                                                                                        % NP changed to display 32 region Ve map for paper: ve=ve-0.05;
   kwi=flow*0.5;
   kwo=(1-Hct)*kwi./ve;
elseif strcmp(modelType,'blind')
   filename=strcat(filename,'.blind.txt')
   [flow ve est_fb t_delay fval delta_t spillover]=readFits(filename)
   kwi=flow*0.5;
   kwo=(1-Hct)*kwi./ve;
elseif strcmp(modelType,'blindKBI')
   filename=strcat(filename,'.blindKBI.txt')
   [flow ve est_fb t_delay fval delta_t spillover]=readFits(filename)
   kwi=flow*0.5;
   kwo=(1-Hct)*kwi./ve;
elseif strcmp(modelType,'full')
   filename=strcat(filename,'.full.txt')
   %[flow ve est_fb t_delay fval delta_t spillover]=readFits(filename)
   % kwi=flow*0.5;  % if wrote out not Ktrans but Ktrans/0.5
%    kwi=flow;    %EVRD 3/09
%    flow=kwi;   %EVRD 3/09
   
   [Ktrans ve est_fb t_delay fval delta_t spillover]=readFits(filename)
   kwi=Ktrans;Ktrans=kwi;
   
   kwo=(1-Hct)*kwi./ve;
elseif strcmp(modelType,'fullTrunc')
   filename=strcat(filename,'.fullTrunc.txt')
   [flow ve est_fb t_delay fval delta_t spillover]=readFits(filename)
   kwi=flow*0.5;
   kwo=(1-Hct)*kwi./ve;
else
   [flow ve est_fb t_delay fval delta_t spillover]=readFits(filename);
   kwi=flow*0.5;
   kwo=(1-Hct)*kwi./ve;
end
%[flow ve est_fb t_delay fval delta_t]=readFits(filename);


% if useIntegralLinearFit~=0
%    if exist('seriesNumAIF')
%       filename=strcat(outpath,'flowvalues.study',int2str(studyNum),'.slice',int2str(sliceNum),'.',int2str(numAzimuthalRegions), '.',int2str(numRadialRegions),'.linearFits','.AIF_',int2str(seriesNumAIF),'_',int2str(sliceNumAIF),'_', num2str(scaleAIF),'.txt')
%    else
%       filename=strcat(outpath,'fitparams.study',int2str(studyNum),'.slice',int2str(sliceNum),'.',int2str(numAzimuthalRegions), '.',int2str(numRadialRegions),'.linearFits.txt')
%    end
% 
% else
% %   filename=strcat(outpath,'fitparams.study',int2str(studyNum),'.slice',int2str(sliceNum),'.pixelwise','.txt')
%    if exist('seriesNumAIF')
%       filename=strcat(outpath,'flowvalues.study',int2str(studyNum),'.slice',int2str(sliceNum),'.',int2str(numAzimuthalRegions), '.',int2str(numRadialRegions),'.AIF_',int2str(seriesNumAIF),'_',int2str(sliceNumAIF),'_', num2str(scaleAIF),'.txt')
%    else
%       filename=strcat(outpath,'fitparams.study',int2str(studyNum),'.slice',int2str(sliceNum),'.',int2str(numAzimuthalRegions), '.',int2str(numRadialRegions),'.txt')
%    end
% 
%    if exist('seriesNumAIF')
% if fixedDelay==99 & fixedVp==99
%    filename=strcat(outpath,'flowvalues.study',int2str(studyNum),'.slice',int2str(sliceNum),'.',int2str(numAzimuthalRegions), '.',int2str(numRadialRegions),'.AIF_',int2str(seriesNumAIF),'_',int2str(sliceNumAIF),'_',num2str(scaleAIF),'.txt');
% elseif fixedDelay==99 & fixedVp~=99
%    filename=strcat(outpath,'flowvalues.study',int2str(studyNum),'.slice',int2str(sliceNum),'.',int2str(numAzimuthalRegions), '.',int2str(numRadialRegions),'.AIF_',int2str(seriesNumAIF),'_',int2str(sliceNumAIF),'_',num2str(scaleAIF),'_fixedVp',num2str(fixedVp),'.txt');
% elseif fixedDelay~=99 & fixedVp==99
%    filename=strcat(outpath,'flowvalues.study',int2str(studyNum),'.slice',int2str(sliceNum),'.',int2str(numAzimuthalRegions), '.',int2str(numRadialRegions),'.AIF_',int2str(seriesNumAIF),'_',int2str(sliceNumAIF),'_',num2str(scaleAIF),'_fixedDelay',num2str(fixedDelay),'.txt');
% elseif fixedDelay~=99 & fixedVp~=99
%    filename=strcat(outpath,'flowvalues.study',int2str(studyNum),'.slice',int2str(sliceNum),'.',int2str(numAzimuthalRegions), '.',int2str(numRadialRegions),'.AIF_',int2str(seriesNumAIF),'_',int2str(sliceNumAIF),'_',num2str(scaleAIF),'_fixedDelay',num2str(fixedDelay),'_fixedVp',num2str(fixedVp),'.txt');
% else
%   disp('error in fit_useAlt writing  out parameters')
% end
% 
%    else
% 
% if fixedDelay==99 & fixedVp==99
%    filename=strcat(outpath,'fitparams.study',int2str(studyNum),'.slice',int2str(sliceNum),'.',int2str(numAzimuthalRegions), '.',int2str(numRadialRegions),'.txt');
% elseif fixedDelay==99 & fixedVp~=99
%    filename=strcat(outpath,'flowvalues.study',int2str(studyNum),'.slice',int2str(sliceNum),'.',int2str(numAzimuthalRegions), '.',int2str(numRadialRegions),'_fixedVp',num2str(fixedVp),'.txt');
% elseif fixedDelay~=99 & fixedVp==99
%    filename=strcat(outpath,'flowvalues.study',int2str(studyNum),'.slice',int2str(sliceNum),'.',int2str(numAzimuthalRegions), '.',int2str(numRadialRegions),'_fixedDelay',num2str(fixedDelay),'.txt');
% elseif fixedDelay~=99 & fixedVp~=99
%    filename=strcat(outpath,'flowvalues.study',int2str(studyNum),'.slice',int2str(sliceNum),'.',int2str(numAzimuthalRegions), '.',int2str(numRadialRegions),'_fixedDelay',num2str(fixedDelay),'_fixedVp',num2str(fixedVp),'.txt');
% else
%   disp('error in fit writing  out parameters')
% end
% 
% end
% 
% 
% 
% end




% %%% Or you can  manually load a flows file processed from
% %%% other analysis methods (ie. Model-Independent analysis) to view here too.
% % %%% ie. stress flows (this example below load pixelwise data and may not work here)
% load /v/raid6/ed/npack/SignalEqnProcessing/TempFlowP110906_13_2.mat
% flow=TempFlowP110906_13_2(3,:);

% % % times   = (0 : nFrames - 1) * delta_t;
% % % kwi=flow*0.5;
% % % kwo=(1-Hct)*kwi./ve;   % divide by zero...
imgParams=zeros(size(tmpImg,1),size(tmpImg,2),5);
% % % 
% % % tmp=0:nTimes-1;
% % % ttimes = delta_t * tmp';
% % % ttimes = tmp';
% % % 
% % % imp=zeros(nTimes,1);
% % % imp(1)=1.0;


ve_Threshold=1;


dAng  = 360 / numAzimuthalRegions;
mask = zeros(nrows, ncols, numAzimuthalRegions);
Angles =  atan2(-(Y - yCenter), X - xCenter) * 180/pi + 180;
for i = 1 : numAzimuthalRegions
   clear RoiX RoiY;
   cn = 0;
   Ang0 = (i - 1) * dAng + start_angle;
   Ang1 = i * dAng       + start_angle;
   for j = 1 : nX
%      if( ((Angles(j) > Ang0) & (Angles(j) <= Ang1)) | ((Angles(j) + 360 > Ang0) & (Angles(j) + 360 <= Ang1)))
     if( ((Angles(j) > Ang0) & (Angles(j) < Ang1)) | ((Angles(j) + 360 > Ang0)& (Angles(j) + 360 < Ang1)) | ((Angles(j)+720 > Ang0) & (Angles(j)+720 < Ang1)))
         cn = cn + 1;
         RoiX(cn) = X(j);
         RoiY(cn) = Y(j);
      end
   end

   for n = 1 : cn
      %mask(RoiY(n), RoiX(n), i) = flow(i);
      %imgParams(RoiY(n), RoiX(n), 1) = flow(i);
      mask(RoiY(n), RoiX(n), i) = Ktrans(i);
      imgParams(RoiY(n), RoiX(n), 1) = Ktrans(i);
      if ve(i) < ve_Threshold
        imgParams(RoiY(n), RoiX(n), 2) = ve(i);
      else
        imgParams(RoiY(n), RoiX(n), 2) = 0;
      end
      imgParams(RoiY(n), RoiX(n), 3) = est_fb(i);
      imgParams(RoiY(n), RoiX(n), 4) = -t_delay(i);
%    imgParams(Y(j), X(j),5) = fval(j);

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
           labelString=sprintf('%d',1);
        else
           labelString=sprintf('%d',i+1);
        end
        ht=text(xx,yy,labelString);
        set(ht,'FontSize',18,'Color',[1 1 1]);

   t = mask(:,:,i);
   N = sum(sum(t));
end

%imgtmp=mask(:,:,1);
%imgtmp=imgtmp*0;
%for ii=1:numAzimuthalRegions
%  imgtmp=imgtmp+mask(:,:,ii);
%  imagesc(imgtmp);
%end



% if (numRadialRegions==2 & ~flagPixelwise)
if (numRadialRegions==2 )
   imgParams=zeros(size(tmpImg,1),size(tmpImg,2),5);
   pixelSizeRatio=1;
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
      [endoRoiX, endoRoiY]= mpi_returnRegionCoords(i, Xendo, Yendo, numAzimuthalRegions, start_angle, pixelSizeRatio, xCenter,yCenter);

      for n = 1 : length(endoRoiX)
         imgParams(endoRoiY(n), endoRoiX(n), 1) = flow(i);
         if ve(i) < ve_Threshold
           imgParams(endoRoiY(n), endoRoiX(n), 2) = ve(i);
         else
           imgParams(endoRoiY(n), endoRoiX(n), 2) = 0;
         end
         imgParams(endoRoiY(n), endoRoiX(n), 3) = est_fb(i);
         imgParams(endoRoiY(n), endoRoiX(n), 4) = -t_delay(i);
      end

      [epiRoiX, epiRoiY]= mpi_returnRegionCoords(i, Xepi, Yepi, numAzimuthalRegions, start_angle, pixelSizeRatio, xCenter,yCenter);
      for n = 1 : length(epiRoiX)
         imgParams(epiRoiY(n), epiRoiX(n), 1) = flow(i+numAzimuthalRegions);
         imgParams(epiRoiY(n), epiRoiX(n), 2) = ve(i+numAzimuthalRegions);
         imgParams(epiRoiY(n), epiRoiX(n), 3) = est_fb(i+numAzimuthalRegions);
         imgParams(epiRoiY(n), epiRoiX(n), 4) = -t_delay(i+numAzimuthalRegions);
      end
   end

end  % if 32


% 
mpi2d stg=6.1
%EVRD, hack this until fix below, not in right order...

% 
% figure(2); clf; % set(gca,'FontSize',16); set(gcf,'Color',[1 1 1]);
% imagesc(imgParams(:,:,1)); axis('image')
% %imagesc(imgParams(22:130,12:130,1)); axis('image') % these bounds were specific to the Ve Paper figures
% title('K^{trans}','FontSize',16)
% spect=load('spect.cmap');
% m=max(max(spect));
% spect=spect/m;
% colormap(spect)
% hh=colorbar;
% set(hh,'FontSize',11); set(hh,'FontWeight','demi');

% %colormap('gray')
% % figure(3); clf
% % colormap(spect)
% % imagesc(imgParams(:,:,2)); axis('image')
% % %imagesc(imgParams(22:130,12:130,2),[0 0.35]); axis('image') % these bounds were specific to the Ve Paper figures
% % title('Ve in 1/min','FontSize',16)
% % hh=colorbar;
% % set(hh,'FontSize',11); set(hh,'FontWeight','demi');
% % figure(4); clf
% % colormap(spect)
% % imagesc(imgParams(:,:,3)); axis('image')
% % title('Vp (unitless)','FontSize',16)
% % colorbar;
% % figure(5); clf
% % colormap(spect)
% % imagesc(imgParams(:,:,4)); axis('image')
% % title('Tdelay in time frames, negative','FontSize',16)
% % colorbar;



% figure(2);  % to bring to front
% figure(6)



return;
