function deltaSI_curves=mpi_si2gdconc(sliceNum,studyNum, outpath, numPreContrast_bld, numPreContrast_tiss, numSkip, scale_to_study1, ncols, nrows, flagUseModel, curves)


tmpImg=zeros(nrows,ncols);
scale_img=ones(nrows,ncols);
init_img=zeros(nrows,ncols);


if nargin==0,
    error('sliceNum argument needed for mpi_si2gdconc.');
end
if ~exist('numPreContrast_bld'),
   numPreContrast_bld=7;
end
if ~exist('numPreContrast_tiss'),
   numPreContrast_tiss=9;
end
if ~exist('numSkip'),
   numSkip=0;
end
if ~exist('scale_to_study1'),   % this not iin place YET!
   scale_to_study1=1;
end

if ~exist('curves'),
   curvefilename=strcat(outpath,'curves.study',int2str(studyNum),'.slice',int2str(sliceNum),'.mat');
   load(curvefilename);
end
   si_curves = curves;
   [nRegs, nTimes]=size(curves)
   bldcurve=si_curves(1,:)';
   nRegs=nRegs-1;
   tisscurves = si_curves(2:nRegs+1,:);
   tisscurve=tisscurves';
   init_bld=sum(bldcurve(numSkip+1:numPreContrast_bld+numSkip))/numPreContrast_bld;
   tiss_avgSpatial=0; counter=0
   for ireg=1:nRegs
      init_tiss(ireg)=sum(tisscurve(numSkip+1:numPreContrast_tiss+numSkip,ireg))/numPreContrast_tiss;
      if(init_tiss(ireg)>0) 
         tiss_avgSpatial=tiss_avgSpatial+init_tiss(ireg);
         counter=counter+1
      end
   end
%   tiss_avgSpatial=mean(init_tiss);
   tiss_avgSpatial=tiss_avgSpatial/counter;




[X,Y,xCenter, yCenter, start_angle,bw1, bw2,bw3]= mpi_loadContours(tmpImg,sliceNum,studyNum, outpath);


if flagUseModel==1
disp('Using fits from before and init_Img to create truth')
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

[I, J]=find(bw1);
   [Ybld, Xbld] = find(bw1);
   for i=1:length(Xbld)
      scale_img(Ybld(i),Xbld(i))=1;
      init_img(Ybld(i),Xbld(i))=0;
   end
nX     = length(X);
if nX==nRegs   
   disp('reading in scale and init images pixelwise')
   for i=1:length(X)
      init_tiss(i)=init_img(Y(i),X(i));
   end
else   % then not doing pixelwise  

 nSect = 8;
 dAng  = 360 / nSect;
 Angles =  atan2(-(Y - yCenter), X - xCenter) * 180/pi + 180;
 for i = 1 : nSect
   clear RoiX RoiY;
   cn = 0;
   Ang0 = (i - 1) * dAng + start_angle;
   Ang1 = i * dAng       + start_angle;
   for j = 1 : nX
      if( ((Angles(j) > Ang0) & (Angles(j) < Ang1)) | ((Angles(j) + 360 > Ang0)& (Angles(j) + 360 < Ang1)))
         cn = cn + 1;
         RoiX(cn) = X(j);
         RoiY(cn) = Y(j);
      end
   end

   for n = 1 : 1  %cn
      init_tiss(i) = init_img(RoiY(n), RoiX(n));
   end

 end
end
tiss_avgSpatial=mean(init_tiss);
end % if flagUseModel


% subtract and divide for coil sens. (?) assuming bld same??
deltaSI_bldcurve= (bldcurve - init_bld);
deltaSI_bldcurve=tiss_avgSpatial*deltaSI_bldcurve./(tiss_avgSpatial*ones(length(bldcurve),1)) ;
ii=find(deltaSI_bldcurve<0);
%deltaSI_bldcurve(ii)=0;


for ireg=1:nRegs   % first subtract off (doesn't matter each region different scales!) Then normalize for coils. Other way, need to normalize init_tiss... 
% hmm, or scale, then re-calc. init_tiss, then subtract off new init_tiss 
% doesn't matter...  but when read in fits and reverse, they are 0 initially, so need to add back init and then scale...
%   deltaSI_tisscurve(:,ireg)=tisscurve(:,ireg);
%   if flagUseModel~=1
      deltaSI_tisscurve(:,ireg)=tisscurve(:,ireg)-init_tiss(ireg);
      deltaSI_tisscurve(:,ireg)=tiss_avgSpatial*deltaSI_tisscurve(:,ireg)./(init_tiss(ireg)*ones(length(tisscurve(:,ireg)),1)) ;
%   end
end
ii=find(deltaSI_tisscurve<0);
%deltaSI_tisscurve(ii)=0;

showcurves([deltaSI_bldcurve'; deltaSI_tisscurve'],'deltaSI scaled for coil');
disp('pause to look at curves - why go away?')
pause
deltaSI_curves =[ deltaSI_bldcurve'; deltaSI_tisscurve'];  % so in same format as curves



if flagUseModel~=1
% next section just to calculate scale and init_images to write out
nX     = length(X);
% should check nX==nRegs  (will for pixelwise, not for large regions!)
% for regions, need code that creates them...
if nX==nRegs   
   disp('writing out scale and init images pixelwise')
   for i=1:length(X)
      scale_img(Y(i),X(i))=init_tiss(i)/tiss_avgSpatial;  % assume in same order
      init_img(Y(i),X(i))=init_tiss(i);
  end

else   % then not doing pixelwise  

% add this 4/23/04 to put scale in 8 regions
nSect = 8;
dAng  = 360 / nSect;
Angles =  atan2(-(Y - yCenter), X - xCenter) * 180/pi + 180;
for i = 1 : nSect
   clear RoiX RoiY;
   cn = 0;
   Ang0 = (i - 1) * dAng + start_angle;
   Ang1 = i * dAng       + start_angle;
   for j = 1 : nX
      if( ((Angles(j) > Ang0) & (Angles(j) < Ang1)) | ((Angles(j) + 360 > Ang0) & (Angles(j) + 360 < Ang1)))
         cn = cn + 1;
         RoiX(cn) = X(j);
         RoiY(cn) = Y(j);
      end
   end

   for n = 1 : cn
      scale_img(RoiY(n), RoiX(n)) = init_tiss(i)/tiss_avgSpatial;  % assume in same order
      init_img(RoiY(n), RoiX(n)) = init_tiss(i); 
   end
end
end   % of if pixelwise


% adding 5/12/04 just to do blasted 32 region case!! 
% need to make this a function!!!!  Used in mpi_re_* and mpi_remake*
% stick in 5/12/04 for model with 32 regions:
 if (nRegs==32)
disp('doing 32 model')
   nSect=16;
   dAng  = 360 / nSect;
   endoFlag=1;
   [Xendo,Yendo,xCenter, yCenter, start_angle,bw1,bw2, bw3]= mpi_loadContours(tmpImg,sliceNum, studyNum,outpath,endoFlag);
   endoFlag=2;
   [Xepi,Yepi,xCenter, yCenter, start_angle,bw1, bw2,bw3]= mpi_loadContours(tmpImg,sliceNum, studyNum,outpath, endoFlag);

   for i = 1 : nSect
   nX=length(Xendo);
   clear RoiX RoiY;
   cn = 0;
   Ang0 = (i - 1) * dAng + start_angle;
   Ang1 = i * dAng       + start_angle;
   Angles =  atan2(-(Yendo - yCenter), Xendo - xCenter) * 180/pi + 180;
   for j = 1 : nX
      if( ((Angles(j) > Ang0) & (Angles(j) < Ang1)) | ((Angles(j) + 360 > Ang0)&(Angles(j) + 360 < Ang1)))
         cn = cn + 1;
         endoRoiX(cn) = Xendo(j);
         endoRoiY(cn) = Yendo(j);
      end
   end
   for n = 1 : cn
      scale_img(endoRoiY(n), endoRoiX(n)) = init_tiss(i)/tiss_avgSpatial;  % assume in same order
      init_img(endoRoiY(n), endoRoiX(n)) = init_tiss(i); 
   end

   nX=length(Xepi);
   clear RoiX RoiY;
   cn = 0;
   Angles =  atan2(-(Yepi - yCenter), Xepi - xCenter) * 180/pi + 180;
   for j = 1 : nX
      if( ((Angles(j) > Ang0) & (Angles(j) < Ang1)) | ((Angles(j) + 360 > Ang0)&(Angles(j) + 360 < Ang1)))
         cn = cn + 1;
         epiRoiX(cn) = Xepi(j);
         epiRoiY(cn) = Yepi(j);
      end
   end
   for n = 1 : cn
      scale_img(epiRoiY(n), epiRoiX(n)) = init_tiss(i)/tiss_avgSpatial;  % assume in same order
      init_img(epiRoiY(n), epiRoiX(n)) = init_tiss(i); 
   end
  end
  end  % if 32




% do blood same way if pixelwise or not:
   [Ybld, Xbld] = find(bw1);
   for i=1:length(Xbld)
      scale_img(Ybld(i),Xbld(i))=1;
      init_img(Ybld(i),Xbld(i))=init_bld;
   end



outfilename=strcat(outpath,'scaleImg.study',int2str(studyNum),'.slice',int2str(sliceNum),'.prebld',int2str(numPreContrast_bld),'.pretiss',int2str(numPreContrast_tiss),'.float');
ff=fopen(outfilename,'w');
fwrite(ff,scale_img','float'); 
fclose(ff);
outfilename=strcat(outpath,'initImg.study',int2str(studyNum),'.slice',int2str(sliceNum),'.prebld',int2str(numPreContrast_bld),'.pretiss',int2str(numPreContrast_tiss),'.float');
ff=fopen(outfilename,'w');
fwrite(ff,init_img','float'); 
fclose(ff);

end  % if flagUseModel


return;
