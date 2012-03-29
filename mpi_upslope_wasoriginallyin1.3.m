%function upslopes=mpi_upslope(sliceNum,studyNum, outpath, filenamePart1, numSkip, newRate, curves)
function upslopes=mpi_upslope(seriesNumAIF, sliceNumAIF, scaleAIF, sliceNum,studyNum, outpath, filenamePart1, delta_t, numAzimuthalRegions, numRadialRegions, flagPixelwise, numSkipUpslope, timeStampFile, curvesFromClusters)

% need to normalize by pre-contrast level  3/1/02

if nargin==0,
    error('sliceNum argument needed for mpi_upslope.');
end
if ~exist('numSkipUpslope'),   % if manually skipping an intital bump desired
   numSkipUpslope=0;
end

%if ~exist('newRate')
%	curvefilename=strcat(outpath,filenamePart1,'.study',int2str(studyNum),'.slice',int2str(sliceNum),'.mat');
%else
%	curvefilename=strcat(outpath,filenamePart1,'.study',int2str(studyNum),'.slice',int2str(sliceNum),'.every',int2str(newRate),'.mat');
%       disp('loading .every file')
%end

if ~exist('curvesFromClusters')
   curvefilename=strcat(outpath,filenamePart1,'.study',int2str(studyNum),'.slice',int2str(sliceNum),'.mat');
   load(curvefilename);
   curves = eval(filenamePart1);   % may want to re-do this so not so tricky
else
   curves=curvesFromClusters;
end

if (seriesNumAIF ~=-1)
   curvefilename=strcat(outpath,filenamePart1,'.study',int2str(studyNum),'.slice',int2str(sliceNum),'.AIF_',int2str(seriesNumAIF),'_',int2str(sliceNumAIF),'_',num2str(scaleAIF),'.mat');
   load(curvefilename);
   curves = eval(filenamePart1); 
end

bldcurve=curves(1,:)';
%gd_curves =[ bldGds'; tissGds'];  % so in same format as curves

nRegs=size(curves,1)-1
tisscurves=curves(2:nRegs+1,:);
tisscurve=tisscurves';
nTimes=length(bldcurve);


%gd_curves =[ bldGds'; tissGds'];  % so in same format as curves

[nRegs, nTimes]=size(curves)
bldcurve=curves(1,:)';

nRegs=nRegs-1;
tisscurves=curves(2:nRegs+1,:);
tisscurve=tisscurves';



% adding 11/15/05
% resample to uniform time points:
flagTimeStamps=0;
 if exist('timeStampFile')
   flagTimeStamps=1;
   load(timeStampFile);   % assumes gives variable timeStamp
   timeStamps=timeStamp(:,1)-timeStamp(1,1);   % slice 1 assumedly
   % units are in seconds after subtract off first time
   origSampleTimes= timeStamps; %normalizes so uniform delta_t sampling will have spacing of one
   if( nTimes < length(timeStamps))  %EVRD xx need to check this, added if 1/14/11
       origSampleTimes=origSampleTimes(2:nTimes+1);
   end
   figure(1); clf; hold on
plot(delta_t*(0:nTimes-1),bldcurve,'-or')
plot(delta_t*(0:nTimes-1),tisscurve(:,1),'-o') % using curve 1 in case only one azimuthal region
plot(origSampleTimes,bldcurve,'-xm')
plot(origSampleTimes,tisscurve(:,1),'-xk') % using curve 1 in case only one azimuthal region
disp('is this next pause off?')
%pause
   bldcurve = interp1(origSampleTimes, bldcurve, delta_t*(0:nTimes-1),'splines',0);
   bldcurve=bldcurve';
   for ii=1:nRegs
      tisscurve(:,ii) = interp1(origSampleTimes,tisscurve(:,ii),delta_t*(0:nTimes-1),'splines',0);
   end
   plot(delta_t*(0:nTimes-1),bldcurve(:),'-dg')
    plot(delta_t*(0:nTimes-1),tisscurve(:,1),'-dc') % using curve 1 in case only one azimuthal region
    plot(delta_t*(0:nTimes-1),tisscurve(:,1),'y') % using curve 1 in case only one azimuthal region
    xlabel('seconds')
    figure(2); clf;
    %plot(origSampleTimes(2:nTimes) - origSampleTimes(1:nTimes-1))
    xlabel('frame number')
%pause
end

%showcurves(curves,curvefilename);

N=3;    % Schwitter2001 used 3 point linear fit in blood and 5 point in tissue
maxslope_bld=0; xvec=1:N; xvec=xvec';
for i=numSkipUpslope:nTimes-N
   bb = find(bldcurve(i:i+N-1,1)~=0);
   if (~isempty(bb))
      [tmpcoeffs,Serr] = lregress(bldcurve(i:i+N-1,1),xvec);
      slope=tmpcoeffs(2,1);
   else
      slope=0;
   end
   if slope > maxslope_bld
       maxslope_bld = slope;
       slope_location_bld=i;   % can also use Serr to check if decent fit
   end
end
figure(1); clf; hold on
   plot(1:nTimes,bldcurve(:,1),'r','Linewidth',2.5)
   tmpcoeffs = mpi_twoplot(bldcurve(slope_location_bld:slope_location_bld+N-1),xvec,slope_location_bld); % to see plots
   slope_location_bld


N=5;
xvec=1:N; xvec=xvec';
maxslope=zeros(1,nRegs);
slope_location=zeros(1,nRegs);
for ireg=1:nRegs
% linear fit of N points:   % see Al Saadi Circ 2000 and Scwhitter Circ 2001
   for i=numSkipUpslope:nTimes/2-N      % make sure nTimes even
%   for i=numSkipUpslope:nTimes-2-N      % make sure nTimes even
%     tmpcoeffs = twoplot(tisscurve(i:i+N-1,ireg),xvec); % to see plots
      bb = find(tisscurve(i:i+N-1,ireg)~=0);
      if (~isempty(bb))
          [tmpcoeffs,Serr] = lregress(tisscurve(i:i+N-1,ireg),xvec);
          slope=tmpcoeffs(2,1);
      else
          slope=0;
      end
      if slope > maxslope(ireg)
         maxslope(ireg) = slope;
         slope_location(ireg)=i;   % can also use Serr to check if decent fit
      end
   end

% this is to show upslope fit
figure(2);clf; hold on
   plot(1:nTimes,tisscurve(:,ireg),'g','Linewidth',2.5)
   offset=slope_location(ireg)
   maxslope(ireg)
   tmpcoeffs = mpi_twoplot(tisscurve(slope_location(ireg):slope_location(ireg)+N-1,ireg),xvec,offset); % to see plots
   disp('region is '),ireg
   set(gca,'FontSize',16); set(gcf,'Color',[1 1 1]);
   xlabel('Time frame'); ylabel('change in signal intensity')
   slope_location(ireg)


% first smooth curve, then find max in first half and call that the peak...
%   tmp=conv(gausswin(nTimes/10),tisscurve(:,ireg));
%   smoothtisscurve(:,ireg)=tmp(1:nTimes);

%   i=numSkipUpslope+3
%   tmp1=tisscurve(i-1,ireg)+tisscurve(i,ireg)+tisscurve(i+1,ireg);
%   tmp2=tisscurve(i-2,ireg)+tisscurve(i-1,ireg)+tisscurve(i,ireg);

%   while tisscurve(i,ireg)> tisscurve(i-1,ireg)
%   while smoothtisscurve(i,ireg) < smoothtisscurve(i-1,ireg)
%       slope_start=i
%       i=i+1;
%   end %while

%   ddt=smoothtisscurve(1+numSkipUpslope:nTimes,ireg)-smoothtisscurve(numSkipUpslope:nTimes-1,ireg)

%   [curvemax slope_top]=max(ddt)
%   slope_top=slope_top+numSkipUpslope;

% to find start
% average up to 10 values just before start... *(from older code...)

% then rise/run:
%   upslope(ireg)=(tisscurve(slope_top,ireg) - tisscurve(slope_start,ireg))/(slope_top-slope_start);

end   % of ireg loop

%gd_curves =[ bldGds'; tissGds'];  % so in same format as curves
upslopes=maxslope
slope_location

if exist('newRate')
	curvefilename=strcat(outpath,filenamePart1,'.study',int2str(studyNum),'.slice',int2str(sliceNum),'.every',int2str(newRate),'.upslopes');
elseif exist('curvesFromClusters') 
	curvefilename=strcat(outpath,filenamePart1,'.study',int2str(studyNum),'.slice',int2str(sliceNum),'.clusters','.upslopes');
elseif (flagPixelwise)
	curvefilename=strcat(outpath,filenamePart1,'.study',int2str(studyNum),'.slice',int2str(sliceNum),'.pixelwise.upslopes');
else
	curvefilename=strcat(outpath,filenamePart1,'.study',int2str(studyNum),'.slice',int2str(sliceNum),'.',int2str(numAzimuthalRegions),'.',int2str(numRadialRegions),'.upslopes');
end


% modiifying output name 9/29/05
clusterFlag=0; useIntegralLinearFit=0; fixedDelay=99; fixedVp=99; 
if exist('seriesNumAIF')
   outfilename=mpi_getFilename(sliceNum, studyNum, outpath, numAzimuthalRegions, numRadialRegions,  clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp, seriesNumAIF, sliceNumAIF, scaleAIF, flagTimeStamps)
else
   outfilename=mpi_getFilename(sliceNum, studyNum, outpath, numAzimuthalRegions, numRadialRegions,  clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp)
end
outfilename=strcat(outfilename,'.upslopes.txt');


fid=fopen(outfilename,'w');
maxslope=maxslope./maxslope_bld;
fprintf(fid,'Created %s ,  delta_t=%f\n',date, delta_t);

for ireg=1:nRegs
   fprintf(fid,'upslope %d   %d     %4.4f  \n',ireg,slope_location(ireg),maxslope(ireg));
end
fprintf(fid,'\nmean  and Std is %4.4f   %4.4f,  coeff var = %6.3f \n',mean(maxslope), std(maxslope), std(maxslope)/mean(maxslope));

fprintf(fid,'\nAbove already normalized by blood slope which is %f  at index %d\n',maxslope_bld, slope_location_bld);
fclose(fid);

%showcurves([bldcurve'; tisscurve'],'DeltaSI or [Gd]');

return ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
