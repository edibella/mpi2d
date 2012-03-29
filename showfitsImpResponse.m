function showfitsImpResponse(curves, studyNum, sliceNum, outpath, numAzimuthalRegions, numRadialRegions, clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp, seriesNumAIF, sliceNumAIF, scaleAIF, modelType, lastFrame, timeStampFile, timeStampFileAIF)

if ~exist('modelType')
    modelType='xxx'
 %   modelType='stdModel'
end
if exist('timeStampFileAIF')  % I think assumes this always exists? haven't thought of will handle a call to this function without it
   flagTimeStamps=1;
   if strcmp(timeStampFileAIF,'') & exist('timeStampFile')
      timeStampFileAIF=timeStampFile;
   end
end

if ~exist('clusterFlag')
   clusterFlag=0;
end
if ~exist('flagPixelwise')
   flagPixelwise=0;
end
if ~exist('useIntegralLinearFit')
   useIntegralLinearFit=0;
end

[nRegs, nTimes]=size(curves)
bldcurve=curves(1,:)';

nRegs=nRegs-1;
tisscurves=curves(2:nRegs+1,:);
tisscurve=tisscurves';

% from older, should do try, catch
%filename=strcat(outpath,'fitparams.study',int2str(studyNum),'.slice',int2str(sliceNum),'.',int2str(numAzimuthalRegions),'.',int2str(numRadialRegions),'.txt');


%flagTimeSamples=0;
%flagTimeSamples=1;
if exist('seriesNumAIF')
   filename=mpi_getFilename(sliceNum, studyNum, outpath, numAzimuthalRegions, numRadialRegions, clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp, seriesNumAIF, sliceNumAIF, scaleAIF, flagTimeStamps)
   % 9/29/05 Also read in saturated AIF for modeling:
   filenamePart1='deltaSIcurves';
 %  curvefilename=strcat(outpath,filenamePart1,'.study',int2str(studyNum),'.slice',int2str(sliceNum),'.mat')
   curvefilename=strcat(outpath,filenamePart1,'.study',int2str(studyNum),'.slice',int2str(sliceNum),'.AIF_',int2str(studyNum),'_',int2str(sliceNum),'_',int2str(numRadialRegions),'.mat')
   load(curvefilename);
   curves = eval(filenamePart1);
   sat_bldcurve=curves(1,1:length(bldcurve))';
else
   filename=mpi_getFilename(sliceNum, studyNum, outpath, numAzimuthalRegions, numRadialRegions,  clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp);
end



Hct=0.45;
if strcmp(modelType(1:4),'ferm')   % so minimum length of modelType is 4 letters
   filename=strcat(filename,'.',modelType,'.txt')
   [flow ve est_fb t_delay T0 fval delta_t spillover]=readFitsFermi(filename)
   kwi=flow;
   kwo=ve;
%elseif strcmp(modelType,'fermiConstrDelay')
%   filename=strcat(filename,'.fermiConstrDelay.txt')
%   [flow ve est_fb t_delay fval delta_t]=readFitsFermi(filename)
%   kwi=flow;
%   kwo=ve;
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
   [flow ve est_fb t_delay fval delta_t spillover]=readFits(filename)
   kwi=flow*0.5;
   kwo=(1-Hct)*kwi./ve;
elseif strcmp(modelType,'full')
   filename=strcat(filename,'.full.txt')
   [flow ve est_fb t_delay fval delta_t spillover]=readFits(filename)
   kwi=flow*0.5;
   kwo=(1-Hct)*kwi./ve;
elseif strcmp(modelType,'fullTrunc')
   filename=strcat(filename,'.fullTrunc.txt')
   [flow ve est_fb t_delay fval delta_t spillover]=readFits(filename)
   kwi=flow*0.5;
   kwo=(1-Hct)*kwi./ve;
elseif strcmp(modelType,'stdModel')
%   filename=strcat(filename,'.txt')
   [flow ve est_fb t_delay fval delta_t spillover]=readFits(filename)
   kwi=flow*0.5;
   kwo=(1-Hct)*kwi./ve;
elseif strcmp(modelType,'fullModel2')
   filename=strcat(filename,'.fullModel2.txt')
   [flow ve est_fb t_delay fval delta_t spillover]=readFits(filename)
   kwi=flow*0.5;
   kwo=(1-Hct)*kwi./ve;
else
   [flow ve est_fb t_delay fval delta_t spillover]=readFits(filename)
   kwi=flow*0.5;
   kwo=(1-Hct)*kwi./ve;
end

delta_t

tmp=0:nTimes-1;
ttimes = delta_t * tmp';
ttimes = tmp';

  if strcmp(modelType,'fermiFullTrunc') || strcmp(modelType,'fullTrunc') || strcmp(modelType,'fermiTrunc')
      nTimes=lastFrame;
      tisscurve=tisscurve(1:lastFrame,:);
      bldcurve=bldcurve(1:lastFrame);
      sat_bldcurve=sat_bldcurve(1:lastFrame);
      tmp=0:nTimes-1;
      ttimes = delta_t * tmp';
      ttimes = tmp';
   end


% resample to uniform time points:
flagTimeStamps=0;
if exist('timeStampFile')
   flagTimeStamps=1;
   bldcurve=interpTimeCurve(timeStampFileAIF, delta_t, bldcurve);
   nTimes=length(bldcurve);

% now do different one for Tissues:
   origtisscurve=tisscurve;
   clear tisscurve
   for ii=1:nRegs
       tisscurve(:,ii)=interpTimeCurve(timeStampFile, delta_t, origtisscurve(:,ii));
   end
 %  pause  % to see last region of tisscurve
   sat_bldcurve=interpTimeCurve(timeStampFile, delta_t, sat_bldcurve);

% add sanity check to make sure AIF and tissue curves are same length, since they
% could come from different timeStamp files
   if (length(bldcurve) > size(tisscurve,1))
       bldcurve=bldcurve(1:size(tisscurve,1));  % tr
       disp('Bld and tiss curves lenghts do not match due to timestamps...\n')
       disp('Truncating bldcurve length to match tisscurve')
       nTimes=size(tisscurve,1);
   end
   if (length(bldcurve) < size(tisscurve,1))
       tisscurve=tisscurve(1:nTimes,:);
       disp('Bld and tiss curves lenghts do not match due to timestamps...\n')
       disp('Truncating tisscurve length to match bldcurve')
   end
end  


numSkip=1
numSkip=0
bldcurve=bldcurve(numSkip+1:nTimes);
nTimes=nTimes-numSkip;
sat_bldcurve=bldcurve;

   tmp=0:nTimes-1;
      ttimes = delta_t * tmp';
      ttimes = tmp';


jitter=0; newRate=1;

errsum = zeros(1,nRegs);
for i=1:nRegs
    nn(:,i)=tisscurve(numSkip+1:nTimes+numSkip,i);
end
clear tisscurve;
tisscurve=nn;
 
for i=1:nRegs
   disp('iReg is '),i
  

bldcurve=0*bldcurve;
bldcurve(10)=1;  % impulse response
sat_bldcurve=0*sat_bldcurve;
sat_bldcurve(10)=1;  % impulse response

   
%   tisscurve(:,i)=nn;
   if strcmp(modelType(1:4),'ferm')
%      [est_curves(:,i), tmpbldcurve]=mpi_fermi_fwdModel(bldcurve, sat_bldcurve, delta_t, kwi(i), kwo(i),t_delay(i), est_fb(i));
      [est_curves(:,i), tmpbldcurve]=mpi_fermi_fwdModel(bldcurve, sat_bldcurve, delta_t, kwi(i), kwo(i), t_delay(i), T0(i), est_fb(i), spillover(i));
   elseif strcmp(modelType,'globaldelay')
      [est_curves(:,i), tmpbldcurve]=mpi_fwdModel(bldcurve, delta_t, kwi(i), kwo(i), t_delay(i), est_fb(i));
   elseif strcmp(modelType,'noBlood')
      [est_curves(:,i), tmpbldcurve]=mpi_sat_fwdModel(bldcurve, sat_bldcurve, delta_t, kwi(i), kwo(i), t_delay(i), est_fb(i), spillover(i));
   elseif strcmp(modelType,'noBlood1')
        est_fb(i)=0;
 t_delayVp(i)=0;
 spillover(i)=0;
      [est_curves(:,i), tmpbldcurve]=mpi_sat_fwdModel(bldcurve, sat_bldcurve, delta_t, kwi(i), kwo(i), t_delay(i), est_fb(i), spillover(i));
   elseif strcmp(modelType,'full')
      [est_curves(:,i), tmpbldcurve]=mpi_sat_fwdModel(bldcurve, sat_bldcurve, delta_t, kwi(i), kwo(i), t_delay(i), est_fb(i), spillover(i));
   elseif strcmp(modelType,'fullModel2')
      [est_curves(:,i), tmpbldcurve]=mpi_sat_fwdModel2(bldcurve, sat_bldcurve, delta_t, kwi(i), kwo(i), t_delay(i), est_fb(i), spillover(i));
   elseif strcmp(modelType,'stdModel')
      [est_curves(:,i), tmpbldcurve]=mpi_fwdModel(bldcurve, delta_t, kwi(i), kwo(i), t_delay(i), est_fb(i));
    else
      [est_curves(:,i), tmpbldcurve]=mpi_sat_fwdModel(bldcurve, sat_bldcurve, delta_t, kwi(i), kwo(i), t_delay(i), est_fb(i), spillover(i));
   end

   est_curves(:,i)=est_curves((jitter+1):newRate:nTimes,i);   % downsampling ste

      err = (tisscurve(:,i) - est_curves(:,i));
      errsum(i) = errsum(i) + sum(err.*err)
      clf; hold on
      set(gca,'FontSize',16); set(gcf,'Color',[1 1 1]);
plot(ttimes,bldcurve,'r','LineWidth',3)
      plot(ttimes,tisscurve(:,i),'x','LineWidth',2.2)
      plot(ttimes,est_curves(1:nTimes,i),'k','LineWidth',3)
  %    plot(ttimes,err,'mo','LineWidth',1)
      disp('numRuns is ')
      numRuns=runtest(err)
      plot(ttimes,zeros(nTimes,1),'g','LineWidth',1)
      legend('Blood Input','Measured Uptake', 'Fit (Model)')
      plot(ttimes,est_curves(1:nTimes,i),':')
      xlabel('Time (sec)')
%      ylabel('Gd concentration or deltaSI')
      ylabel('deltaSI')
%      if studyNum>100
          pause
%      end
      

end
errsum
sum(errsum)


if studyNum<100   % only write out for simulating  (series  99x are the simulations)
      outfilename=strcat(outpath,'fitCurves.study',int2str(studyNum),'.slice',int2str(sliceNum),'.',int2str(numAzimuthalRegions), '.',int2str(numRadialRegions),'.AIF_',int2str(seriesNumAIF),'_',int2str(sliceNumAIF),'_', num2str(scaleAIF),'.model1.mat')
      deltaSIcurves=[bldcurve'; est_curves'];
      save(outfilename, 'deltaSIcurves');
%      save outfilename bldcurve -append;
end



return
