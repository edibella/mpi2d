function showfits(curves, studyNum, sliceNum, outpath, numAzimuthalRegions, numRadialRegions, clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp, seriesNumAIF, sliceNumAIF, scaleAIF, modelType, lastFrame, timeStampFile, timeStampFileAIF)

figure(22);
if ~exist('modelType')
    modelType='xxx'
 %   modelType='stdModel'
end
%if exist('timeStampFileAIF')  % I think assumes this always exists? haven't thought of will handle a call to this function without it
if strcmp(timeStampFileAIF,'')  % I think assumes this always exists? haven't thought of will handle a call to this function without it
   flagTimeStamps=1;
   if strcmp(timeStampFileAIF,'') & exist('timeStampFile')
      timeStampFileAIF=timeStampFile;
   end
end
flagTimeStamps=1
if strcmp(timeStampFile,'') 
   flagTimeStamps=0;
   clear timeStampFile; clear timeStampFileAIF;
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
nRegs=size(curves,1)-2;
% need to figure out good way to pass in saturated bld curve. 
tisscurves=curves(2:nRegs+1,:);
tisscurve=tisscurves';
sat_bldcurve=curves(nRegs+2,:)';

% from older, should do try, catch
%filename=strcat(outpath,'fitparams.study',int2str(studyNum),'.slice',int2str(sliceNum),'.',int2str(numAzimuthalRegions),'.',int2str(numRadialRegions),'.txt');


%flagTimeSamples=0;
%flagTimeSamples=1;
if exist('seriesNumAIF')
   filename=mpi_getFilename(sliceNum, studyNum, outpath, numAzimuthalRegions, numRadialRegions, clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp, seriesNumAIF, sliceNumAIF, scaleAIF, flagTimeStamps)
   % 9/29/05 Also read in saturated AIF for modeling:
   filenamePart1='deltaSIcurves';
   curvefilename=strcat(outpath,filenamePart1,'.study',int2str(studyNum),'.slice',int2str(sliceNum),'.mat')
   %curvefilename=strcat(outpath,filenamePart1,'.study',int2str(studyNum),'.slice',int2str(sliceNum),'.AIF_',int2str(seriesNumAIF),'_',int2str(sliceNumAIF),'_',num2str(scaleAIF),'.mat');
 
   %curvefilename=strcat(outpath,filenamePart1,'.study',int2str(studyNum),'.slice',int2str(sliceNum),'.AIF_',int2str(studyNum),'_',int2str(sliceNum),'_',int2str(numRadialRegions),'.mat')
   %curvefilename=strcat(outpath,filenamePart1,'.study',int2str(studyNum),'.slice',int2str(sliceNum),'.AIF_',int2str(seriesNumAIF),'_',int2str(sliceNumAIF),'_',int2str(scaleAIF),'.mat') % changed below to display fits when alternate AIF is used NP 10/24/06
   
   load(curvefilename);
   curves = eval(filenamePart1);
   %sat_bldcurve=curves(1,1:length(bldcurve))';
else
   filename=mpi_getFilename(sliceNum, studyNum, outpath, numAzimuthalRegions, numRadialRegions,  clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp);
end




Hct=0.45;
if length(modelType)>3
   if strcmp(modelType(1:4),'ferm')   % so minimum length of modelType is 4 letters
      filename=strcat(filename,'.',modelType,'.txt')
      [FLOW_0 Scalar_F ve est_fb t_delay T0 fval delta_t spillover]=readFitsFermi(filename)
      kwi=Scalar_F;
      kwo=ve;
   end
end
if strcmp(modelType,'jw')
   filename=strcat(filename,'.jw.txt')
   [flow ve ext t_delay TC fval delta_t]=readFitsJW(filename)
   kwi=flow;
   kwo=(1-Hct)*kwi./ve;
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
  % kwi=flow*0.5;  % if wrote out not Ktrans but Ktrans/0.5
    kwi=flow;
   kwo=(1-Hct)*kwi./ve;
elseif strcmp(modelType,'fullTrunc')
   filename=strcat(filename,'.fullTrunc.txt')
   [flow ve est_fb t_delay fval delta_t spillover]=readFits(filename)
   kwi=flow*0.5;
   kwo=(1-Hct)*kwi./ve;
elseif strcmp(modelType,'fullNoP')
   filename=strcat(filename,'.fullNoP.txt')
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
%else
%   [flow ve est_fb t_delay fval delta_t spillover]=readFits(filename)
%   kwi=flow*0.5;
 %  kwo=(1-Hct)*kwi./ve;
end

delta_t

tmp=0:nTimes-1;
ttimes = delta_t * tmp';
ttimes = tmp';


%  if strcmp(modelType,'fermiFullTrunc') || strcmp(modelType,'fullTrunc') || strcmp(modelType,'fermiTrunc')
  if strcmp(modelType,'noBlood1') || strcmp(modelType,'fullTrunc') || strcmp(modelType,'fermiTrunc')
      nTimes=lastFrame;
      tisscurve=tisscurve(1:lastFrame,:);
      bldcurve=bldcurve(1:lastFrame);
      sat_bldcurve=sat_bldcurve(1:lastFrame);
      tmp=0:nTimes-1;
      ttimes = delta_t * tmp';
      ttimes = tmp';
   end



%[timeDiff_AIF, bldcurve]=alignCurves(bldcurve, sat_bldcurve);

numSkip=1
numSkip=0
bldcurve=bldcurve(numSkip+1:nTimes);
nTimes=nTimes-numSkip;
sat_bldcurve=bldcurve;

   tmp=0:nTimes-1;
      ttimes = delta_t * tmp';
      %ttimes = tmp';


jitter=0; newRate=1;

errsum = zeros(1,nRegs);
for i=1:nRegs
    nn(:,i)=tisscurve(numSkip+1:nTimes+numSkip,i);
end
clear tisscurve;
tisscurve=nn;
clf; hold on
hh=gcf;
figure(77); clf; hold on; set(gcf,'position',[100 100 800 400])
figure(hh);
for i=1:nRegs
   disp('iReg is '),i
  
   
%   tisscurve(:,i)=nn;
 if length(modelType) > 3
   if strcmp(modelType(1:4),'ferm')
%      [est_curves(:,i), tmpbldcurve]=mpi_fermi_fwdModel(bldcurve, sat_bldcurve, delta_t, kwi(i), kwo(i),t_delay(i), est_fb(i));
      [est_curves(:,i), tmpbldcurve]=mpi_fermi_fwdModel(bldcurve, sat_bldcurve, delta_t, kwi(i), kwo(i), t_delay(i), T0(i), est_fb(i), spillover(i));
   end
 end
 if strcmp(modelType,'jw')
      [est_curves(:,i), tmpbldcurve]=mpi_jw_fwdModel(x,bldcurve, sat_bldcurve, tisscurve,delta_t,kwi(i), kwo(i), ext(i), t_delay(i), TC(i))
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
%    else
%      [est_curves(:,i), tmpbldcurve]=mpi_sat_fwdModel(bldcurve, sat_bldcurve, delta_t, kwi(i), kwo(i), t_delay(i), est_fb(i), spillover(i));
   end

   est_curves(:,i)=est_curves((jitter+1):newRate:nTimes,i);   % downsampling ste

      err = (tisscurve(:,i) - est_curves(:,i));
      errsum(i) = errsum(i) + sum(err.*err)
      perf=kwi(i)   % just to print out each time
      %clf; hold on
      set(gca,'FontSize',16); set(gcf,'Color',[1 1 1]);
%keyboard
plot(ttimes,bldcurve,'r','LineWidth',3)
      plot(ttimes,tisscurve(:,i),'x','LineWidth',2.2)
      plot(ttimes,est_curves(1:nTimes,i),'k','LineWidth',3)
  %    plot(ttimes,err,'mo','LineWidth',1)
   %   disp('numRuns is ')
      numRuns=runtest(err);  % runs of errors of same sign... or similar
      plot(ttimes,zeros(nTimes,1),'g','LineWidth',1)
      legend('Blood Input','Measured Uptake', 'Fit (Model)')
      plot(ttimes,est_curves(1:nTimes,i),':')
      xlabel('Time (sec)')
%      ylabel('Gd concentration or deltaSI')
      ylabel('deltaSI')
      hold on;
%      if studyNum>100
          pause(.1);
%      end
    
       fignum=gcf;
       firstTime=1;
    figure(77); hold on; set(gcf,'position',[100 100 800 400])
    ColOrd = get(gca,'ColorOrder');
    plot(ttimes,est_curves(1:nTimes,i),'Color',ColOrd(i,:),'LineWidth',2)
    plot(ttimes,tisscurve(:,i),'x--','Color',ColOrd(i,:),'LineWidth',1.3)
    %plot(ttimes,tisscurve(:,i),'Color',ColOrd(i,:),'linestyle','x--','LineWidth',1.3)
    
%       [AX,H1,H2] = plotyy(ttimes,bldcurve,ttimes,est_curves,'plot');
%       set(H1,'Color','r','linewidth',2.2)
%       set(H2,'Color','k','linewidth',2.2)
%       if firstTime
%           ylim=get(AX(2),'YLim');
%           firstTime=0;
%       else
%           set(AX(2),'YLim',ylim);
%       end
%       [AX,H1,H2] = plotyy(ttimes,bldcurve,ttimes,tisscurve(:,i),'plot');
%       set(H2,'Color','b','linestyle','o','linewidth',2.2)
      pause
figure(fignum)

end
errsum
sum(errsum)

% set(get(AX(1),'Ylabel'),'String','AIF') 
%       set(get(AX(2),'Ylabel'),'String','Tissue') 
% 
if studyNum<100   % only write out for simulating  (series  99x are the simulations)
      outfilename=strcat(outpath,'fitCurves.study',int2str(studyNum),'.slice',int2str(sliceNum),'.',int2str(numAzimuthalRegions), '.',int2str(numRadialRegions),'.AIF_',int2str(seriesNumAIF),'_',int2str(sliceNumAIF),'_', num2str(scaleAIF),'.',modelType,'.mat')
      deltaSIcurves=[bldcurve'; est_curves'];
      save(outfilename, 'deltaSIcurves');
%      save outfilename bldcurve -append;
end

% debugging 5/19/06
%i=1;
%%   if strcmp(modelType,'full')
%kwi(1)=7.514/2.0; ve(1)=0.216; t_delay(1)=0.355;
%   kwo(1)=(1-Hct)*kwi(1)./ve(1);
%
%load Output/deltaSIcurves.study998.slice3.AIF_998_3_1.mat
%%bldcurve=deltaSIcurves(1,:);
%%sat_bldcurve=bldcurve';
%est_curves(:,1)=deltaSIcurves(2,:)';
%
%      [test_est_curves(:,i), tmpbldcurve]=mpi_sat_fwdModel(bldcurve, sat_bldcurve, delta_t, kwi(i), kwo(i), t_delay(i), est_fb(i), spillover(i));
%      aa=test_est_curves(:,1) - est_curves(:,1);
%      max(aa)
%      min(aa)

%keyboard
%   end


return
