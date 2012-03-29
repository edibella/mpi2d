function fitparams = fit_modelFull(seriesNumAIF, sliceNumAIF, scaleAIF, sliceNum,studyNum,outpath,filenamePart1,delta_t, numAzimuthalRegions, numRadialRegions, flagPixelwise, numSkip, curvesFromClusters, fixedVp, fixedDelay, timeStampFile, timeStampFileAIF)

fitparams=1;
if nargin==0,
    error('sliceNum argument needed for this file');
end
if ~exist('numSkip'),   % if manually skipping an intital bump desired
   numSkip=0;
end
if exist('timeStampFileAIF')
   if strcmp(timeStampFileAIF,'') & exist('timeStampFile')  
      timeStampFileAIF=timeStampFile;
   end
end

useIntegralLinearFit=0;
clusterFlag=0;
   curvefilename=strcat(outpath,filenamePart1,'.study',int2str(studyNum),'.slice',int2str(sliceNum),'.AIF_',int2str(seriesNumAIF),'_',int2str(sliceNumAIF),'_',num2str(scaleAIF),'.mat');
   load(curvefilename);
   curves = eval(filenamePart1);   % may want to re-do this so not so tricky

bldcurve=curves(1,:)';
nRegs=size(curves,1)-1
tisscurves=curves(2:nRegs+1,:);
tisscurve=tisscurves';
nTimes=length(bldcurve);

% 9/19/05 Also read in saturated AIF for modeling:
   curvefilename=strcat(outpath,filenamePart1,'.study',int2str(studyNum),'.slice',int2str(sliceNum),'.mat')
   load(curvefilename);
   curves = eval(filenamePart1);
   sat_bldcurve=curves(1,:)';

oldoptions=optimset('fmincon')             
options=optimset('fmincon')             
options=optimset(oldoptions,'TolFun',0.001,'TolX',0.001,'MaxIter',1e3) ;
options=optimset(oldoptions,'TolFun',0.0001,'TolX',0.0001,'MaxIter',1e4) ;

X0(1)=0.5; LB(1)=0.01*X0(1); UB(1)=20*X0(1); 
X0(2)=0.5; LB(2)=0.01*X0(2); UB(2)=20*X0(2); % these initial guesses and bounds are for kwo
%X0(2)=0.44; LB(2)=0.2; UB(2)=1.8; % changed to directly estimate the Partition Coefficient (NP 10/10/06)
%%%% The initial guesses and bounds above match those in "fit_model_2comp.m"

if fixedDelay==99
   X0(3)=-0.4; LB(3)=-15.0; UB(3)=0.8;  % changed 8/25/05 evrd to do t_delay
else
  if fixedVp==99
      X0(3)=0.05; LB(3)=0.0;  UB(3)=0.15;  % changed from .7UB 10/3/05
  end
end
if fixedVp==99 & fixedDelay==99
   X0(4)=0.05; LB(4)=0.0;  UB(4)=0.15;   % changed 9/13/04 to try
                      % changed from .7UB to .15 10/3/05 and start guess from 0 tempted to force LB to .05
end
X0(5)=0.0; LB(5)=0.0; UB(5)=0.4;   % spillover

A=zeros(1,length(X0));
B=0;
A=[];
B=[];
nonlcon=[];


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


[estGlobalDelay, bldcurve]=alignCurves(bldcurve, sat_bldcurve);  % shifts bldcurve to match sat_bldcurve, integer shifts. So can use single delay in fits.
%[estGlobalDelay, shiftedcurve]=alignCurves(bldcurve, mean(tisscurve,2));  % need to check is doing mean of 2nd
% find delay by fitting to mean tissue curve:

[x,fval,exitflag,output,lambda,grad,hessian]=fmincon('sat_func_fitNoBlood_orig',X0,A,B,A,B,LB,UB,nonlcon,options,bldcurve(numSkip+1:nTimes),sat_bldcurve(numSkip+1:nTimes),mean(tisscurve,2),delta_t, 1, fixedDelay, fixedVp);  
estGlobalDelay=x(3)



modelType='full'
fvalues=[]; clear xx;
for ii=1:nRegs
   disp('Fitting - doing iterations for noblood1 model')
   [x,fval,exitflag,output,lambda,grad,hessian]=fmincon('sat_func_fit',X0,A,B,A,B,LB,UB,nonlcon,options,bldcurve(numSkip+1:nTimes),sat_bldcurve(numSkip+1:nTimes),tisscurve(:,ii),delta_t, 1, fixedDelay, fixedVp, estGlobalDelay);  
   xx(ii,:)= x;   % concat regions
   fvalues=[fvalues fval];
end
printFitParams(xx,fvalues, sliceNum, studyNum, outpath, delta_t, numAzimuthalRegions, numRadialRegions,  clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp, seriesNumAIF, sliceNumAIF, scaleAIF,flagTimeStamps, modelType)


return
