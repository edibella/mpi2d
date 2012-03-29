function fitparams = fit_model_2comp(seriesNumAIF, sliceNumAIF, scaleAIF, sliceNum,studyNum,outpath,filenamePart1,delta_t, numAzimuthalRegions, numRadialRegions, flagPixelwise, numSkip, curvesFromClusters, fixedVp, fixedDelay, modelType, framesToSelect, ranget, timeStampFile, timeStampFileAIF )

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
nRegs=size(curves,1)-1;
tisscurves=curves(2:nRegs+1,:);
tisscurve=tisscurves';
nTimes=length(bldcurve);

% 9/19/05 Also read in saturated AIF for modeling:
   curvefilename=strcat(outpath,filenamePart1,'.study',int2str(studyNum),'.slice',int2str(sliceNum),'.mat');
   load(curvefilename);
   curves = eval(filenamePart1);
   sat_bldcurve=curves(1,:)';
   sat_bldcurve=sat_bldcurve(1:length(bldcurve));

oldoptions=optimset('fmincon')  ;           
options=optimset('fmincon')   ;          
options=optimset(oldoptions,'TolFun',0.001,'TolX',0.001,'MaxIter',1e3) ;
options=optimset(oldoptions,'TolFun',0.0001,'TolX',0.0001,'MaxIter',1e4) ;

% 5/19/06
options=optimset(oldoptions,'TolFun',0.000001,'Display','off') ;

X0(1)=1.0; LB(1)=0.01*X0(1); UB(1)=20*X0(1);
X0(2)=2.5; LB(2)=0.01*X0(2); UB(2)=20*X0(2); % changed bounds from below on 10/23/06 NP, to try and constrain Ve to 0-100%
% X0(1)=0.5; LB(1)=0.01*X0(1); UB(1)=20*X0(1);
% X0(2)=0.5; LB(2)=0.01*X0(2); UB(2)=20*X0(2); % these initial guesses and bounds are for K1 and k2 (Ktrans and Kep) 

%X0(2)=0.44; LB(2)=0.2; UB(2)=1.8; % changed to directly estimate the Partition Coefficient assuming K1/k2=Ve/(1-Hct) (NP 10/10/06)
%%%% The initial guesses and bounds above match those in "fit_modelFull.m"
%%%% Also, the guesses below--for multiple independent fits--match these bounds here

% for this part, want to fit for time delay!! % fised 5/1/06
orig_fixedDelay=fixedDelay;
fixedDelay=99;

if fixedDelay==99
   X0(3)=-0.4; LB(3)=-15.0; UB(3)=0.8;  % changed 8/25/05 evrd to do t_delay
else
  if fixedVp==99
      X0(3)=0.05; LB(3)=0.0;  UB(3)=0.15;  % changed from .7UB 10/3/05
  end
end
if fixedVp==99 & fixedDelay==99
   X0(4)=0.05; LB(4)=0.0;  UB(4)=0.15;   % changed 9/13/04 to try
                   % changed from .7UB to .15 10/3/05 and start guess from 0
% tempted to force LB to .05
end
%X0(5)=0.0; LB(5)=0.0; UB(5)=0.4;   % spillover

A=zeros(1,length(X0));
B=0;
A=[];
B=[];
nonlcon=[];


% resample to uniform time points:
flagTimeStamps=0;
if exist('timeStampFile')
   flagTimeStamps=1;
   load(timeStampFile);   % assumes gives variable timeStamp
   if(size(timeStamp,1) == 1)
       %the timestamp was not saved in the format mpi2d expects
       timeStamp = timeStamp';
       save( timeStampFile,'timeStamp');
   end
    counter=1;
        for j=framesToSelect  %ranget    % EVRD 4/26/11, to crop to ranget. Assume subset of shifts. 
            if (j <=max(ranget)) && (j >=min(ranget))               
                selectedTimeStamps(counter,1) = timeStamp(j);
                counter=counter+1;
            end
        end
    
   bldcurve=interpTimeCurve(selectedTimeStamps, delta_t, bldcurve);
   nTimes=length(bldcurve);

% now do different one for Tissues:
   origtisscurve=tisscurve;
   clear tisscurve

   for ii=1:nRegs
       try
           tisscurve(:,ii)=interpTimeCurve(timeStampFile, delta_t, origtisscurve(:,ii));
       catch
           disp('Likely a problem with a tissue region having Nans. Take a look with 3.21 or 3.11, redraw ')
       end
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
estGlobalDelay=x(3);
fixedDelay=orig_fixedDelay;  % now put original fixedDelay back

if fixedDelay==0
  fixedDelay=estGlobalDelay;
end
%%%%%%%% CHANGED NP 071807 to force deltaT to be zero (also in sat_func_fit.m)
% fixedDelay=0;
% estGlobalDelay=0;

fvalues=99*10^9*ones(1,nRegs); clear xx;
% adding  4/7/06 to do multiple different starts to try and get lowest chi-sq.
nFitRuns=2;
for jj=1:nFitRuns  
   clear X0
   X0(1)=0.5; LB(1)=0.01*X0(1); UB(1)=20*X0(1);
   X0(2)=0.5; LB(2)=0.01*X0(2); UB(2)=20*X0(2);
   %X0(2)=0.44; LB(2)=0.2; UB(2)=1.8; % changed to directly estimate the Partition Coefficient (NP 10/10/06)


   X0(1)=0.5*jj; LB(1)=0.01*X0(1); UB(1)=20*X0(1);
   X0(2)=0.5*jj; LB(2)=0.01*X0(2); UB(2)=20*X0(2);
   %X0(2)=0.44; LB(2)=0.2; UB(2)=1.8; % changed to directly estimate the Partition Coefficient (NP 10/10/06)

% 5/18/06 debugging
%   X0(1)=7.514/2.0;
%   X0(2)=9.5664;

%%%% other realistic guesses NP071807 to result in a K1=0.6 ml/g/min and k2=1.65/(1-Hct) min^-1
% X0(1)=0.6;  % to test with known guesses
% X0(2)=1.65/(1-0.45);

   if jj==4
       X0(1)=2; LB(1)=0.2; UB(1)=6;
       X0(2)=2; LB(2)=0.21; UB(2)=10*X0(2);
       %X0(2)=0.44; LB(2)=0.2; UB(2)=1.8; % changed to directly estimate the Partition Coefficient (NP 10/10/06)
   end
   
if fixedDelay==99
%   X0(3)=-0.4; LB(3)=-15.0; UB(3)=0.8;  % changed 8/25/05 evrd to do t_delay
   X0(3)=-1; LB(3)=-15.0; UB(3)=0;  % changed 8/25/05 evrd to do t_delay
   if jj==1
      X0(3)=-0.4; LB(3)=-15.0; UB(3)=0.8;  % changed 8/25/05 evrd to do t_delay
   elseif jj==2
      X0(3)=-7; LB(3)=-15.0; UB(3)=0.9;  % changed 8/25/05 evrd to do t_delay
   elseif jj==3
      X0(3)=-2; LB(3)=-15.0; UB(3)=0.9;  % changed 8/25/05 evrd to do t_delay
   end
else
  if fixedVp==99
   X0(3)=0.05; LB(3)=0;  UB(3)=0.15; 
   if jj==1
      X0(3)=0.05; LB(3)=0;  UB(3)=0.15;
   elseif jj==2
      X0(3)=0; LB(3)=0; UB(3)=0.5;  % changed 8/25/05 evrd to do t_delay
   elseif jj==3
      X0(3)=0.15; LB(3)=0.02; UB(3)=0.2;  % changed 8/25/05 evrd to do t_delay
   end
  end
end
if fixedVp==99 & fixedDelay==99
   X0(4)=0.05; LB(4)=0.0;  UB(4)=0.15;  
   if jj==1
      X0(4)=0.05; LB(4)=0.0;  UB(4)=0.15;   % changed 9/13/04 to try
   elseif jj==2
      X0(4)=0; LB(4)=0; UB(4)=0.5;  % changed 8/25/05 evrd to do t_delay
   elseif jj==3
      X0(4)=0.15; LB(4)=0.02; UB(4)=0.2;  % changed 8/25/05 evrd to do t_delay
   end
end
if (strcmp(modelType,'full')==1)  || strcmp(modelType,'fullTrunc')==1   % will always have fixedVp=99
   if fixedDelay==99
      paramNum=5;
   else
      paramNum=4;
   end
   X0(paramNum)=0.0; LB(paramNum)=0.0; UB(paramNum)=0.0;   % spillover set to zero here by constraining bounds to zero explicitly..
   if jj==1
      X0(paramNum)=0.0; LB(paramNum)=0.0; UB(paramNum)=0.0;   % spillover
   elseif jj==2
      X0(paramNum)=0.0; LB(paramNum)=0.0; UB(paramNum)=0.0;   % spillover
   elseif jj==3
      X0(paramNum)=0.0; LB(paramNum)=0.0; UB(paramNum)=0.0;  % changed 8/25/05 evrd to do t_delay
   end
end


A=zeros(1,length(X0));
B=0; A=[]; B=[]; nonlcon=[];

disp(['fit number is ' num2str(jj)])

for ii=1:nRegs
   disp(['Fitting - doing iterations for modeltype:' modelType])
   warning('off', 'all')
   
   [x,fval,exitflag,output,lambda,grad,hessian]=fmincon('sat_func_fit',X0,A,B,A,B,LB,UB,nonlcon,options,bldcurve(numSkip+1:nTimes),sat_bldcurve(numSkip+1:nTimes),tisscurve(:,ii),delta_t, 1, fixedDelay, fixedVp, estGlobalDelay);  
% -- To get the least value of chi square -- %
    warning('on', 'all')
   if fval < fvalues(ii)
      fvalues(ii)=fval;
      fitNumForThisRegion(ii)=jj;   % index of region # and fit # (dont actually need) 
      xx(ii,:)= x;   % concat regions
   end
%      fvalues=[fvalues fval];
end
%keyboard
fitNumForThisRegion;
fvalues;
end % jj looop

for ii=1:nRegs
if estGlobalDelay~=99 && fixedDelay==99
     % generate new fval, without lambda=.33*errsum;
     estGlobalDelay=99;  % so wont' use lambda
     fvalues(ii)=sat_func_fit(xx(ii,:),bldcurve, sat_bldcurve, tisscurve,delta_t,1,fixedDelay, fixedVp, estGlobalDelay);
     
%         set(gca,'FontSize',16); set(gcf,'Color',[1 1 1]);
%       plot(ttimes,bldcurve,'r','LineWidth',3)
%       plot(ttimes,tisscurve(:,i),'x','LineWidth',2.2)
%       plot(ttimes,est_curves(1:nTimes,i),'k','LineWidth',3)
%       numRuns=runtest(err);  % runs of errors of same sign... or similar
%       legend('Blood Input','Measured Uptake', 'Fit (Model)')
%       plot(ttimes,est_curves(1:nTimes,i),':')
%       xlabel('Time (sec)')
% %      ylabel('Gd concentration or deltaSI')
%       ylabel('deltaSI')
%       hold on;
% %      if studyNum>100
%           pause(.1);
          
end
end


printFitParams(xx,fvalues, sliceNum, studyNum, outpath, delta_t, numAzimuthalRegions, numRadialRegions,  clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp, seriesNumAIF, sliceNumAIF, scaleAIF,flagTimeStamps, modelType)

fitNumForThisRegion;   % index of region # and fit # (dont actually need) 
%disp('pausing to see fits used for eadch regoin')

return
