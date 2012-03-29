function fitparams = fit_modelFermi(seriesNumAIF, sliceNumAIF, scaleAIF, sliceNum,studyNum,outpath,filenamePart1,delta_t, numAzimuthalRegions, numRadialRegions, flagPixelwise, numSkip, curvesFromClusters, fixedVp, fixedDelay, modelType, timeStampFile, timeStampFileAIF)

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

%  curvefilename=strcat(outpath,filenamePart1,'.study',int2str(studyNum),'.slice',int2str(sliceNum),'.AIF_',int2str(seriesNumAIF),'_',int2str(sliceNumAIF),'_',num2str(scaleAIF),'.mat');
   curvefilename=strcat(outpath,filenamePart1,'.study',int2str(studyNum),'.slice',int2str(sliceNum),'.AIF_',int2str(seriesNumAIF),'_',int2str(sliceNumAIF),'_',num2str(scaleAIF),'.mat');
   load(curvefilename);
   curves = eval(filenamePart1);   % may want to re-do this so not so tricky
% % %%%%% below was only used to create a mean single tis curve for display in model comparison paper   
% % tmp=mean(curves(2:end,:));curves=[curves(1,:);tmp];

bldcurve=curves(1,:)';
nRegs=size(curves,1)-1
tisscurves=curves(2:nRegs+1,:);
tisscurve=tisscurves';
nTimes=length(bldcurve);

% 9/19/05 Also read in saturated AIF for modeling:
%   curvefilename=strcat(outpath,filenamePart1,'.study',int2str(studyNum),'.slice',int2str(sliceNum),'.mat')
   curvefilename=strcat(outpath,filenamePart1,'.study',int2str(studyNum),'.slice',int2str(sliceNum),'.AIF_',int2str(seriesNumAIF),'_',int2str(sliceNumAIF),'_',num2str(scaleAIF),'.mat');
   load(curvefilename);
   curves = eval(filenamePart1);
% %    %%%%% below was only used to create a mean single tis curve for display in model comparison paper   
% % tmp=mean(curves(2:end,:));curves=[curves(1,:);tmp];

   sat_bldcurve=curves(1,:)';
   sat_bldcurve=sat_bldcurve(1:length(bldcurve));

oldoptions=optimset('fmincon')             
options=optimset('fmincon')             
options=optimset(oldoptions,'TolFun',0.001,'TolX',0.001,'MaxIter',1e3) ;
options=optimset(oldoptions,'TolFun',0.0001,'TolX',0.0001,'MaxIter',1e4) ;

X0(1)=0.5; LB(1)=0.01*X0(1); UB(1)=20*X0(1); 
X0(2)=0.5; LB(2)=0.01*X0(2); UB(2)=20*X0(2);
X0(2)=0.5; LB(2)=0.01*X0(2); UB(2)=80*X0(2); % changed to much larger kwo values when data is truncated NP 08/22/08
if fixedDelay==99
   X0(3)=-0.4; LB(3)=-15.0; UB(3)=0.8;  % changed 8/25/05 evrd to do t_delay
   X0(3)=-0.4; LB(3)=-15.0; UB(3)=15;  % changed TO values when data is truncated NP 08/22/08
else
  if fixedVp==99
      X0(3)=0.05; LB(3)=0.0;  UB(3)=0.15;  % changed from .7UB 10/3/05
      X0(3)=0.05; LB(3)=0.0;  UB(3)=15;  % changed TO values when data is truncated NP 08/22/08
  end
end
if fixedVp==99 & fixedDelay==99
   X0(4)=0.05; LB(4)=0.0;  UB(4)=0.15;   % changed 9/13/04 to try
                   % changed from .7UB to .15 10/3/05 and start guess from 0
% tempted to force LB to .05
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

sat_bldcurve=sat_bldcurve(1:length(bldcurve));
[estGlobalDelay, bldcurve]=alignCurves(bldcurve, sat_bldcurve);  % shifts bldcurve to match sat_bldcurve, integer shifts. So can use single delay in fits.
% find delay by fitting to mean tissue curve:
estGlobalDelay=99;

if (strcmp(modelType, 'fermiConstrDelay')==1)
%   [x,fval,exitflag,output,lambda,grad,hessian]=fmincon('sat_func_fitNoBlood_orig',X0,A,B,A,B,LB,UB,nonlcon,options,bldcurve(numSkip+1:nTimes),sat_bldcurve(numSkip+1:nTimes),mean(tisscurve,2),delta_t, 1, fixedDelay, fixedVp);  
%   [x,fval,exitflag,output,lambda,grad,hessian]=fmincon('fermi_func_fit',X0,A,B,A,B,LB,UB,nonlcon,options,bldcurve(numSkip+1:nTimes),sat_bldcurve(numSkip+1:nTimes),mean(tisscurve,2),delta_t, 1, fixedDelay, fixedVp);  
%   estGlobalDelay=x(3)
end

clear X0;
% range for Fermi function 9/21/05
   X0(1)=0.5; LB(1)=0.01*X0(1); UB(1)=20*X0(1); 
   X0(2)=4.0; LB(2)=0.001*X0(2); UB(2)=10*X0(2); % changed kwo values when data is truncated NP 08/22/08
   X0(3)=0.0; LB(3)=0.0; UB(3)=25.0  ;  % T0 term 
   X0(3)=2.0; LB(3)=1.0; UB(3)=25.0  ;  % T0 term 
   X0(3)=2.0; LB(3)=0.0; UB(3)=15.0  ;  % T0 term %%% this range varies when the data is truncated for first-pass analysis only NP 08/22/08
  

if fixedDelay==99
   X0(4)=-0.4; LB(4)=-15.0; UB(4)=0.8;  
else
  if fixedVp==99
      X0(4)=0.05; LB(4)=0.0;  UB(4)=0.15;  % changed from .7UB 10/3/05
  end
end
if fixedVp==99 & fixedDelay==99
   X0(5)=0.05; LB(5)=0.0;  UB(5)=0.15; 
end
if (strcmp(modelType, 'fermiFull')==1 || strcmp(modelType, 'fermiFullTrunc')==1)  % so can have more than fermiFull
    X0(6)=0.0; LB(6)=0.0; UB(6)=0.4;   % spillover
end

A=zeros(1,length(X0));
B=0; A=[]; B=[]; nonlcon=[];

if (strcmp(modelType, 'fermiConstrDelay')==1)
%   [x,fval,exitflag,output,lambda,grad,hessian]=fmincon('sat_func_fitNoBlood_orig',X0,A,B,A,B,LB,UB,nonlcon,options,bldcurve(numSkip+1:nTimes),sat_bldcurve(numSkip+1:nTimes),mean(tisscurve,2),delta_t, 1, fixedDelay, fixedVp);  
   [x,fval,exitflag,output,lambda,grad,hessian]=fmincon('fermi_func_fit',X0,A,B,A,B,LB,UB,nonlcon,options,bldcurve(numSkip+1:nTimes),sat_bldcurve(numSkip+1:nTimes),mean(tisscurve,2),delta_t, 1, fixedDelay, fixedVp, estGlobalDelay);  
   estGlobalDelay=x(4)
end


fvalues=99999*ones(1,nRegs); clear xx;
% adding  4/7/06 to do multiple different starts to try and get lowest chi-sq.
nFitRuns=4;
for jj=1:nFitRuns  
   clear X0
   X0(1)=0.5; LB(1)=0.01*X0(1); UB(1)=20*X0(1); 
   X0(2)=4.0; LB(2)=0.001*X0(2); UB(2)=10*X0(2); % changed kwo value when the data is trucated at first-pass 08/22/08 NP
   X0(3)=12.0; LB(3)=1.0; UB(3)=25.0  ;  % T0 term 
   X0(3)=1.0; LB(3)=0.5; UB(3)=15.0  ;  % changed T0 term when the data is trucated at first-pass 08/22/08 NP 
   if jj==2
      X0(2)=8.0; LB(2)=0.1*X0(2); UB(2)=2*X0(2); 
      X0(2)=8.0; LB(2)=0.01*X0(2); UB(2)=5*X0(2); % changed kwo value when the data is trucated at first-pass 08/22/08 NP
      X0(3)=1.0; LB(3)=1.0; UB(3)=25.0  ;  % T0 term 
      X0(3)=1.0; LB(3)=0.0; UB(3)=15.0  ;  % changed to a much larger T0 term when the data is trucated at first-pass 08/22/08 NP 
  end
   if jj==3
      X0(2)=8.0; LB(2)=0.1*X0(2); UB(2)=2*X0(2); 
      X0(2)=8.0; LB(2)=0.01*X0(2); UB(2)=5*X0(2); % changed kwo value when the data is trucated at first-pass 08/22/08 NP
      X0(3)=4.0; LB(3)=0.5; UB(3)=25.0  ;  % T0 term 
      X0(3)=1.0; LB(3)=0.05; UB(3)=15.0  ;  % changed T0 term when the data is trucated at first-pass 08/22/08 NP 
   end
   if jj==4
      X0(2)=2.0; LB(2)=0.1*X0(2); UB(2)=8*X0(2); 
      X0(2)=2.0; LB(2)=0.1*X0(2); UB(2)=40*X0(2); % changed to a much larger value when the data is trucated at first-pass 08/22/08 NP
      X0(3)=1.0; LB(3)=0; UB(3)=12.0  ;  % T0 term 
      X0(3)=1.0; LB(3)=0; UB(3)=15.0  ;  % changed to a much larger T0 term when the data is trucated at first-pass 08/22/08 NP 
   end
   
   if fixedDelay==99
      X0(4)=-0.4; LB(4)=-15.0; UB(4)=0.0  %0.8;
      if jj==2
         X0(4)=0; LB(4)=-15.0; UB(4)=-0.5   %0.9;
      elseif jj==3
         X0(4)=-2; LB(4)=-15.0; UB(4)=-0.1   %0.9;
      end
   else
     if fixedVp==99
          X0(4)=0.05; LB(4)=0;  UB(4)=0.15;
        if jj==1
          X0(4)=0.05; LB(4)=0;  UB(4)=0.15;  % changed from .7UB 10/3/05
        elseif jj==2
          X0(4)=0; LB(4)=0; UB(4)=0.5;  % changed 8/25/05 evrd to do t_delay
        elseif jj==3
          X0(4)=0.15; LB(4)=0.02; UB(4)=0.2;  % changed 8/25/05 evrd to do t_delay
        end
     end
   end
if fixedVp==99 & fixedDelay==99
   X0(5)=0.05; LB(5)=0.0;  UB(5)=0.15; 
   if jj==2
      X0(5)=0; LB(5)=0; UB(5)=0.5;  % changed 8/25/05 evrd to do t_delay
   elseif jj==3
      X0(5)=0.15; LB(5)=0.02; UB(5)=0.2;  % changed 8/25/05 evrd to do t_delay
   end
end
if (strcmp(modelType, 'fermiFull')==1 || strcmp(modelType, 'fermiFullTrunc')==1)  % so can have more than fermiFull
 % need to fix this 4/26/06! will be X0(5) if delay is fixed
 % also fix above that jj==4, undefined! have been getting best results with jj=4!
   if fixedDelay==99
     paramNum=6;
   else
     paramNum=5;
   end
    X0(paramNum)=0.0; LB(paramNum)=0.0; UB(paramNum)=0.4;   % spillover
    if jj==2
      X0(paramNum)=0.05; LB(paramNum)=0.05; UB(paramNum)=0.3;   % spillover
    elseif jj==3
      X0(paramNum)=0.15; LB(paramNum)=0; UB(paramNum)=0.2;
    end
end

   if jj==4
       X0(1)=2; LB(1)=0.2; UB(1)=6;
       X0(2)=2; LB(2)=0.21; UB(2)=10*X0(2);
       X0(2)=2; LB(2)=0.021; UB(2)=20*X0(2);  % changed kwo range when the data is trucated at first-pass 08/22/08 NP
    end


A=zeros(1,length(X0));
B=0; A=[]; B=[]; nonlcon=[];
if(jj==1)
    % do on average region to get time delay to fix (if t_delay=0), or
    % constrain
    X0tmp=X0; LBtmp=LB; UBtmp=UB;
    X0tmp(4)=-0.4; LBtmp(4)=-15.0; UBtmp(4)=-0.2  %0.8;
    if fixedVp==99
       X0tmp(5)=0.05; LBtmp(5)=0.0;  UBtmp(5)=0.15; 
    end
    %%% below fixedDelay is always set to 99 in order to estimate the global time delay which can then be optionally used in the regional curve fits.
    [x,fval,exitflag,output,lambda,grad,hessian]=fmincon('fermi_func_fit',X0tmp,A,B,A,B,LBtmp,UBtmp,nonlcon,options,bldcurve(numSkip+1:nTimes),sat_bldcurve(numSkip+1:nTimes),mean(tisscurve,2),delta_t, 1, 99, fixedVp, estGlobalDelay);
    estGlobalDelay= x(4)
 
    if fixedDelay==0
        fixedDelay=estGlobalDelay
    end
end
  
disp('fit number is '),jj


for ii=1:nRegs
   disp('Fitting - doing iterations for modeltype'), modelType
   ii
   warning('off', 'all')
   
   [x,fval,exitflag,output,lambda,grad,hessian]=fmincon('fermi_func_fit',X0,A,B,A,B,LB,UB,nonlcon,options,bldcurve(numSkip+1:nTimes),sat_bldcurve(numSkip+1:nTimes),tisscurve(:,ii),delta_t, 1, fixedDelay, fixedVp, estGlobalDelay);  
    warning('on','all')
   

   if fval < fvalues(ii)
      fvalues(ii)=fval;
      fitNumForThisRegion(ii)=jj;   % index of region # and fit # (dont actually need) 
      xx(ii,:)= x;   % concat regions
   end
%      fvalues=[fvalues fval];
end

fitNumForThisRegion
fvalues
end % jj looop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if estGlobalDelay~=99 && fixedDelay==99
  for ii=1:nRegs
     % generate new fval, without lambda=.33*errsum;
     estGlobalDelay=99;  % so wont' use lambda
     fvalues(ii)=fermi_func_fit(xx(ii,:),bldcurve, sat_bldcurve, tisscurve(:,ii),delta_t,1,fixedDelay, fixedVp, estGlobalDelay);
%  matlab function to be called by fmincon
 end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ii=1:nRegs
   Scalar_F(ii) = xx(ii,1);
   kwo(ii)= xx(ii,2);
   Tzero(ii)=xx(ii,3);
   est_fb(ii)=fixedVp;

   hr=60/delta_t;
   FLOW_0(ii) = (hr*(delta_t/60)*xx(ii,1)./(exp(xx(ii,2)*(0-xx(ii,3)/60.0))+1) ); %% added 09/10/08 NP to find flow as hFermi(0) not simply the max[hFermi(t)]
   
   if fixedDelay==0
     t_delay(ii)=estGlobalDelay  
   else
     t_delay(ii)=fixedDelay ; 
   end
   if fixedDelay==99
      t_delay(ii) = xx(ii,4);
   else
     if fixedVp==99
        est_fb(ii)= xx(ii,4);
     end
   end
   if fixedVp==99 & fixedDelay==99
      est_fb(ii)= xx(ii,5);
   end
   if (strcmp(modelType,'fermiFull')==1 || strcmp(modelType, 'fermiFullTrunc')==1)
      if fixedDelay==99
        spillover(ii)=xx(ii,6);
        spillover(ii)=0;  % will force to zero here when using truncated data NP 08/22/08
  
      else
        spillover(ii)=xx(ii,5);
        spillover(ii)=0;  % will force to zero here when using truncated data NP 08/22/08
      end
   else
       spillover(ii)=0;
     end
% hacked in 5/16/06, also fixed in fermi_func_fit.m
   est_fb(ii)=0;
end


if fixedDelay~=99
     fixedDelay=0;
end
% above is ugly hack for writing out filename 4/26/06

if exist('seriesNumAIF')
   outfilename=mpi_getFilename(sliceNum, studyNum, outpath, numAzimuthalRegions, numRadialRegions,  clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp, seriesNumAIF, sliceNumAIF, scaleAIF, flagTimeStamps)
else
   outfilename=mpi_getFilename(sliceNum, studyNum, outpath, numAzimuthalRegions, numRadialRegions,  clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp)
end

outfilename=strcat(outfilename,'.',modelType,'.txt')
%outfilename=strcat(outfilename,'.fermi.constrDelay.txt')

fid=fopen(outfilename,'w');
index=1:length(FLOW_0);
fprintf(fid,'Created %s ,  delta_t=%f\n',date, delta_t);

out=[index; FLOW_0]; %%  added 09/10/08 NP ... this is the actual flow estimate, hFermi(t=0), scaled by the avg. hr.
fprintf(fid,'FLOW_0  %d    %6.3f\n',out); 
fprintf(fid,'\n');

fprintf(fid,'flowMean and Std  %6.2f %6.2f , coeff var= %6.2f\n',mean(FLOW_0), std(FLOW_0), std(FLOW_0)/mean(FLOW_0));
fprintf(fid,'\n');


out=[index; Scalar_F];  %% this is simply the scalar value of hFermi (ie. max[hFermi(t)]) and is NOT equal to flow unless the max of hFermi(t) occurs at time t=0
fprintf(fid,'Scalar_F  %d    %6.3f\n',out); 
fprintf(fid,'\n');

Hct=0.45; 
out=[index; kwo];  
fprintf(fid,'kwo    %d     %6.3f\n',out);
fprintf(fid,'\n');

out=[index; Tzero]; 
fprintf(fid,'Tzero    %d     %6.3f\n',out);
fprintf(fid,'\n');

out=[index; est_fb];  
fprintf(fid,'est_fb   %d     %6.3f\n',out);
fprintf(fid,'\n');

out=[index; t_delay]; 
fprintf(fid,'t0    %d     %6.3f\n',out);
fprintf(fid,'\n');

out=[index; spillover]; 
fprintf(fid,'spillover    %d     %6.3f\n',out);
fprintf(fid,'\n');

out=[index; fvalues]; 
fprintf(fid,'fval  %d     %6.3f\n',out);
fclose(fid);

return;
