function fitparams = fit_modelNoBlood1(seriesNumAIF, sliceNumAIF, scaleAIF, sliceNum,studyNum,outpath,filenamePart1,delta_t, numAzimuthalRegions, numRadialRegions, flagPixelwise, numSkip, curvesFromClusters, fixedVp, fixedDelay, timeStampFile, timeStampFileAIF)

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
X0(2)=0.5; LB(2)=0.01*X0(2); UB(2)=20*X0(2); 
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



modelType='noBlood1'
clear X0
X0(1)=0.5; LB(1)=0.01*X0(1); UB(1)=20*X0(1); 
X0(2)=0.5; LB(2)=0.01*X0(2); UB(2)=20*X0(2); 
if fixedDelay==99
   X0(3)=-0.4; LB(3)=-15.0; UB(3)=0.8;  % changed 8/25/05 evrd to do t_delay
else
  if fixedVp==99
      X0(3)=0.05; LB(3)=0.0;  UB(3)=0.15;  % changed from .7UB 10/3/05
  end
end
if fixedVp==99 & fixedDelay==99
   X0(4)=0.05; LB(4)=0.0;  UB(4)=0.15;   % changed 9/13/04 to try
end

A=zeros(1,length(X0));
B=0;
A=[];
B=[];
nonlcon=[];
fvalues=[]; clear xx;
for ii=1:nRegs
   disp('Fitting - doing iterations for noblood1 model')
   [x,fval,exitflag,output,lambda,grad,hessian]=fmincon('sat_func_fitNoBlood',X0,A,B,A,B,LB,UB,nonlcon,options,bldcurve(numSkip+1:nTimes),sat_bldcurve(numSkip+1:nTimes),tisscurve(:,ii),delta_t, 1, fixedDelay, fixedVp, estGlobalDelay);  
   xx(ii,:)= x;   % concat regions
   fvalues=[fvalues fval];
end
printFitParams(xx,fvalues, sliceNum, studyNum, outpath, delta_t, numAzimuthalRegions, numRadialRegions,  clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp, seriesNumAIF, sliceNumAIF, scaleAIF,flagTimeStamps, modelType)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   X0(1)=0.5; LB(1)=0.01*X0(1); UB(1)=20*X0(1); 
   X0(2)=0.5; LB(2)=0.01*X0(2); UB(2)=20*X0(2); 
if fixedDelay==99
   X0(3)=0.0; LB(3)=-10.0; UB(3)=0.8;  % all this section 10/24/05
else
  if fixedVp==99
      X0(3)=0.04; LB(3)=0.01;  UB(3)=0.15; 
  end
end
if fixedVp==99 & fixedDelay==99
   X0(4)=0.04; LB(4)=0.01;  UB(4)=0.15;
end
X0(5)=0.03; LB(5)=0.0; UB(5)=0.5;   % spillover

% temp 11/18/05
X0(1)=0.5; LB(1)=0.01*X0(1); UB(1)=20*X0(1);
X0(2)=0.5; LB(2)=0.01*X0(2); UB(2)=20*X0(2);  
X0(5)=0.03; LB(5)=0.0; UB(5)=0.5;   % spillover
X0(5)=0.0; LB(5)=0.0; UB(5)=0.5;   % spillover

modelType='funcfit'
fvalues=[]; clear xx;
for ii=1:nRegs
   [x,fval,exitflag,output,lambda,grad,hessian]=fmincon('func_fit',X0,A,B,A,B,LB,UB,nonlcon,options,bldcurve(numSkip+1:nTimes),tisscurve(:,ii),delta_t, 1, fixedDelay, fixedVp);  
   xx(ii,:)= x;   % concat regions
   fvalues=[fvalues fval];
end

printFitParams(xx,fvalues, sliceNum, studyNum, outpath, delta_t, numAzimuthalRegions, numRadialRegions,  clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp, seriesNumAIF, sliceNumAIF, scaleAIF,flagTimeStamps, modelType)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


modelType='fermi';   % only put in fermi now, need to make how writing out a function of this.

% range for Fermi function 9/21/05
%   X0(1)=0.5; LB(1)=0.01*X0(1); UB(1)=20*X0(1); 
   X0(2)=4.0; LB(2)=0.001*X0(2); UB(2)=100*X0(2); 
   X0(4)=0.0; LB(4)=0.0;  UB(4)=25.0;

% new starting guesses Nov 11 2005
clear X0
   X0(1)=0.5; LB(1)=0.01*X0(1); UB(1)=20*X0(1); 
   X0(2)=4.0; LB(2)=0.001*X0(2); UB(2)=100*X0(2); 
   X0(3)=2.0; LB(3)=0.0;  UB(3)=25.0;
% fix time delay
fixedDelayOrig=fixedDelay;
fixedDelay=estGlobalDelay;

for ii=1:nRegs
   disp('Fitting - doing iterations for Fermi model!')
   [x,fval,exitflag,output,lambda,grad,hessian]=fmincon('fermi_func_fit',X0,A,B,A,B,LB,UB,nonlcon,options,bldcurve(numSkip+1:nTimes),tisscurve(:,ii),delta_t, 1, fixedDelay, fixedVp);  
   flow(ii) = x(1);
   kwo(ii)= x(2);
   est_fb(ii)=fixedVp;
   if fixedDelay==99
      t_delay(ii) = x(3);
   else
     if fixedVp==99
        est_fb(ii)= x(3);
     end
     t_delay(ii)=fixedDelay;
   end
   if fixedVp==99 & fixedDelay==99
      est_fb(ii)= x(4);
   end

   Tzero=est_fb;
   fvalues(ii)=fval;
end

fixedDelay=fixedDelayOrig;  % set it back

if exist('seriesNumAIF')
   outfilename=mpi_getFilename(sliceNum, studyNum, outpath, numAzimuthalRegions, numRadialRegions,  clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp, seriesNumAIF, sliceNumAIF, scaleAIF, flagTimeStamps)
else
   outfilename=mpi_getFilename(sliceNum, studyNum, outpath, numAzimuthalRegions, numRadialRegions,  clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp)
end

outfilename=strcat(outfilename,'.fermi2.txt')

fid=fopen(outfilename,'w');
index=1:length(flow);
fprintf(fid,'Created %s ,  delta_t=%f\n',date, delta_t);
out=[index; flow]; 

fprintf(fid,'flow  %d    %6.3f\n',out); 
fprintf(fid,'\n');

fprintf(fid,'flowMean and Std  %6.2f %6.2f , coeff var= %6.2f\n',mean(flow), std(flow), std(flow)/mean(flow));
fprintf(fid,'\n');

Hct=0.45; 
out=[index; kwo];  
fprintf(fid,'kwo    %d     %6.3f\n',out);
fprintf(fid,'\n');

out=[index; Tzero]; 
fprintf(fid,'Tzero    %d     %6.3f\n',out);
fprintf(fid,'\n');

out=[index; t_delay]; 
fprintf(fid,'t0    %d     %6.3f\n',out);
fprintf(fid,'\n');

out=[index; fvalues]; 
fprintf(fid,'fval  %d     %6.3f\n',out);
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% now do multiple regions with single delay fit:

globalVp=0;
for ii=1:nRegs
   X0(ii)=0.5; LB(ii)=0.01*X0(ii); UB(ii)=20*X0(ii);
   X0(nRegs+ii)=0.5; LB(nRegs+ii)=0.01*X0(nRegs+ii); UB(nRegs+ii)=20*X0(nRegs+ii);
end
if globalVp
   X0(2*nRegs+1)=0.0; LB(2*nRegs+1)=0.0; UB(2*nRegs+1)=0.4;
   X0(2*nRegs+2)=-0.4; LB(2*nRegs+2)=-15.0; UB(2*nRegs+2)=0.8;
else
   X0(2*nRegs+ii)=0.05; LB(2*nRegs+ii)=0.0;  UB(2*nRegs+ii)=0.2;  % changed Nov 11 2005
   X0(3*nRegs+1)=-0.4; LB(3*nRegs+1)=-15.0; UB(3*nRegs+1)=0.8;
end

   disp('Fitting - doing iterations on all regions at once')
% 8/10/05

modelType='globaldelay2';
   [x,fval,exitflag,output,lambda,grad,hessian]=fmincon('funcMultipleRegions_fit',X0,A,B,A,B,LB,UB,nonlcon,options,bldcurve(numSkip+1:nTimes),tisscurve,delta_t, nRegs);  % last 3 are params passed to func(x,p1,p2,p3)

%  covMatrix=inv(hessian);
%  varKwi(ii)=covMatrix(1,1);

kwi= x(1:nRegs);
kwo= x(nRegs+1:2*nRegs);
est_fb= x(2*nRegs+1:3*nRegs);  % est_fb
est_fb=zeros(1,nRegs);    % changed Nov 11, 2005 - not fitting in funcMult*
t_delay = x(3*nRegs+1);
t_delay=t_delay*ones(1,nRegs);
% want all fval to be zero except last region 
for ii=1:nRegs
   fvalues(ii)=0;
end
fvalues(nRegs)=fval(1);

if exist('seriesNumAIF')
   outfilename=mpi_getFilename(sliceNum, studyNum, outpath, numAzimuthalRegions, numRadialRegions,  clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp, seriesNumAIF, sliceNumAIF, scaleAIF, flagTimeStamps)
else
   outfilename=mpi_getFilename(sliceNum, studyNum, outpath, numAzimuthalRegions, numRadialRegions,  clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp)
end
outfilename=strcat(outfilename,'.globaldelay2.txt')

fid=fopen(outfilename,'w');
flow=kwi/0.5;
index=1:length(flow);
fprintf(fid,'Created %s ,  delta_t=%f\n',date, delta_t);
out=[index; flow]; 

fprintf(fid,'flow  %d    %6.3f\n',out); 
fprintf(fid,'\n');

fprintf(fid,'flowMean and Std  %6.2f %6.2f , coeff var= %6.2f\n',mean(flow), std(flow), std(flow)/mean(flow));
fprintf(fid,'\n');
Hct=0.45; 
out=[index; (1-Hct)*kwi./kwo];   % becuase kwi for arterial
fprintf(fid,'ve    %d     %6.3f\n',out);
fprintf(fid,'\n');
out=[index; est_fb]; 
fprintf(fid,'est_fb    %d     %6.3f\n',out);
fprintf(fid,'\n');
out=[index; t_delay]; 
fprintf(fid,'t0    %d     %6.3f\n',out);
fprintf(fid,'\n');
out=[index; fvalues]; 
fprintf(fid,'fval  %d     %6.3f\n',out);
fclose(fid);


%%%%% now do multiple regions with single delay fit AND single Vp term:

globalVp=1;
for ii=1:nRegs
   X0(ii)=0.5; LB(ii)=0.01*X0(ii); UB(ii)=20*X0(ii);
   X0(nRegs+ii)=0.5; LB(nRegs+ii)=0.01*X0(nRegs+ii); UB(nRegs+ii)=20*X0(nRegs+ii);
end
if globalVp
   X0(2*nRegs+1)=0.05; LB(2*nRegs+1)=0.02; UB(2*nRegs+1)=0.2;  % changed Nov 11 2005
   X0(2*nRegs+2)=-0.4; LB(2*nRegs+2)=-15.0; UB(2*nRegs+2)=0.8;
else
   X0(2*nRegs+ii)=0.0; LB(2*nRegs+ii)=0.0;  UB(2*nRegs+ii)=0.4;
   X0(3*nRegs+1)=-0.4; LB(3*nRegs+1)=-15.0; UB(3*nRegs+1)=0.8;
end

   disp('Fitting - doing iterations on all regions at once')
% 8/10/05
   [x,fval,exitflag,output,lambda,grad,hessian]=fmincon('funcMultipleRegions_fit',X0,A,B,A,B,LB,UB,nonlcon,options,bldcurve(numSkip+1:nTimes),tisscurve,delta_t, nRegs);  % figures out if fit for all Vp based on number of uknowns in X0

%  covMatrix=inv(hessian);
%  varKwi(ii)=covMatrix(1,1);

kwi= x(1:nRegs);
kwo= x(nRegs+1:2*nRegs);
if globalVp
   est_fb= x(2*nRegs+1);
   est_fb=est_fb*ones(1,nRegs);
   t_delay = x(2*nRegs+2);
   t_delay=t_delay*ones(1,nRegs);
else
   est_fb= x(2*nRegs+1:3*nRegs);  % est_fb
   t_delay = x(3*nRegs+1);
   t_delay=t_delay*ones(1,nRegs);
end
% want all fval to be zero except last region 
for ii=1:nRegs
   fvalues(ii)=0;
end
fvalues(nRegs)=fval(1);

if exist('seriesNumAIF')
   outfilename=mpi_getFilename(sliceNum, studyNum, outpath, numAzimuthalRegions, numRadialRegions,  clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp, seriesNumAIF, sliceNumAIF, scaleAIF, flagTimeStamps)
else
   outfilename=mpi_getFilename(sliceNum, studyNum, outpath, numAzimuthalRegions, numRadialRegions,  clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp)
end
outfilename=strcat(outfilename,'.globaldelayAndVp2.txt')

fid=fopen(outfilename,'w');
flow=kwi/0.5;
index=1:length(flow);
fprintf(fid,'Created %s ,  delta_t=%f\n',date, delta_t);
out=[index; flow]; 

fprintf(fid,'flow  %d    %6.3f\n',out); 
fprintf(fid,'\n');

fprintf(fid,'flowMean and Std  %6.2f %6.2f , coeff var= %6.2f\n',mean(flow), std(flow), std(flow)/mean(flow));
fprintf(fid,'\n');
Hct=0.45; 
out=[index; (1-Hct)*kwi./kwo];   % becuase kwi for arterial
fprintf(fid,'ve    %d     %6.3f\n',out);
fprintf(fid,'\n');
out=[index; est_fb]; 
fprintf(fid,'est_fb    %d     %6.3f\n',out);
fprintf(fid,'\n');
out=[index; t_delay]; 
fprintf(fid,'t0    %d     %6.3f\n',out);
fprintf(fid,'\n');
out=[index; fvalues]; 
fprintf(fid,'fval  %d     %6.3f\n',out);
fclose(fid);

return;
