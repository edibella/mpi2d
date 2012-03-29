function fitMultipleRegions(sliceNum,studyNum,outpath,filenamePart1,delta_t, numAzimuthalRegions, numRadialRegions, flagPixelwise, numSkip, curvesFromClusters)

if nargin==0,
    error('sliceNum argument needed for fitMultipleRegions');
end
if ~exist('numSkip'),   % if manually skipping an intital bump desired
   numSkip=0;
end

if ~exist('curvesFromClusters')
   curvefilename=strcat(outpath,filenamePart1,'.study',int2str(studyNum),'.slice',int2str(sliceNum),'.mat')
   load(curvefilename);
   curves = eval(filenamePart1);   % may want to re-do this so not so tricky
else
   curves=curvesFromClusters;
end

bldcurve=curves(1,:)';
%gd_curves =[ bldGds'; tissGds'];  % so in same format as curves

nRegs=size(curves,1)-1
tisscurves=curves(2:nRegs+1,:);
tisscurve=tisscurves';
nTimes=length(bldcurve);

% measured = fv*bld + (1-fv)*( conv( bld , k1*exp(-k2*(t-t0))) )

%options=optimset('TolX',1)             
%options=optimset('fminsearch')             
oldoptions=optimset('fmincon')             
options=optimset('fmincon')             
%options=optimset(oldoptions,'TolFun',0.01,'TolX',0.01,'MaxIter',1e3) ;
options=optimset(oldoptions,'TolFun',0.001,'TolX',0.001,'MaxIter',1e3) ;
options=optimset(oldoptions,'TolFun',0.0001,'TolX',0.0001,'MaxIter',1e4) ;
%options=optimset(oldoptions,'GradObj','on')             
%options=optimset(oldoptions,'MaxFunEvals', 1e6,'MaxIter',1e4)             
%options=optimset(oldoptions,'LevenbergMarquardt','on','GradObj','on','MaxFunEvals', 1e6,'MaxIter',1e4,'TolFun',0.005)             
%options=optimset(oldoptions,'LevenbergMarquardt','on','MaxFunEvals', 1e6,'MaxIter',1e4)             



%X0=0.1*ones(1,4*nRegs);
%nFitRuns=1;
%for ii=1:nFitRuns
% revamp this 11/16/03
%   X0(1)=0.5; LB(1)=0.01*X0(1); UB(1)=20*X0(1); 
%   X0(2)=0.5; LB(2)=0.01*X0(2); UB(2)=20*X0(2); 
%   X0(4)=-0.4; LB(4)=-15.0; UB(4)=0.8;  % changed 8/02/05 evrd

for ii=1:nRegs
   X0(ii)=0.5; LB(ii)=0.01*X0(ii); UB(ii)=20*X0(ii); 
   X0(nRegs+ii)=0.5; LB(nRegs+ii)=0.01*X0(nRegs+ii); UB(nRegs+ii)=20*X0(nRegs+ii); 
%   X0(2*nRegs+ii)=0.0; LB(2*nRegs+ii)=0.0;  UB(2*nRegs+ii)=0.6; 
%   X0(3*nRegs+ii)=-0.4; LB(3*nRegs+ii)=-15.0; UB(3*nRegs+ii)=0.8;  % changed 8/02/05 evrd
end
globalVp=1
if globalVp
   X0(2*nRegs+1)=0.0; LB(2*nRegs+1)=0.0; UB(2*nRegs+1)=0.4; 
   X0(2*nRegs+2)=-0.4; LB(2*nRegs+2)=-15.0; UB(2*nRegs+2)=0.8; 
else
   X0(3)=0.0; LB(3)=0.0;  UB(3)=0.6;
   X0(3*nRegs+1)=-0.4; LB(3*nRegs+1)=-15.0; UB(3*nRegs+1)=0.8;
end

whos X0

%   LB=0.001*X0; UB=200*X0;
%   LB=0.1*X0; UB=20*X0;
%   LB(4)=-10;

%UB(2*nRegs+1:3*nRegs)=0.7;
A=zeros(1,length(X0));
B=0;
A=[];
B=[];
nonlcon=[];

   disp('Fitting - doing iterations on all regions at once')
%for ii=1:nRegs
% 8/10/05

   [x,fval,exitflag,output,lambda,grad,hessian]=fmincon('funcMultipleRegions_fit',X0,A,B,A,B,LB,UB,nonlcon,options,bldcurve(numSkip+1:nTimes),tisscurve,delta_t, nRegs);  % last 3 are params passed to func(x,p1,p2,p3)

%  covMatrix=inv(hessian);
%  varKwi(ii)=covMatrix(1,1);
%   ii
%   save ftmp ii

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

clusterFlag=0; flagPixelwise=0; useIntegralLinearFit=0;
fixedDelay=99; fixedVp=99;
if exist('seriesNumAIF')
   outfilename=mpi_getFilename(sliceNum, studyNum, outpath, numAzimuthalRegions, numRadialRegions,  clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp, seriesNumAIF, sliceNumAIF, scaleAIF)
else
   outfilename=mpi_getFilename(sliceNum, studyNum, outpath, numAzimuthalRegions, numRadialRegions,  clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp)
end
outfilename=strcat(outfilename,'.globaldelayAndVp.txt')

fid=fopen(outfilename,'w');
flow=kwi/0.5;
index=1:length(flow);
fprintf(fid,'Created %s ,  delta_t=%f\n',date, delta_t);
out=[index; flow];

fprintf(fid,'flow  %d    %6.3f\n',out);
fprintf(fid,'\n');

fprintf(fid,'flowMean and Std  %6.2f %6.2f , coeff var= %6.2f\n',mean(flow), std(flow),std(flow)/mean(flow));
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
