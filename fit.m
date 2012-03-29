function fitparams = fit(sliceNum,studyNum,outpath,filenamePart1,delta_t, numAzimuthalRegions, numRadialRegions, flagPixelwise, numSkip, curvesFromClusters, fixedVp, fixedDelay)
% need some wort of vararg!!!!
%flagPixelwise=0;

if nargin==0,
    error('sliceNum argument needed for mpi_fitCurves');
end
if ~exist('numSkip'),   % if manually skipping an intital bump desired
   numSkip=0;
end


useIntegralLinearFit=0;
clusterFlag=0;
%if ~exist('curvesFromClusters')
   curvefilename=strcat(outpath,filenamePart1,'.study',int2str(studyNum),'.slice',int2str(sliceNum),'.mat')
   load(curvefilename);
   curves = eval(filenamePart1);   % may want to re-do this so not so tricky
%else
%   curves=curvesFromClusters;
%end

bldcurve=curves(1,:)';
%gd_curves =[ bldGds'; tissGds'];  % so in same format as curves

nRegs=size(curves,1)-1
tisscurves=curves(2:nRegs+1,:);
tisscurve=tisscurves';
nTimes=length(bldcurve);

% measured = fv*bld + (1-fv)*( conv( bld , k1*exp(-k2*(t-t0))) )

%options=optimset('TolX',1)             
%options=optimset('fminsearch')             
oldoptions=optimset('fmincon');             
options=optimset('fmincon')   ;          
%options=optimset(oldoptions,'TolFun',0.01,'TolX',0.01,'MaxIter',1e3) ;
options=optimset(oldoptions,'TolFun',0.001,'TolX',0.001,'MaxIter',1e3) ;
options=optimset(oldoptions,'TolFun',0.0001,'TolX',0.0001,'MaxIter',1e4) ;
%options=optimset(oldoptions,'GradObj','on')             
%options=optimset(oldoptions,'MaxFunEvals', 1e6,'MaxIter',1e4)             
%options=optimset(oldoptions,'LevenbergMarquardt','on','GradObj','on','MaxFunEvals', 1e6,'MaxIter',1e4,'TolFun',0.005)             
%options=optimset(oldoptions,'LevenbergMarquardt','on','MaxFunEvals', 1e6,'MaxIter',1e4)             

 options=optimset(options,'Display', 'off');
warning('off', 'all');

%X0=0.1*ones(1,4*nRegs);
%nFitRuns=1;
%for ii=1:nFitRuns
%   X0=ii*0.2*ones(1,4*nTissRegs);
% revamp this 11/16/03
   X0(1)=0.5; LB(1)=0.01*X0(1); UB(1)=20*X0(1); 
   X0(2)=0.5; LB(2)=0.01*X0(2); UB(2)=20*X0(2); 

% try this 9/15/05
   X0(1)=1.5; LB(1)=0.001*X0(1); UB(1)=90*X0(1); 
   X0(2)=0.5; LB(2)=0.001*X0(2); UB(2)=90*X0(2); 
%   X0(3)=0.22222;  
%   X0(3)=0.1; LB(3)=0.0;  UB(3)=0.6;  

%   X0(3)=0.0; LB(3)=0.0;  UB(3)=0.6;   % changed 9/13/04 to try

%   X0(4)=0.1
%   X0(4)=12; LB(4)=11; UB(4)=12;  % positive t0 moves fit to right
%   X0(4)=0; LB(4)=-.125; UB(4)=6;

%   X0(4)=-0.4; LB(4)=-15.0; UB(4)=0.8;  % changed 8/02/05 evrd

if fixedDelay==99
   X0(3)=-0.4; LB(3)=-15.0; UB(3)=0.8;  % changed 8/25/05 evrd to do t_delay thrid
else
  if fixedVp==99
      X0(3)=0.0; LB(3)=0.0;  UB(3)=0.7;
  end
end
if fixedVp==99 & fixedDelay==99
   X0(4)=0.0; LB(4)=0.0;  UB(4)=0.7;   % changed 9/13/04 to try
end

% change for Fermi function 9/21/05
   X0(4)=0.0; LB(4)=-5.0;  UB(4)=9.0;  

%UB(2*nRegs+1:3*nRegs)=0.7;
A=zeros(1,length(X0));
B=0;
A=[];
B=[];
nonlcon=[];


for ii=1:nRegs
   %[x,options]=fminsearch('funcCM',X0);
   %[x,options]=fmincon('func_nonPCA',X0,A,B,A,B,LB,UB,'dummyfunc',options);
 %  disp('Fitting - doing iterations')
   %[x,fval,exitflag,output]=fmincon('func_fit',X0,A,B,A,B,LB,UB,nonlcon,options,curves(:,numSkip+1:nTimes),delta_t, nRegs);  % last 3 are params passed to func(x,p1,p2,p3)
%  if ii>=600 
%  if ii~=996 && ii~=997 && ii~=1093 
%if ii~=142 && ii~=143
%   [x,fval,exitflag,output]=fmincon('func_fit',X0,A,B,A,B,LB,UB,nonlcon,options,bldcurve(numSkip+1:nTimes),tisscurve(:,ii),delta_t, 1);  % last 3 are params passed to func(x,p1,p2,p3)
% 8/10/05
%  [x,fval,exitflag,output,lambda,grad,hessian]=fmincon('func_fit',X0,A,B,A,B,LB,UB,nonlcon,options,bldcurve(numSkip+1:nTimes),tisscurve(:,ii),delta_t, 1);  % last 3 are params passed to func(x,p1,p2,p3)
 [x,fval,exitflag,output,lambda,grad,hessian]=fmincon('func_fit',X0,A,B,A,B,LB,UB,nonlcon,options,bldcurve(numSkip+1:nTimes),tisscurve(:,ii),delta_t, 1, fixedDelay, fixedVp);  % last 3 are params passed to func(x,p1,p2,p3)
% [x,fval,exitflag,output,lambda,grad,hessian]=fmincon('fermi_func_fit',X0,A,B,A,B,LB,UB,nonlcon,options,bldcurve(numSkip+1:nTimes),tisscurve(:,ii),delta_t, 1, fixedDelay, fixedVp);  % last 3 are params passed to func(x,p1,p2,p3)
%  end
%end
   %[x,fval]=fminsearch('func_nonPCA',X0);  % simplex
   %[x,options]=constr('funcCM',X0);

  covMatrix=inv(hessian);
  varKwi(ii)=covMatrix(1,1);
   ii
%   save ftmp ii
   kwi(ii) = x(1);
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

%est_fb(ii)=0.05;
%est_fb(ii)=0.00;
%t_delay(ii) = -3.0;
%t_delay(ii) = -0.5;
%   est_fb(ii)= 0;
%   t_delay(ii) = 0;
   fvalues(ii)=fval;

[est_curves, tmpbldcurve]=mpi_fwdModel(bldcurve, delta_t, x(1), x(2), t_delay(ii), est_fb(ii));

end

warning('on', 'all');

fitparams(1:nRegs) = kwi;
%if alternateAIF_flag==1
%   kwi=kwi/scalefactor;
%end

fitparams(nRegs+1:2*nRegs) = kwo ;
fitparams(2*nRegs+1:3*nRegs) = est_fb;
fitparams(3*nRegs+1:4*nRegs) = t_delay;
fitparams(4*nRegs+1:5*nRegs) = fvalues;

% write out in ascii
%if exist('curvesFromClusters')
%if curvesFromClusters~=0
%   outfilename=strcat(outpath,'fitparams.study',int2str(studyNum),'.slice',int2str(sliceNum),'.clusters','.txt');
%elseif (flagPixelwise==1)
%   outfilename=strcat(outpath,'fitparams.study',int2str(studyNum),'.slice',int2str(sliceNum),'.pixelwise','.txt');
%else

if exist('seriesNumAIF')
   outfilename=mpi_getFilename(sliceNum, studyNum, outpath, numAzimuthalRegions, numRadialRegions,  clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp, seriesNumAIF, sliceNumAIF, scaleAIF)
else
   outfilename=mpi_getFilename(sliceNum, studyNum, outpath, numAzimuthalRegions, numRadialRegions,  clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp)
end

%if fixedDelay==99 & fixedVp==99
%   outfilename=strcat(outpath,'flowvalues.study',int2str(studyNum),'.slice',int2str(sliceNum),'.',int2str(numAzimuthalRegions), '.',int2str(numRadialRegions),'.txt');
%elseif fixedDelay==99 & fixedVp~=99 
%   outfilename=strcat(outpath,'flowvalues.study',int2str(studyNum),'.slice',int2str(sliceNum),'.',int2str(numAzimuthalRegions), '.',int2str(numRadialRegions),'_fixedVp',num2str(fixedVp),'.txt');
%elseif fixedDelay~=99 & fixedVp==99
%   outfilename=strcat(outpath,'flowvalues.study',int2str(studyNum),'.slice',int2str(sliceNum),'.',int2str(numAzimuthalRegions), '.',int2str(numRadialRegions),'_fixedDelay',num2str(fixedDelay),'.txt');
%elseif fixedDelay~=99 & fixedVp~=99
%   outfilename=strcat(outpath,'flowvalues.study',int2str(studyNum),'.slice',int2str(sliceNum),'.',int2str(numAzimuthalRegions), '.',int2str(numRadialRegions),'_fixedDelay',num2str(fixedDelay),'_fixedVp',num2str(fixedVp),'.txt');
%else
%  disp('error in fit writing  out parameters')
%end

   outfilename1=strcat(outpath,'flowvalues.study',int2str(studyNum),'.slice',int2str(sliceNum),'.',int2str(numAzimuthalRegions), '.',int2str(numRadialRegions),'.txt');
   outfilename2=strcat(outpath,'othervalues.study',int2str(studyNum),'.slice',int2str(sliceNum),'.',int2str(numAzimuthalRegions), '.',int2str(numRadialRegions),'.txt');


fid=fopen(outfilename,'w');
%fid1=fopen(outfilename1,'w');
%fid2=fopen(outfilename2,'w');
index=1:length(kwi);
%fprintf(fid,'Created %s ,  delta_t=%f   Order below is kwi, kwo, fv, t_delay, fval\n',date, delta_t);
fprintf(fid,'Created %s ,  delta_t=%f\n',date, delta_t);
out=[index; kwi./0.5]; 

fprintf(fid,'flow  %d    %6.3f\n',out);  % assume E=0.5, kwi=E*F
fprintf(fid,'\n');
out=[index; sqrt(varKwi)]; 
%fprintf(fid,'flowStd  %d    %6.3f\n',out);  % assume E=0.5, kwi=E*F
out=[index; sqrt(varKwi)]; 
fprintf(fid,'flowMean and Std  %6.2f %6.2f \n',mean(kwi./0.5), std(kwi./0.5));  % assume E=0.5, kwi=E*F
fprintf(fid,'\n');

out1=kwi./0.5;
%fprintf(fid1,'%6.3f ',out1);
%fclose(fid1);



Hct=0.45; 
out=[index; (1-Hct)*kwi./kwo];   % becuase kwi for arterial
fprintf(fid,'ve    %d     %6.3f\n',out);
fprintf(fid,'\n');

out1=(1-Hct)*kwi./kwo;
%fprintf(fid2,'%6.3f ',out1);


out=[index; est_fb]; 
fprintf(fid,'fv    %d     %6.3f\n',out);
fprintf(fid,'\n');

out1=est_fb;
%fprintf(fid2,'%6.3f ',out1);

out=[index; t_delay]; 
fprintf(fid,'t0    %d     %6.3f\n',out);
fprintf(fid,'\n');

out1=t_delay;
%fprintf(fid2,'%6.3f ',out1);
%fclose(fid2);



out=[index; fvalues]; 
fprintf(fid,'fval  %d     %6.3f\n',out);

%for ireg=1:nRegs
%   fprintf(fid,'%d		%6.3f\n',ireg,kwi(ireg));
%end
%fprintf(fid,'\n');
%for ireg=1:nRegs
%%   fprintf(fid,'%d		%6.3f     kwo\n',ireg,kwo(ireg));
%   fprintf(fid,'%d		%6.3f\n',ireg,kwo(ireg));
%end
%fprintf(fid,'\n');

fclose(fid);


%showfits(curves, fitparams, delta_t);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% this is in a separate file. Commenting out 10/31/02 so don't get confused. (but can still get back if needed)
%function showfits(curves, studyNum, sliceNum)
%
%bldcurve=curves(1,:)';
%
%nRegs=8;
%tisscurves=curves(2:nRegs+1,:);
%tisscurve=tisscurves';
%nTimes=length(bldcurve);
%
%% read 1st line
%filename=strcat('fitparams.study',int2str(studyNum),'.slice',int2str(sliceNum),'.txt');
%fid=fopen(filename ,'r' );
%a =fscanf(fid,'%s',3);
%a =fscanf(fid,'%s',1);
%delta_t=str2num(a(9:14));
%
%% Start processing.
%counter=0;
%while(counter <= rowstart)
%    line = fgetl(fid);
%    counter=counter+1;
%end
%
%while ~isempty(a)
%  a =fscanf(fid,'%s',1);
%  switch(a)
%     case {'kwi'}
%      	disp('kwi')
%    	a =fscanf(fid,'%s',1);
%        aa=str2num(a);
%        kwi(i)=aa;
%     case {'kwo'}
%     	disp('kwo')
%    	a =fscanf(fid,'%s',1);
%        aa=str2num(a);
%        kwo(i)=aa;
%     case {'fv'}
%    	a =fscanf(fid,'%s',1);
%        aa=str2num(a);
%        est_fb(i)=aa;
%     case {'t0'}
%    	a =fscanf(fid,'%s',1);
%        aa=str2num(a);
%        t_delay(i)=aa;
%     otherwise
% 	disp('not listed')
%  end
%
%end
%
%%[hh,stuff]=myhdrload(infile,7)
%%if ~exist('ntimes'),
%%   ntimes=length(stuff);
%%end
%%time=stuff(:,1);
%
%
%tmp=0:nTimes-1;
%ttimes = delta_t * tmp';
%ttimes = tmp';
%
%kwi = fitparams(1:nRegs);
%kwo = fitparams(nRegs+1:2*nRegs);
%est_fb = fitparams(2*nRegs+1:3*nRegs);
%t_delay = fitparams(3*nRegs+1:4*nRegs);
%
%imp=zeros(nTimes,1);
%imp(1)=1.0;
%
%errsum = 0.0;
%clf
%hold on
%for i=1:nRegs
%disp('iReg is '),i
%      est_curves = conv(bldcurve, delta_t/60*((1-est_fb(i))*kwi(i)*exp(-kwo(i)*delta_t/60*(ttimes-t_delay(i)*ones(nTimes,1))) + est_fb(i)*imp) );
%      err = (tisscurve(:,i) - est_curves(1:nTimes));
%      errsum = errsum + sum(err.*err)
%   clf; hold on
%   set(gca,'FontSize',16); set(gcf,'Color',[1 1 1]);
%   plot(ttimes,bldcurve,'r','LineWidth',3)
%   plot(ttimes,tisscurve(:,i),'x','LineWidth',2.2)
%   plot(ttimes,est_curves(1:nTimes),'k','LineWidth',3)
%   legend('Blood Input','Measured Uptake', 'Fit (Model)')
%   plot(ttimes,est_curves(1:nTimes),':')
%   xlabel('Time')
%   ylabel('Gd concentration (mmol)')
%   pause
%end
%
%return
