



load PETcurvesOSEM_stress   %bld1 tiss1 tiss2

% these are from rest?  stress?  FBP?  
% see PET



nRegs=2
nTimes=length(bld1);

fixedDelay=99
fixedVp=99
delta_t=5   % seconds


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
 
   extrapval=0;
   origSampleTimes=[5*ones(1,14) 10*ones(1,6) 20*ones(1,3) 30*ones(1,4)];
   for ii=1:length(origSampleTimes)
       sumTimes(ii)=sum(origSampleTimes(1:ii))
   end
   
totalTime=14*5 + 6*10 + 3* 20 + 4*30
%from email, seems 14x5, 6x10, 3x20, 4x30
nTimes=round(totalTime/delta_t);

   newbld = interp1(sumTimes, bld1, delta_t*(0:nTimes-1),'cubic',extrapval);
% NOTE: interp1 seems different in matlab6! Using matlab7 here!
   newbld=newbld';



   newtiss1 = interp1(sumTimes, tiss1, delta_t*(0:nTimes-1),'cubic',extrapval);
% NOTE: interp1 seems different in matlab6! Using matlab7 here!
   newtiss1=newtiss1';
newtiss2 = interp1(sumTimes, tiss2, delta_t*(0:nTimes-1),'cubic',extrapval);
% NOTE: interp1 seems different in matlab6! Using matlab7 here!
   newtiss2=newtiss2';

   
% % add sanity check to make sure AIF and tissue curves are same length, since they
% % could come from different timeStamp files
%    if (length(bldcurve) > size(tisscurve,1))
%        bldcurve=bldcurve(1:size(tisscurve,1));  % tr
%        disp('Bld and tiss curves lenghts do not match due to timestamps...\n')
%        disp('Truncating bldcurve length to match tisscurve')
%        nTimes=size(tisscurve,1);
%    end
%    if (length(bldcurve) < size(tisscurve,1))
%        tisscurve=tisscurve(1:nTimes,:);
%        disp('Bld and tiss curves lenghts do not match due to timestamps...\n')
%        disp('Truncating tisscurve length to match bldcurve')
%    end
%  



%[x,fval,exitflag,output,lambda,grad,hessian]=fmincon('sat_func_fitNoBlood_orig',X0,A,B,A,B,LB,UB,nonlcon,options,bldcurve(numSkip+1:nTimes),sat_bldcurve(numSkip+1:nTimes),mean(tisscurve,2),delta_t, 1, fixedDelay, fixedVp);  
estGlobalDelay=0  %x(3)


%keyboard

modelType='full'
fvalues=[]; clear xx;
for ii=1:1
   disp('Fitting - doing iterations for noblood1 model')
   [x,fval,exitflag,output,lambda,grad,hessian]=fmincon('func_fit_PETwater',X0,A,B,A,B,LB,UB,nonlcon,options,newbld,newbld,newtiss1,delta_t, 1, fixedDelay, fixedVp, estGlobalDelay);  
   xx1(ii,:)= x;   % concat regions
   fvalues1=[fvalues fval];
end
%printFitParams(xx,fvalues, sliceNum, studyNum, outpath, delta_t, numAzimuthalRegions, numRadialRegions,  clusterFlag, flagPixelwise, useIntegralLinearFit, fixedDelay, fixedVp, seriesNumAIF, sliceNumAIF, scaleAIF,flagTimeStamps, modelType)


[x,fval,exitflag,output,lambda,grad,hessian]=fmincon('func_fit_PETwater',X0,A,B,A,B,LB,UB,nonlcon,options,newbld,newbld,newtiss2,delta_t, 1, fixedDelay, fixedVp, estGlobalDelay);  
   xx2(ii,:)= x;   % concat regions
   fvalues2=[fvalues fval];
   

