function [Fsum]=funcMultipleRegions_fit(x,bldcurve, tisscurve,delta_t,nRegs)
% matlab function to be called by fmincon

%disp('iteration ')
%ii=ii+1

if ~exist('modelType')
    modelType='xxx';
end

%bldcurve=curves(1,:)';
%tisscurves=curves(2:nRegs+1,:);
%tisscurve=tisscurve';
nTimes=length(bldcurve);

%tisscurve=tisscurve(1:nTimes);

tmp=0:nTimes-1;
%ttimes = delta_t * tmp';
ttimes = tmp';

kwi = x(1:nRegs);
kwo = x(nRegs+1:2*nRegs);

globalVp=0;
if length(x) == 2*nRegs+2
   globalVp=1;
   est_fb = x(2*nRegs+1);
   t_delay = x(2*nRegs+2);
else
   est_fb = x(2*nRegs+1:3*nRegs);
   t_delay = x(3*nRegs+1);
   est_fb = zeros(1,nRegs); 
end
%t_delay = 1.0*ones(1,nRegs);
%t_delay = 0.0*ones(1,nRegs);
%est_fb = zeros(1,nRegs);  
%est_fb = 0.1*ones(1,M);

jitter=0; newRate=1;    % for downsampling ISMRM 2003 (not used here)

imp=zeros(nTimes,1);
imp(1)=1.0;

errsum = 0.0;


if strcmp(modelType,'fermi')
   for i=1:nRegs
      [est_curves, tmpbldcurve]=mpi_fermi_fwdModel(bldcurve, delta_t, kwi(i), kwo(i), t_delay, est_fb(i));
      est_curves=est_curves((jitter+1):newRate:nTimes);   % downsampling step
      err = (tisscurve(:,i) - est_curves);
      errsum = errsum + sum(err.*err);
   end
   Fsum = errsum;
else
   for i=1:nRegs
      if globalVp 
         [est_curves, tmpbldcurve]=mpi_fwdModel(bldcurve, delta_t, kwi(i), kwo(i), t_delay, est_fb);
      else
         [est_curves, tmpbldcurve]=mpi_fwdModel(bldcurve, delta_t, kwi(i), kwo(i), t_delay, est_fb(i));
      end
      est_curves=est_curves((jitter+1):newRate:nTimes);   % downsampling step
      err = (tisscurve(:,i) - est_curves);
      errsum = errsum + sum(err.*err);
%clf 
%hold on
%plot(bldcurve,'r')
%plot(tmpbldcurve,'m')
%plot(est_curves(1:nTimes))
%%plot(tisscurve(:,i),'g')
%plot(tisscurve(:,i),'kx')
%drawnow
%t_delay
%pause
   end
   Fsum = errsum;
end

return

