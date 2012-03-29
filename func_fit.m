function [Fsum]=func_fit(x,bldcurve, tisscurve,delta_t,nRegs,fixedDelay, fixedVp)
%  matlab function to be called by fmincon


%disp('iteration ')
%ii=ii+1


%bldcurve=curves(1,:)';
%tisscurves=curves(2:nRegs+1,:);
%tisscurve=tisscurve';
nTimes=length(bldcurve);


tmp=0:nTimes-1;
%ttimes = delta_t * tmp';
ttimes = tmp';

kwi = x(1:nRegs);
kwo = x(nRegs+1:2*nRegs);

t_delay = fixedDelay*ones(1,nRegs);
est_fb = fixedVp*ones(1,nRegs);
if (length(x)>2*nRegs)
   t_delay = fixedDelay*ones(1,nRegs);
   if fixedDelay==99
       t_delay = x(2*nRegs+1:3*nRegs);
   else
       est_fb = x(2*nRegs+1:3*nRegs);
   end
end
if (length(x)>3*nRegs)
   est_fb = x(3*nRegs+1:4*nRegs);
end

%t_delay = x(3*nRegs+1:4*nRegs);
%t_delay = 1.0*ones(1,nRegs);
%t_delay = 0.0*ones(1,nRegs);
%t_delay = -0.5*ones(1,nRegs);
%est_fb = 0.05*ones(1,nRegs);
%est_fb = 0.0*ones(1,nRegs);
%est_fb = zeros(1,nRegs);  

jitter=0; newRate=1;    % for downsampling ISMRM 2003 (not used here)

imp=zeros(nTimes,1);
imp(1)=1.0;


errsum = 0.0;
for i=1:nRegs
   [est_curves, tmpbldcurve]=mpi_fwdModel(bldcurve, delta_t, kwi(i), kwo(i), t_delay(i), est_fb(i));
%      tmpbldcurve = interp1(1:nTimes,bldcurve,1-t_delay(i):nTimes-t_delay(i),'linear',0);
% try this 4/1/04
%      tmpbldcurve = interp1(1:nTimes,bldcurve,1-t_delay(i):nTimes-t_delay(i),'v5cubic',0);
%%%      tmpbldcurve = interp1(1:nTimes,bldcurve,1-t_delay(i):nTimes-t_delay(i),'splines',0);
%      est_curves = conv(tmpbldcurve, delta_t/60*((1-est_fb(i))*kwi(i)*exp(-kwo(i)*delta_t/60*(ttimes-t_delay(i)*ones(nTimes,1))) + est_fb(i)*imp) );
%%%      est_curves = conv(tmpbldcurve, ((1-est_fb(i))*(delta_t/60)*kwi(i)*exp(-kwo(i)*delta_t/60*(ttimes)) + est_fb(i)*imp) );
      est_curves=est_curves((jitter+1):newRate:nTimes);   % downsampling step
      err = (tisscurve(1:nTimes,i) - est_curves);
%      est_curves = conv(bldcurve, delta_t/60*((1-est_fb(i))*kwi(i)*exp(-kwo(i)*delta_t/60*(ttimes-t_delay(i)*ones(nTimes,1))) + est_fb(i)*imp) );
%      err = (tisscurve(:,i) - est_curves(1:nTimes));
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

end

Fsum = errsum;
%disp('ERROR is '),Fsum
%t_delay
%clf; plot(est_curves,'k')
%hold on
%plot(tisscurve(:,1),'x')



% constraint:   G < zeros(G)
% constrant a function of x (tissue values) to have norm of one? 
% try just all curves greater than 0 
% or say greater than 1

% not currently using
%G=-10;
% and constrain  fb
%G = -10*est_fb;

